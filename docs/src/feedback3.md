# Feedback 3



```@example sir3
using DifferentialEquations
using LinearAlgebra
using Random
using Plots
using LaTeXStrings
using Colors

default(fontfamily = "Computer Modern", framestyle = :box)
```

```@example sir3
β       = 0.8
γ       = 0.2
S0      = 0.92
I0      = 0.04
Qbudget = 16.0
tf      = 100.0
nothing
```

The herd-immunity threshold is

```math
S_h=\frac{\gamma}{\beta}.
```

```@example sir3
Sh = γ / β
```

### Target levels

As in Feedback 1 and Feedback 2, the reference trajectory uses the
optimal target infection level computed from the true initial conditions.
For estimated initial conditions $(\hat\theta_1,\hat\theta_2)$, define

```math
I_h
=
\hat\theta_2+\hat\theta_1-S_h
-S_h\log\left(\frac{\hat\theta_1}{S_h}\right).
```

The corresponding critical intervention budget is

```math
Q_{\mathrm{crit}}
=
\frac{I_h-\hat\theta_2}
{\beta S_h\hat\theta_2}.
```

The reference target value is therefore

```math
I^\star=
\begin{cases}
\dfrac{I_h}{1+\beta S_hQ},
&
Q<Q_{\mathrm{crit}},\\[2mm]
\hat\theta_2,
&
Q\geq Q_{\mathrm{crit}}.
\end{cases}
```

Feedback 3 instead uses a threshold that is recomputed from the current
estimated state and the remaining intervention budget. In the numerical
implementation, this quantity is

```math
I^\dagger(\hat I,\hat S,C)
=
\frac{
\hat I+\hat S-S_h-S_h\log\left(\dfrac{\hat S}{S_h}\right)
}{
1+\beta S_h C
}.
```

Thus, unlike Feedbacks 1 and 2, the switching threshold is not fixed by
the initial condition or by a constant lower bound: it evolves together
with $(\hat S,\hat I,C)$.

```@example sir3
function ISTAR(θ̂1::Real, θ̂2::Real)
    Ih = θ̂2 + θ̂1 - Sh - Sh * log(θ̂1 / Sh)

    Qcrit =
        (Ih - θ̂2) /
        (β * Sh * max(θ̂2, 1e-12))

    return Qbudget < Qcrit ?
           Ih / (1.0 + β * Sh * Qbudget) :
           θ̂2
end

Ibar = ISTAR(S0, I0)

safe_ = 1e-10

function IDAG(Ŝ::Real, Î::Real, C::Real)
    Ssafe = max(Ŝ, safe_)

    Ih =
        Î + Ssafe - Sh -
        Sh * log(Ssafe / Sh)

    return Ih /
           (1.0 + β * Sh * max(C, 0.0))
end
nothing
```

## Feedback laws

Two feedback laws are considered.

The first one is the adaptive feedback control based on the current
threshold $I^\dagger(\hat I,\hat S,C)$:

```math
u_{3}(\hat I,\hat S,C):=
\begin{cases}
1-\dfrac{S_h}{\hat S},
& \text{if } \hat I \geq I^\dagger(\hat I,\hat S,C)
\text{ and } \hat S>S_h,\\[6pt]
0, & \text{otherwise}.
\end{cases}
```

The modified feedback contains the same additional multiplicative
correction used for the other feedback strategies:

```math
\widetilde u_{3}(\hat I,\hat S,C):=
\begin{cases}
1-\dfrac{S_h}{\hat S}
\dfrac{I^\dagger(\hat I,\hat S,C)}{\hat I},
& \begin{aligned}
  \text{if } \hat I \geq I^\dagger(\hat I,\hat S,C) \\
  \text{and } \hat S > S_h,
\end{aligned} \\[8pt]
0, & \text{otherwise}.
\end{cases}
```

The ratio $I^\dagger/\hat I$ strengthens the intervention when the
estimated infected population lies above the adaptive threshold. Since
$I^\dagger$ depends on the remaining budget $C$, both the switching
surface and the modified control value evolve during the closed-loop
trajectory.

```@example sir3
@inline incidence(S, I, u) =
    β * S * I * (1.0 - u)

CONTROL_LAW = :u3
```

Possible values for `CONTROL_LAW`:
- `:u3`
- `:u3_tilde`

```@example sir3
u3_exact(Ŝ, Î, Idag) =
    (Ŝ ≤ Sh || Î < Idag) ?
    0.0 :
    1.0 - Sh / Ŝ

u3_tilde_exact(Ŝ, Î, Idag) =
    (Ŝ ≤ Sh || Î < Idag) ?
    0.0 :
    1.0 - (Sh / Ŝ) * (Idag / Î)

ft_control(Ŝ, Î, Idag) =
    CONTROL_LAW === :u3 ?
    u3_exact(Ŝ, Î, Idag) :
    u3_tilde_exact(Ŝ, Î, Idag)

nothing
```

### Smooth control used for the EKF

For the Extended Kalman Filter simulation, the discontinuous switching
conditions are regularized.

We introduce

```math
\sigma_\varepsilon(z)
=
\frac12
\left(
1+\tanh\left(\frac{z}{\varepsilon}\right)
\right).
```

```@example sir3
ϵcontrol = 1e-6

@inline function sigmoid_stable(z, ϵ)

    r = z / ϵ

    if r <= -40.0
        return 0.0
    elseif r >= 40.0
        return 1.0
    else
        return 0.5 * (1.0 + tanh(r))
    end
end
nothing
```

For $u_3$,

```@example sir3
function nsn_u3(Ŝ, Î, C)

    Ssafe = max(Ŝ, safe_)

    Idag =
        IDAG(Ssafe, Î, C)

    σI = sigmoid_stable(
        Î - Idag,
        ϵcontrol
    )

    σS = sigmoid_stable(
        Ssafe - Sh,
        ϵcontrol
    )

    return clamp(
        (1.0 - Sh / Ssafe) *
        σI *
        σS,
        0.0,
        1.0
    )
end
nothing
```

For $\widetilde u_3$,

```@example sir3
function nsn_u3_tilde(Ŝ, Î, C)

    Ssafe = max(Ŝ, safe_)
    Isafe = max(Î, safe_)

    Idag =
        IDAG(Ssafe, Î, C)

    σI = sigmoid_stable(
        Î - Idag,
        ϵcontrol
    )

    σS = sigmoid_stable(
        Ssafe - Sh,
        ϵcontrol
    )

    return clamp(
        (1.0 - Sh / Ssafe) *
        (Idag / Isafe) *
        σI *
        σS,
        0.0,
        1.0
    )
end

nsn_control(Ŝ, Î, C) =
    CONTROL_LAW === :u3 ?
    nsn_u3(Ŝ, Î, C) :
    nsn_u3_tilde(Ŝ, Î, C)

nothing
```

# Reference trajectory

The reference solution is divided into three phases.

### Phase 1 — Uncontrolled epidemic

Initially,

```math
u(t)=0.
```

The system therefore evolves according to

```math
\begin{aligned}
\dot S &= -\beta SI,\\
\dot I &= \beta SI-\gamma I.
\end{aligned}
```

until the infection reaches the prescribed level $I^\star$.

```@example sir3
tgrid = 0.0:0.05:tf

cb1 =
    ContinuousCallback(
        (v, t, integrator) ->
            v[2] - Ibar,
        terminate!
    )

sol1 = solve(
    ODEProblem(
        (dv, v, p, t) -> begin

            y = β * v[1] * v[2]

            dv[1] = -y
            dv[2] =  y - γ * v[2]
            dv[3] =  0.0

        end,

        [S0, I0, Qbudget],
        (0.0, tf)
    ),

    saveat = tgrid,
    callback = cb1,
    reltol = 1e-8,
    abstol = 1e-8
)
nothing
```

### Phase 2 — Intervention arc

Once $I=I^\star$, the feedback

```math
u(t)=1-\frac{S_h}{S(t)}
```

is applied.

The budget variable satisfies

```math
\dot C(t)=-u(t).
```

The intervention ends either when the budget is exhausted or when
$S=S_h$.

```@example sir3
cb2_C =
    ContinuousCallback(
        (v, t, integrator) -> v[3],
        terminate!
    )

cb2_S =
    ContinuousCallback(
        (v, t, integrator) ->
            v[1] - Sh,
        nothing,
        terminate!
    )

cb2 = CallbackSet(cb2_C, cb2_S)

sol2 = solve(
    ODEProblem(
        (dv, v, p, t) -> begin

            S, I, C =
                v[1], v[2], v[3]

            u = 1.0 - Sh / S

            y = β * S * I

            dv[1] = -y * (1 - u)
            dv[2] =  y * (1 - u) - γ * I
            dv[3] = -u

        end,

        sol1.u[end],
        (sol1.t[end], tf)
    ),

    saveat = tgrid,
    callback = cb2,
    reltol = 1e-8,
    abstol = 1e-8
)
nothing
```

### Phase 3 — Final uncontrolled trajectory

After the controlled phase,

```math
u(t)=0
```

again.

```@example sir3
sol3 = solve(
    ODEProblem(
        (dv, v, p, t) -> begin

            y = β * v[1] * v[2]

            dv[1] = -y
            dv[2] =  y - γ * v[2]
            dv[3] =  0.0

        end,

        sol2.u[end],
        (sol2.t[end], tf)
    ),

    saveat = tgrid,
    reltol = 1e-8,
    abstol = 1e-8
)
nothing
```

The three phases are concatenated as follows:

```@example sir3
t_ref = [
    sol1.t;
    sol2.t[2:end];
    sol3.t[2:end]
]

S_ref = [
    sol1[1, :];
    sol2[1, 2:end];
    sol3[1, 2:end]
]

I_ref = [
    sol1[2, :];
    sol2[2, 2:end];
    sol3[2, 2:end]
]

C_ref = [
    sol1[3, :];
    sol2[3, 2:end];
    sol3[3, 2:end]
]

n1 = length(sol1.t)
n2 = length(sol2.t) - 1
n3 = length(sol3.t) - 1

u_ref = [
    zeros(n1);
    1.0 .- Sh ./ S_ref[n1+1:n1+n2];
    zeros(n3)
]
nothing
```

# Extended Kalman Filter

We next compare the finite-time observer with a continuous-time
Extended Kalman Filter.

Let

```math
x=
\begin{pmatrix}
S\\
I
\end{pmatrix}.
```

The measured quantity is the epidemic incidence

```math
h(x,u)
=
\beta SI(1-u).
```

The EKF is written as

```math
\dot{\hat x}
=
f(\hat x,u)
+
PH^\top R^{-1}
\left(
y-h(\hat x,u)
\right),
```

with covariance equation

```math
\dot P
=
AP+PA^\top
-
PH^\top R^{-1}HP
+
Q.
```

### EKF parameters

```@example sir3
R =
    1e-4

Qproc =
    Matrix(
        Diagonal(
            [1e-2, 1e-2]
        )
    )

P0 =
    Matrix(
        Diagonal(
            [1e-2, 1e-2]
        )
    )

Ŝ0, Î0 =
    1.0, 0.0

meas_noise_amp =
    0.1
nothing
```

### Measurement noise

A deterministic multi-frequency signal is used to generate the
measurement perturbation.

```@example sir3
const NOISE_W =
    (
        0.17,
        0.41,
        0.93,
        2.13,
        3.71,
        7.29,
        11.93,
        19.07
    )

const NOISE_P =
    (
        0.00,
        1.13,
        2.47,
        0.62,
        4.01,
        5.28,
        1.97,
        3.55
    )

@inline function meas_noise(t)

    s = 0.0

    @inbounds for k in eachindex(NOISE_W)

        s += sin(
            NOISE_W[k] * t +
            NOISE_P[k]
        )

    end

    return meas_noise_amp *
           sqrt(2.0 / length(NOISE_W)) *
           s
end
nothing
```

The noisy measurement is

```math
y(t)
=
\beta S(t)I(t)
(1-u(t))
\left(
1+\nu(t)
\right).
```

### EKF dynamics

```@example sir3
function ekf_rhs!(dz, z, _, t)

    S =
        z[1]

    I =
        z[2]

    Ŝ =
        z[3]

    Î =
        z[4]

    C_hat =
        z[9]

    P =
        reshape(
            @view(z[5:8]),
            2,
            2
        )

    u_raw =
        nsn_control(Ŝ, Î, C_hat)

    u =
        C_hat > 0.0 ?
        u_raw :
        0.0

    h_true =
        incidence(S, I, u)

    dS =
        -h_true

    dI =
        h_true - γ * I

    y =
        h_true *
        (1.0 + meas_noise(t))

    S_eval =
        max(Ŝ, safe_)

    I_eval =
        max(Î, safe_)

    h_hat =
        incidence(
            S_eval,
            I_eval,
            u
        )

    innovation =
        y - h_hat

    A = [
        -β * I_eval * (1.0 - u)    -β * S_eval * (1.0 - u);
         β * I_eval * (1.0 - u)     β * S_eval * (1.0 - u) - γ
    ]

    H = [
        β * I_eval * (1.0 - u)     β * S_eval * (1.0 - u)
    ]

    K =
        (P * H') / R

    f_hat = [
        -β * S_eval * I_eval * (1.0 - u),

         β * S_eval * I_eval * (1.0 - u) -
         γ * I_eval
    ]

    d_hat =
        f_hat +
        vec(K) * innovation

    if Ŝ <= safe_ &&
       d_hat[1] < 0.0

        d_hat[1] = 0.0
    end

    if Î <= 0.0 &&
       d_hat[2] < 0.0

        d_hat[2] = 0.0
    end

    dP =
        A * P +
        P * A' -
        (P * H' * H * P) / R +
        Qproc

    dz[1] = dS
    dz[2] = dI

    dz[3] = d_hat[1]
    dz[4] = d_hat[2]

    dz[5:8] .= vec(dP)

    dz[9] = -u

    return nothing
end
nothing
```

The intervention is disabled after budget exhaustion:

```@example sir3
cb_budget_ekf =
    ContinuousCallback(
        (z, t, integrator) ->
            z[9],

        integrator ->
            (integrator.u[9] = 0.0)
    )
nothing
```

The complete EKF simulation is

```@example sir3
z0_ekf = [
    S0,
    I0,
    Ŝ0,
    Î0,
    vec(P0)...,
    Qbudget
]

prob_ekf =
    ODEProblem(
        ekf_rhs!,
        z0_ekf,
        (0.0, tf)
    )

sol_ekf =
    solve(
        prob_ekf,
        AutoTsit5(
            Rosenbrock23()
        );

        dt = 1e-4,
        dtmax = 5e-3,

        reltol = 1e-6,
        abstol = 1e-6,

        saveat = 0.05,

        maxiters = 10^6,

        callback =
            cb_budget_ekf
    )
nothing
```

The resulting trajectories are extracted using

```@example sir3
t_ekf =
    sol_ekf.t

S_ekf =
    [z[1] for z in sol_ekf.u]

I_ekf =
    [z[2] for z in sol_ekf.u]

Shat_ekf =
    [z[3] for z in sol_ekf.u]

Ihat_ekf =
    [z[4] for z in sol_ekf.u]

C_ekf =
    [z[9] for z in sol_ekf.u]

u_ekf =
    ifelse.(
        C_ekf .> 0.0,
        nsn_control.(
            Shat_ekf,
            Ihat_ekf,
            C_ekf
        ),
        0.0
    )

e_S_ekf =
    abs.(
        S_ekf .-
        Shat_ekf
    )

e_I_ekf =
    abs.(
        I_ekf .-
        Ihat_ekf
    )
nothing
```

# Finite-time observer

We now use a finite-time parameter estimator based on Dynamic
Regressor Extension and Mixing (DREM).

The transformed states are reconstructed according to

```math
\hat S
=
\hat\theta_1+\chi_1,
```

and

```math
\hat I
=
\chi_2+
e^{-\gamma t}
\hat\theta_2.
```

## FT Observer parameters

```@example sir3
λ_id = 1e-4

ρ_id = 1e6

α_id = 0.3

ε_id = 1e-6

Ŝ0_ft = 1.0

Î0_ft = 0.0

Y0 = β * S0 * I0
nothing
```

### DREM dynamics

Let

```math
\psi_1
=
(1-u)\beta\chi_2,
```

and

```math
\psi_2
=
(1-u)\beta e^{-\gamma t}\chi_1.
```

The filtered regression variables generate the determinant

```math
\Delta
=
\Omega_{11}\Omega_{22}
-
\Omega_{12}\Omega_{21}.
```

The finite-time parameter estimator has the structure

```math
\dot{\hat\theta}_i
=
\rho
\Delta
\frac{
Y_i-\Delta\hat\theta_i
}{
\left|
Y_i-\Delta\hat\theta_i
\right|^{1-\alpha}
+
\varepsilon
}.
```

The complete implementation is

```@example sir3
function sir_drem!(dv, v, p, t)

    egt =
        exp(-γ * t)

    Ŝ =
        v[11] +
        v[3]

    Î =
        v[4] +
        egt * v[12]

    Î_dag =
        IDAG(Ŝ, Î, v[13])

    u_raw =
        ft_control(
            Ŝ,
            Î,
            Î_dag
        )

    u =
        v[13] > 0.0 ?
        u_raw :
        0.0

    y =
        β *
        v[1] *
        v[2] *
        (1.0 - u)

    dv[1] =
        -y

    dv[2] =
        y -
        γ * v[2]

    dv[3] =
        -y

    dv[4] =
        -γ * v[4] +
        y

    ψ1 =
        (1.0 - u) *
        β *
        v[4]

    ψ2 =
        (1.0 - u) *
        β *
        egt *
        v[3]

    z =
        y -
        (1.0 - u) *
        (
            egt * Y0 +
            β * v[3] * v[4]
        )

    dv[5] =
        ψ1 * z -
        λ_id * v[5]

    dv[6] =
        ψ2 * z -
        λ_id * v[6]

    dv[7] =
        ψ1^2 -
        λ_id * v[7]

    dv[8] =
        ψ1 * ψ2 -
        λ_id * v[8]

    dv[9] =
        ψ2 * ψ1 -
        λ_id * v[9]

    dv[10] =
        ψ2^2 -
        λ_id * v[10]

    Δ =
        v[7] * v[10] -
        v[8] * v[9]

    Y1 =
        v[10] * v[5] -
        v[8] * v[6]

    Y2 =
        -v[9] * v[5] +
        v[7] * v[6]

    dv[11] =
        ρ_id *
        Δ *
        (Y1 - Δ * v[11]) /
        (
            abs(
                Y1 -
                Δ * v[11]
            )^(1 - α_id)
            +
            ε_id
        )

    dv[12] =
        ρ_id *
        Δ *
        (Y2 - Δ * v[12]) /
        (
            abs(
                Y2 -
                Δ * v[12]
            )^(1 - α_id)
            +
            ε_id
        )

    dv[13] =
        -u

    return nothing
end
nothing
```

### Initial condition and numerical integration

```@example sir3
v0_ft =
    zeros(13)

v0_ft[1] =
    S0

v0_ft[2] =
    I0

v0_ft[11] =
    Ŝ0_ft

v0_ft[12] =
    Î0_ft

v0_ft[13] =
    Qbudget
nothing
```

Budget exhaustion is handled using

```@example sir3
cb_budget_ft =
    ContinuousCallback(
        (v, t, integrator) ->
            v[13],

        integrator ->
            (integrator.u[13] = 0.0)
    )
nothing
```

The system is integrated with a stiff solver:

```@example sir3
prob_ft =
    ODEProblem(
        sir_drem!,
        v0_ft,
        (0.0, tf)
    )

sol_ft =
    solve(
        prob_ft,
        RadauIIA5(
            autodiff = false
        );

        reltol = 1e-9,
        abstol = 1e-11,

        dtmax = 0.1,

        saveat = 0.05,

        maxiters =
            Int(1e7),

        callback =
            cb_budget_ft
    )
nothing
```

### Estimated trajectories

```@example sir3
t_ft =
    sol_ft.t

S_ft =
    sol_ft[1, :]

I_ft =
    sol_ft[2, :]

C_ft =
    sol_ft[13, :]

χ1_ft =
    sol_ft[3, :]

χ2_ft =
    sol_ft[4, :]

θ1hat_ft =
    sol_ft[11, :]

θ2hat_ft =
    sol_ft[12, :]
nothing
```

The state estimates are reconstructed using

```@example sir3
egt_v_ft =
    exp.(
        -γ .* t_ft
    )

Shat_ft =
    θ1hat_ft .+
    χ1_ft

Ihat_ft =
    χ2_ft .+
    egt_v_ft .* θ2hat_ft

Idag_ft =
    IDAG.(
        Shat_ft,
        Ihat_ft,
        C_ft
    )

u_ft =
    ifelse.(
        C_ft .> 0.0,
        ft_control.(
            Shat_ft,
            Ihat_ft,
            Idag_ft
        ),
        0.0
    )
nothing
```

The estimation errors are

```@example sir3
e_S_ft =
    abs.(
        S_ft .-
        Shat_ft
    )

e_I_ft =
    abs.(
        I_ft .-
        Ihat_ft
    )
nothing
```

# Performance indicators

The maximum estimation errors are computed for both observers.

For the EKF,

```@example sir3
println(
    "EKF max |S-Ŝ| = ",
    maximum(e_S_ekf)
)

println(
    "EKF max |I-Î| = ",
    maximum(e_I_ekf)
)
```

For the finite-time observer,

```@example sir3
max_I_ft, idx_ft =
    findmax(I_ft)

println(
    "FT max I(t) = ",
    max_I_ft,
    " at t = ",
    t_ft[idx_ft]
)

println(
    "FT max |S-Ŝ| = ",
    maximum(e_S_ft)
)

println(
    "FT max |I-Î| = ",
    maximum(e_I_ft)
)
```

### Time spent above the infection threshold

An additional performance indicator measures how long the infection
trajectory remains above $I^\star$.

```@example sir3
function crossing_stats(
    tvec,
    Ivec,
    thresh
)

    n =
        length(tvec)

    entry =
        findfirst(
            k ->
                Ivec[k] <= thresh &&
                Ivec[k+1] > thresh,
            1:n-1
        )

    entry === nothing &&
        return (
            nothing,
            nothing,
            0.0
        )

    interp(k) =
        tvec[k] +
        (
            thresh -
            Ivec[k]
        ) *
        (
            tvec[k+1] -
            tvec[k]
        ) /
        (
            Ivec[k+1] -
            Ivec[k]
        )

    t_in =
        interp(entry)

    exitidx =
        findfirst(
            k ->
                k > entry &&
                Ivec[k] >= thresh &&
                Ivec[k+1] < thresh,
            1:n-1
        )

    t_out =
        exitidx === nothing ?
        tvec[end] :
        interp(exitidx)

    return (
        t_in,
        t_out,
        t_out - t_in
    )
end
nothing
```

Then,

```@example sir3
_, _, dt_ekf =
    crossing_stats(
        t_ekf,
        I_ekf,
        Ibar
    )

_, _, dt_ft =
    crossing_stats(
        t_ft,
        I_ft,
        Ibar
    )

println(
    "EKF: I(t) above Ibar for Δt ≈ ",
    dt_ekf
)

println(
    "FT : I(t) above Ibar for Δt ≈ ",
    dt_ft
)
```

# Closed loop with the modified law

Since `ft_control` and `nsn_control` dispatch on `CONTROL_LAW` at call
time, the same closed-loop simulations can be repeated with the
modified law simply by switching the symbol and re-solving the same
problems. No other code is modified.

```@example sir3
CONTROL_LAW = :u3_tilde
nothing
```

The EKF-based closed loop with $\widetilde u_3$:

```@example sir3
sol_ekf_ut =
    solve(
        prob_ekf,
        AutoTsit5(
            Rosenbrock23()
        );

        dt = 1e-4,
        dtmax = 5e-3,

        reltol = 1e-6,
        abstol = 1e-6,

        saveat = 0.05,

        maxiters = 10^6,

        callback =
            cb_budget_ekf
    )
nothing
```

```@example sir3
t_ekf_ut =
    sol_ekf_ut.t

S_ekf_ut =
    [z[1] for z in sol_ekf_ut.u]

I_ekf_ut =
    [z[2] for z in sol_ekf_ut.u]

Shat_ekf_ut =
    [z[3] for z in sol_ekf_ut.u]

Ihat_ekf_ut =
    [z[4] for z in sol_ekf_ut.u]

C_ekf_ut =
    [z[9] for z in sol_ekf_ut.u]

u_ekf_ut =
    ifelse.(
        C_ekf_ut .> 0.0,
        nsn_control.(
            Shat_ekf_ut,
            Ihat_ekf_ut,
            C_ekf_ut
        ),
        0.0
    )

e_S_ekf_ut =
    abs.(
        S_ekf_ut .-
        Shat_ekf_ut
    )

e_I_ekf_ut =
    abs.(
        I_ekf_ut .-
        Ihat_ekf_ut
    )
nothing
```

The finite-time-observer-based closed loop with $\widetilde u_3$:

```@example sir3
sol_ft_ut =
    solve(
        prob_ft,
        RadauIIA5(
            autodiff = false
        );

        reltol = 1e-9,
        abstol = 1e-11,

        dtmax = 0.1,

        saveat = 0.05,

        maxiters =
            Int(1e7),

        callback =
            cb_budget_ft
    )
nothing
```

```@example sir3
t_ft_ut =
    sol_ft_ut.t

S_ft_ut =
    sol_ft_ut[1, :]

I_ft_ut =
    sol_ft_ut[2, :]

C_ft_ut =
    sol_ft_ut[13, :]

χ1_ft_ut =
    sol_ft_ut[3, :]

χ2_ft_ut =
    sol_ft_ut[4, :]

θ1hat_ft_ut =
    sol_ft_ut[11, :]

θ2hat_ft_ut =
    sol_ft_ut[12, :]

egt_v_ft_ut =
    exp.(
        -γ .* t_ft_ut
    )

Shat_ft_ut =
    θ1hat_ft_ut .+
    χ1_ft_ut

Ihat_ft_ut =
    χ2_ft_ut .+
    egt_v_ft_ut .* θ2hat_ft_ut

Idag_ft_ut =
    IDAG.(
        Shat_ft_ut,
        Ihat_ft_ut,
        C_ft_ut
    )

u_ft_ut =
    ifelse.(
        C_ft_ut .> 0.0,
        ft_control.(
            Shat_ft_ut,
            Ihat_ft_ut,
            Idag_ft_ut
        ),
        0.0
    )

e_S_ft_ut =
    abs.(
        S_ft_ut .-
        Shat_ft_ut
    )

e_I_ft_ut =
    abs.(
        I_ft_ut .-
        Ihat_ft_ut
    )
nothing
```

Performance indicators for the modified law:

```@example sir3
println(
    "EKF[u3_tilde] max |S-Ŝ| = ",
    maximum(e_S_ekf_ut)
)

println(
    "EKF[u3_tilde] max |I-Î| = ",
    maximum(e_I_ekf_ut)
)

max_I_ft_ut, idx_ft_ut =
    findmax(I_ft_ut)

println(
    "FT[u3_tilde] max I(t) = ",
    max_I_ft_ut,
    " at t = ",
    t_ft_ut[idx_ft_ut]
)

println(
    "FT[u3_tilde] max |S-Ŝ| = ",
    maximum(e_S_ft_ut)
)

println(
    "FT[u3_tilde] max |I-Î| = ",
    maximum(e_I_ft_ut)
)

_, _, dt_ekf_ut =
    crossing_stats(
        t_ekf_ut,
        I_ekf_ut,
        Ibar
    )

_, _, dt_ft_ut =
    crossing_stats(
        t_ft_ut,
        I_ft_ut,
        Ibar
    )

println(
    "EKF[u3_tilde]: I(t) above Ibar for Δt ≈ ",
    dt_ekf_ut
)

println(
    "FT[u3_tilde]: I(t) above Ibar for Δt ≈ ",
    dt_ft_ut
)
```

# Numerical comparison

The following plots compare

1. the susceptible population,
2. the infected population,
3. the feedback control,
4. the remaining intervention budget.

The reference trajectory is shown together with the EKF-based and
finite-time-observer-based closed-loop trajectories.

```@example sir3
COLOR_REF =
    :black

COLOR_EKF =
    RGB(
        0.14,
        0.44,
        0.64
    )

COLOR_FT =
    RGB(
        0.75,
        0.22,
        0.17
    )

COLOR_SH =
    RGB(
        0.15,
        0.68,
        0.38
    )
nothing
```

Common plotting options are

```@example sir3
common = (
    guidefontsize = 13,
    tickfontsize = 10,
    legendfontsize = 12,
    framestyle = :box,
    grid = false,
    xticks = 0:10:tf,
    xlims = (0, tf),
    left_margin = 0Plots.mm,
    right_margin = 1Plots.mm,
    dpi = 500
)
nothing
```

### Susceptible population

```@example sir3
pS =
    plot(
        t_ref,
        S_ref;

        lw = 3.0,
        color = COLOR_REF,

        label =
            L"S_{\mathrm{ref}}",

        legend = :topright,
        legend_column = 3,

        top_margin =
            3Plots.mm,

        bottom_margin =
            1Plots.mm,

        common...
    )

plot!(
    pS,
    t_ekf,
    S_ekf;

    lw = 2.0,
    color = COLOR_EKF,
    ls = :solid,

    label =
        L"S_{\mathrm{EKF,true}}"
)

plot!(
    pS,
    t_ekf,
    Shat_ekf;

    lw = 2.0,
    color = COLOR_EKF,
    ls = :dot,

    label =
        L"\hat S_{\mathrm{EKF}}"
)

plot!(
    pS,
    t_ft,
    S_ft;

    lw = 2.0,
    color = COLOR_FT,
    ls = :solid,

    label =
        L"S_{\mathrm{FT,true}}"
)

plot!(
    pS,
    t_ft,
    Shat_ft;

    lw = 2.0,
    color = COLOR_FT,
    ls = :dot,

    label =
        L"\hat S_{\mathrm{FT}}"
)

hline!(
    pS,
    [Sh];

    lw = 1.2,
    color = COLOR_SH,
    ls = :dashdot,

    label =
        L"S_h"
)
```

### Infected population

```@example sir3
pI =
    plot(
        t_ref,
        I_ref;

        lw = 3.0,
        color = COLOR_REF,

        label =
            L"I_{\mathrm{ref}}",

        legend = :topright,
        legend_column = 3,

        top_margin =
            1Plots.mm,

        bottom_margin =
            1Plots.mm,

        common...
    )

plot!(
    pI,
    t_ekf,
    I_ekf;

    lw = 2.0,
    color = COLOR_EKF,
    ls = :solid,

    label =
        L"I_{\mathrm{EKF,true}}"
)

plot!(
    pI,
    t_ekf,
    Ihat_ekf;

    lw = 2.0,
    color = COLOR_EKF,
    ls = :dot,

    label =
        L"\hat I_{\mathrm{EKF}}"
)

plot!(
    pI,
    t_ft,
    I_ft;

    lw = 2.0,
    color = COLOR_FT,
    ls = :solid,

    label =
        L"I_{\mathrm{FT,true}}"
)

plot!(
    pI,
    t_ft,
    Ihat_ft;

    lw = 2.0,
    color = COLOR_FT,
    ls = :dot,

    label =
        L"\hat I_{\mathrm{FT}}"
)

plot!(
    pI,
    t_ft,
    [
        C_ft[k] > 0.0 ?
        Idag_ft[k] :
        NaN
        for k in eachindex(t_ft)
    ];

    lw = 1.6,
    color = COLOR_SH,
    ls = :dash,

    label =
        L"I^{\dagger}(\hat I,\hat S,C)"
)
```

### Control

```@example sir3
pU =
    plot(
        t_ref,
        u_ref;

        lw = 3.0,
        color = COLOR_REF,

        label =
            L"u_{\mathrm{ref}}",

        ylims =
            (-0.02, 1.02),

        legend =
            :topright,

        top_margin =
            1Plots.mm,

        bottom_margin =
            1Plots.mm,

        common...
    )

plot!(
    pU,
    t_ekf,
    u_ekf;

    lw = 2.0,
    color = COLOR_EKF,
    ls = :solid,

    label =
        L"u_{\mathrm{EKF}}"
)

plot!(
    pU,
    t_ft,
    u_ft;

    lw = 2.0,
    color = COLOR_FT,
    ls = :solid,

    label =
        L"u_{\mathrm{FT}}"
)
```

### Remaining budget

```@example sir3
pC =
    plot(
        t_ref,
        C_ref;

        lw = 3.0,
        color = COLOR_REF,

        label =
            L"C_{\mathrm{ref}}",

        xlabel =
            L"t",

        ylims =
            (-0.5, Qbudget + 1),

        legend =
            :topright,

        top_margin =
            1Plots.mm,

        bottom_margin =
            4Plots.mm,

        common...
    )

plot!(
    pC,
    t_ekf,
    C_ekf;

    lw = 2.0,
    color = COLOR_EKF,
    ls = :solid,

    label =
        L"C_{\mathrm{EKF}}"
)

plot!(
    pC,
    t_ft,
    C_ft;

    lw = 2.0,
    color = COLOR_FT,
    ls = :solid,

    label =
        L"C_{\mathrm{FT}}"
)

hline!(
    pC,
    [0.0];

    lw = 1.0,
    color = :gray,
    ls = :dot,

    label =
        false
)
```

### Modified feedback trajectories

The closed-loop trajectories obtained with $\widetilde u_3$ are added
to the same panels using dashed line styles.

```@example sir3
plot!(
    pS,
    t_ekf_ut,
    S_ekf_ut;

    lw = 2.0,
    color = COLOR_EKF,
    ls = :dash,

    label =
        L"S_{\mathrm{EKF,true},\widetilde u_3}"
)

plot!(
    pS,
    t_ekf_ut,
    Shat_ekf_ut;

    lw = 2.0,
    color = COLOR_EKF,
    ls = :dashdot,

    label =
        L"\hat S_{\mathrm{EKF},\widetilde u_3}"
)

plot!(
    pS,
    t_ft_ut,
    S_ft_ut;

    lw = 2.0,
    color = COLOR_FT,
    ls = :dash,

    label =
        L"S_{\mathrm{FT,true},\widetilde u_3}"
)

plot!(
    pS,
    t_ft_ut,
    Shat_ft_ut;

    lw = 2.0,
    color = COLOR_FT,
    ls = :dashdot,

    label =
        L"\hat S_{\mathrm{FT},\widetilde u_3}"
)
```

```@example sir3
plot!(
    pI,
    t_ekf_ut,
    I_ekf_ut;

    lw = 2.0,
    color = COLOR_EKF,
    ls = :dash,

    label =
        L"I_{\mathrm{EKF,true},\widetilde u_3}"
)

plot!(
    pI,
    t_ekf_ut,
    Ihat_ekf_ut;

    lw = 2.0,
    color = COLOR_EKF,
    ls = :dashdot,

    label =
        L"\hat I_{\mathrm{EKF},\widetilde u_3}"
)

plot!(
    pI,
    t_ft_ut,
    I_ft_ut;

    lw = 2.0,
    color = COLOR_FT,
    ls = :dash,

    label =
        L"I_{\mathrm{FT,true},\widetilde u_3}"
)

plot!(
    pI,
    t_ft_ut,
    Ihat_ft_ut;

    lw = 2.0,
    color = COLOR_FT,
    ls = :dashdot,

    label =
        L"\hat I_{\mathrm{FT},\widetilde u_3}"
)

plot!(
    pI,
    t_ft_ut,
    [
        C_ft_ut[k] > 0.0 ?
        Idag_ft_ut[k] :
        NaN
        for k in eachindex(t_ft_ut)
    ];

    lw = 1.4,
    color = COLOR_SH,
    ls = :dot,

    label =
        L"I^{\dagger}_{\widetilde u_3}"
)
```

```@example sir3
plot!(
    pU,
    t_ekf_ut,
    u_ekf_ut;

    lw = 2.0,
    color = COLOR_EKF,
    ls = :dash,

    label =
        L"u_{\mathrm{EKF},\widetilde u_3}"
)

plot!(
    pU,
    t_ft_ut,
    u_ft_ut;

    lw = 2.0,
    color = COLOR_FT,
    ls = :dash,

    label =
        L"u_{\mathrm{FT},\widetilde u_3}"
)
```

```@example sir3
plot!(
    pC,
    t_ekf_ut,
    C_ekf_ut;

    lw = 2.0,
    color = COLOR_EKF,
    ls = :dash,

    label =
        L"C_{\mathrm{EKF},\widetilde u_3}"
)

plot!(
    pC,
    t_ft_ut,
    C_ft_ut;

    lw = 2.0,
    color = COLOR_FT,
    ls = :dash,

    label =
        L"C_{\mathrm{FT},\widetilde u_3}"
)
```
