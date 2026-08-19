# FT_Observer.jl

Documentation for `FT_Observer.jl`.

## Problem setting

`FT_Observer.jl` contains the numerical implementations accompanying the paper

*Observer synthesis and peak reduction for the SIR model with output feedback under budget-constrained interventions*.

The objective is to reduce the epidemic peak in a controlled SIR model when the susceptible and infected populations are not directly available for feedback. Instead, only the incidence rate

```math
y(t)=\beta S(t)I(t)(1-u(t))
```

is assumed to be measured, while the intervention satisfies a prescribed budget constraint.

We consider the controlled SIR model

```math
\begin{aligned}
\dot S(t)
    &= -\beta S(t)I(t)\bigl(1-u(t)\bigr),\\
\dot I(t)
    &= \beta S(t)I(t)\bigl(1-u(t)\bigr)-\gamma I(t),
\end{aligned}
```

where

- $S(t)$ is the susceptible population,
- $I(t)$ is the infected population,
- $\beta>0$ is the transmission rate,
- $\gamma>0$ is the recovery rate,
- $u(t)\in[0,1]$ denotes the intervention intensity.

The intervention is subject to the budget constraint

```math
\int_0^{t_f}u(t)\,dt\leq Q.
```

## Finite-time state estimation

The package implements a finite-time observer that reconstructs the susceptible and infected populations from incidence measurements. The observer reformulates state reconstruction as an estimation problem for the unknown initial conditions and uses Dynamic Regressor Extension and Mixing (DREM) to obtain finite-time, monotone parameter convergence.

The reconstructed states are then used in closed loop with the epidemic-control strategies studied in the paper.

## Common quantities

The herd-immunity threshold is

```math
S_h=\frac{\gamma}{\beta}.
```

The nominal optimal infection threshold $I^\star$ is computed from the estimated initial conditions $(\hat\theta_1,\hat\theta_2)$ and the available budget $Q$ as described in the paper.

The reference trajectory under perfect state information serves as a benchmark for the output-feedback strategies.

## Numerical comparison framework

Each feedback section contains:

- **Reference trajectory**: the optimal trajectory under perfect state information, computed in three phases (uncontrolled epidemic, intervention arc, final uncontrolled trajectory).
- **Finite-time observer**: closed-loop simulation using the DREM-based state estimator.
- **Extended Kalman Filter**: a continuous-time EKF under noisy incidence measurements for comparison.
- **Performance indicators**: estimation errors, epidemic peak, time spent above the infection threshold, control signal, and remaining budget.

The simulations illustrate the evolution of the true and estimated SIR states, the control signal, the remaining budget, estimation errors, and the effect of the different feedback constructions on epidemic-peak reduction and budget consumption.

## Feedback strategies

The numerical documentation is organized around three feedback laws:

- **[Feedback 1](feedback1.md)** — uses the threshold
  $I^\star(\hat\theta_1,\hat\theta_2,Q)$, obtained from the estimated initial conditions.

- **[Feedback 2](feedback2.md)** — uses the constant threshold
  $I_l$, which is independent of the estimated initial conditions and provides a simple suboptimal feedback strategy.

- **[Feedback 3](feedback3.md)** — uses the adaptive threshold
  $I^\dagger(\hat I,\hat S,C)$, which depends on the estimated state and the remaining intervention budget.

For each strategy, both the original feedback law $u_i$ and its modified version $\widetilde u_i$ are considered. The modified laws include an additional multiplicative correction designed to keep the infected trajectory closer to the corresponding target level and to improve the use of the available intervention budget.