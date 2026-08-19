module Common

using DifferentialEquations
using LinearAlgebra
using Random
using Plots
using LaTeXStrings
using Colors

export incidence, sigmoid_stable, ISTAR, crossing_stats
export R_ekf, Qproc_ekf, P0_ekf, Ŝ0_ekf, Î0_ekf
export λ_id, ρ_id, α_id, ε_id, Ŝ0_ft, Î0_ft, Y0
export COLOR_REF, COLOR_EKF, COLOR_FT, COLOR_SH, common

# Common numerical parameters
β = 0.8
γ = 0.2
S0 = 0.92
I0 = 0.04
Qbudget = 16.0
tf = 100.0

# Derived quantities
Sh = γ / β

# Safe threshold
safe_ = 1e-10

# Plotting defaults
default(fontfamily = "Computer Modern", framestyle = :box)

# Common colors
COLOR_REF = :black
COLOR_EKF = RGB(0.14, 0.44, 0.64)
COLOR_FT = RGB(0.75, 0.22, 0.17)
COLOR_SH = RGB(0.15, 0.68, 0.38)

# Common plotting options
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

# Incidence function
@inline incidence(S, I, u) = β * S * I * (1.0 - u)

# Measurement noise parameters
const NOISE_W = (
    0.17, 0.41, 0.93, 2.13, 3.71, 7.29, 11.93, 19.07
)

const NOISE_P = (
    0.00, 1.13, 2.47, 0.62, 4.01, 5.28, 1.97, 3.55
)

meas_noise_amp = 0.1

@inline function meas_noise(t)
    s = 0.0
    @inbounds for k in eachindex(NOISE_W)
        s += sin(NOISE_W[k] * t + NOISE_P[k])
    end
    return meas_noise_amp * sqrt(2.0 / length(NOISE_W)) * s
end

# Stable sigmoid function
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

# ISTAR function for computing optimal threshold from estimated initial conditions
function ISTAR(θ̂1::Real, θ̂2::Real)
    Ih = θ̂2 + θ̂1 - Sh - Sh * log(θ̂1 / Sh)
    Qcrit = (Ih - θ̂2) / (β * Sh * max(θ̂2, 1e-12))
    return Qbudget < Qcrit ? Ih / (1.0 + β * Sh * Qbudget) : θ̂2
end

# Crossing statistics function
function crossing_stats(tvec, Ivec, thresh)
    n = length(tvec)
    
    entry = findfirst(k -> Ivec[k] <= thresh && Ivec[k+1] > thresh, 1:n-1)
    entry === nothing && return (nothing, nothing, 0.0)
    
    interp(k) = tvec[k] + (thresh - Ivec[k]) * (tvec[k+1] - tvec[k]) / (Ivec[k+1] - Ivec[k])
    
    t_in = interp(entry)
    
    exitidx = findfirst(k -> k > entry && Ivec[k] >= thresh && Ivec[k+1] < thresh, 1:n-1)
    t_out = exitidx === nothing ? tvec[end] : interp(exitidx)
    
    return (t_in, t_out, t_out - t_in)
end

# EKF parameters
R_ekf = 1e-4
Qproc_ekf = Matrix(Diagonal([1e-2, 1e-2]))
P0_ekf = Matrix(Diagonal([1e-2, 1e-2]))
Ŝ0_ekf, Î0_ekf = 1.0, 0.0

# FT Observer parameters
λ_id = 1e-4
ρ_id = 1e6
α_id = 0.3
ε_id = 1e-6
Ŝ0_ft = 1.0
Î0_ft = 0.0
Y0 = β * S0 * I0

end # module
