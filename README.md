# Observer Synthesis and Peak Reduction for the SIR Model

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://anasxbouali.github.io/FT_Observer/dev/)
[![HAL Preprint](https://img.shields.io/badge/HAL-05699171-blue)](https://hal.science/hal-05699171)

## Documentation

Comprehensive documentation, including mathematical details and usage examples, is available at:
- [Development Documentation](https://anasxbouali.github.io/FT_Observer/dev/)

The repository implements observers and controllers for the standard SIR model with a control input $u(t) \in [0,1]$ representing intervention intensity:

$$
\begin{aligned}
\dot{S} &= -\beta(1-u)SI \\
\dot{I} &= \beta(1-u)SI - \gamma I
\end{aligned}
$$

The observer estimates the state $(\hat{S}, \hat{I})$ from incidence rate measurements only, with guaranteed finite-time convergence.

### Control Strategies Implemented

The repository provides implementations of two main categories of Null-Singular-Null (NSN) feedback controls:

#### 1. Standard NSN Feedback Controls
- **Control $u_1$**: Threshold-based control with adaptive parameters $\hat{\theta}_1, \hat{\theta}_2$.
- **Control $u_2$**: Fixed threshold control with incidence rate monitoring.
- **Control $u_3$**: Budget-aware control with cumulative cost tracking.

#### 2. Modified NSN Feedback Controls
- **Modified Control $\tilde{u}_1$**: Improved version with dynamic scaling for enhanced sub-optimality.
- **Modified Control $\tilde{u}_2$**: Enhanced threshold control with incidence normalization.
- **Modified Control $\tilde{u}_3$**: Optimized budget consumption strategy.

*(Detailed mathematical definitions of these controls are available in the [documentation](https://anasxbouali.github.io/FT_Observer/dev/feedback1/) and the associated preprint).*

## Repository Structure

```text
FT_Observer/
├── docs/                    # Documenter.jl documentation source
├── examples/                # Jupyter notebooks and simulation scripts
├── src/                     # Core Julia implementations (observer, controls, simulation)
├── test/                    # Validation and test scripts
├── Project.toml             # Julia project dependencies
└── README.md                # This file
```

## Getting Started

### Prerequisites
- Julia 1.8 or higher
- Jupyter Notebook (optional, for interactive exploration of the `examples/` directory)

### Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/AnasXbouali/FT_Observer.git
   cd FT_Observer
   ```

2. Instantiate the Julia environment to install dependencies:
   ```julia
   using Pkg
   Pkg.instantiate()
   ```



## Citation

If you use this code or refer to the methods implemented in this repository, please cite the following preprint:

```bibtex
@inproceedings{bouali2026observer,
  title={Observer synthesis and peak reduction for the {SIR} model with output feedback under budget-constrained interventions},
  author={Bouali, Anas and Patelski, Rados{\l}aw and Rapaport, Alain and Efimov, Denis and Ushirobira, Rosane},
  booktitle={65th IEEE Conference on Decision and Control (CDC)},
  year={2026},
  hal_id={hal-05699171},
  url={https://hal.science/hal-05699171}
}
```
