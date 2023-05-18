# RPKOpt
## Reactor Point Kinetics parameters determination by ODE-constrained optimization

This collection of scripts solves the Reactor Point Kinetics equations


$$\frac{\text{d} N^{j}\left(t\right)}{\text{d} t} =\left(\varrho^{j}\left(t\right)-\beta_{\text{eff}}\right)\frac{N^{j}\left(t\right)}{\varLambda}+\sum_{i=1}^{m}\lambda_{i}c_{i}^{j}\left(t\right),$$

$$\frac{\text{d} c_{i}^{j}\left(t\right)}{\text{d} t} =\frac{\beta_{\text{eff},i}}{\varLambda}N^{j}\left(t\right)-\lambda_{i}c_{i}^{j}\left(t\right)\text{ for }i=1,\ldots,m,$$

$$N^{j}\left(0\right) =N_{\text{ini}}^{j},$$

$$c^{j}\left(0\right) =c_{i,\text{ini}}^{j}$$

and optimizes the reactor kinetic parameters $\beta_{\text{eff},i}$ by minimizing the $L_{2}$ error between the reactor power output $N^{j}$ predicted by numerical
solution of the model and the measurements $N^{j}_{\text{exp}}$

$$\underset{\vec{\beta}\in W_{\text{ad}}}{\min\,} J\left(N,c;\vec{\beta}\right)=\frac{1}{2}\sum_{j\in S}\int_{0}^{T_{j}}\left(N_{\text{exp}}^{j}\left(t\right)-N^{j}\left(t\right)\right)^{2}\text{d} t\$$

$S$ denotes a set of experiments with reactivity and power-output data loaded from the `MATLAB/datafiles` directory.

Numerical optimization procedure is performed by gradient descent (GD) or its variants, by using direct o adjoint gradient computation.

## Scripts and their purpose

The algorithm is implemented in MATLAB and is intended for interactive use. Hence, all scripts share a set of global variables also available from the
MATLAB workspace.

| Script | Purpose |
| ----- | ----- |
| `reactor_init` | initializes the computation. All parameters are set here. All variables except persistent variables are cleared. |
| `reactor_solve` | solves the primary problem for the current settings of $\beta_{\text{eff},i}$ (usually called from within other scripts) |
| `reactor_optimize` | peforms a given number of GD epochs, updates the plots of the results after each epoch |
| `reactor_randomize` | performs optimization by random shooting in the neighborhood of the current $\vec{\beta}$  |
| `reactor_seed` | performs multiple optimization runs (instances) from randomized starting points |
| `reactor_plot_loss` | plots multiple loss curves into one graph (results of `reactor_seed`) |
| `reactor_plot_paths` | plots the evolution of all $\beta_{\text{eff},i}$ in a composite figure |

`reactor_optimize` or `reactor_randomize` may be called multiple times after `reactor_init`

The script `reactor_study` uses persistent variables (named `PERSIST_....`) to perform a parametric study (the current version
iterates over different settings of solver accuracy, gradient computation method, and subset of experiments.
