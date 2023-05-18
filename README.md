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

Numerical optimization procedure is performed by gradient descent or its variants, by using direct o adjoint gradient computation.

## Scripts and their purpose

| Script | Purpose |
-----

