# General Log-Linear Model (GLLM)

Let us now assume a log-linear variant of the GLM, where

$$
\mathbf{Y} = \exp.\left(\mathbf{X}\mathbf{B}\right) + \mathbf{U} ~,
$$

where we use $\exp .$ (Julia's convention) to denote that the exponential
is applied independently to each element of the matrix (_i.e._, it is
_not_ a matrix exponential).

In the rest of this section, we will write the _model fit_
$\mathbf{Z} = \exp.\left(\mathbf{X}\mathbf{B}\right)$ for conciseness.

The problem with this model is that no closed-form solution for the
maximum-likelihood $\hat{\mathbf{B}}$ (or posterior distribution
$p(\mathbf{B}\mid\mathbf{Y})$) exists.

We propose two algorithms to estimate the posterior distribution
over $\hat{\mathbf{B}}$, assuming a flat prior.

## Laplace approximation

The first method uses Laplace's approximation.

:::{admonition} Laplace approximation
In the general case of a density function $p\left(\mathbf{x}\right)$, we would
compute the mode of the distribution and the curvature of its negative
log-density evaluated at the mode:

$$
\begin{align*}
\hat{\boldsymbol{\mu}} & = \arg\min_{\mathbf{x}} \left\{ -\ln p\left(\mathbf{x}\right) \right\} ~,\\
\hat{\boldsymbol{\Lambda}} & = \left.\frac{\partial}{\partial\mathbf{x}}\right|_{\hat{\boldsymbol{\mu}}} \left\{-\ln p\left(\mathbf{x}\right)\right\}
~ .
\end{align*}
$$

The posterior distribution is then approximated by the Gaussian distribution

$$
p\left(\mathbf{x}\right) \approx p\left(\mathbf{x} ~\middle|~ \hat{\boldsymbol{\mu}}, \hat{\boldsymbol{\Lambda}}^{-1}\right)
$$
:::

In the case of our log-linear model, we need an iterative prcedure to find
the mode. We will use a second order optimization scheme which will also provides
us with the curvature of the log-likelihood, which is used in our approximation.

The negative log-likelihood of our model is

$$
\begin{align*}
\mathrm{NLL}
& = -\ln p\left(\mathbf{Y} ~\middle|~ \mathbf{Z}\right) \\
& = \frac{1}{2}\left[\operatorname{Tr}\left(
    \left(\mathbf{Z} - \mathbf{Y} \right)
    \mathbf{S}^{-1}
    \left(\mathbf{Z}- \mathbf{Y} \right)^\mathrm{T}
    \mathbf{C}^{-1}
\right)
+ N\ln\left|\mathbf{C}\right|
+ M\ln\left|\mathbf{S}\right|
+ NM\ln 2\pi
\right]
~.
\end{align*}
$$

Differentiating with respect to $\mathbf{B}$, we find

$$
\mathbf{G} = \frac{\partial\mathrm{NLL}}{\partial\mathbf{B}}
=
\mathbf{X}^\mathrm{T}\left[
    \mathbf{Z}
    \odot
    \left(
    \mathbf{C}^{-1}
    \left(\mathbf{Z} - \mathbf{Y} \right)
    \mathbf{S}^{-1}
    \right)
\right]
~,
$$

and (assuming $\mathbf{S}=\mathbf{I}$ for simplicity)

$$
\mathbf{H}_n =
\frac{\partial\mathrm{NLL}}{\partial\mathbf{b}_n\partial\mathbf{b}_n^\mathrm{T}}
 =
\mathbf{X}^\mathrm{T}\left[
    \left(
    \mathbf{z}_n
    \mathbf{z}_n^\mathrm{T}
    \right)
    \odot
    \mathbf{C}^{-1}
    +
    \operatorname{diag}\left(
        \mathbf{z}_n
        \odot
        \left[
        \mathbf{C}^{-1}
        \left(\mathbf{z}_n - \mathbf{y}_n\right)
        \right]
    \right)
\right]
\mathbf{X}
~.
$$

:::{admonition} Symbolic verification
:class: dropdown
```matlab
% ----------------------------------------------------------------------
% Domain size
% -----------
N   = 1;                        % Number of voxels
M   = 3;                        % Number of model variables
K   = 2;                        % Number of model parameters
% ----------------------------------------------------------------------
% Input matrices
% --------------
X   = sym('X', [M K], 'real');  % Design matrix
B   = sym('B', [K N], 'real');  % Model parameters
Y   = sym('Y', [M N], 'real');  % Observations
A   = sym('C', [M M], 'real');  % Inverse covariance matrix
A   = A + A';                   % (symmetric)
% ----------------------------------------------------------------------
% Forward model
% -------------
Z   = exp(X * B);               % Model fit
R   = Z - Y;                    % Residuals
NLL = 0.5 * trace(R*R'*A);      % Negative log-likelihood (terms in B)

% ----------------------------------------------------------------------
% Compute gradient
% ----------------
G = sym([]);
for k=1:K
    for n=1:N
        G(k,n) = diff(NLL, B(k,n));
    end
end

% Analytical gradient
% -------------------
GG = X' * (Z .* (A * R));

assert(~any(any(simplify(G-GG, 100))));

% ----------------------------------------------------------------------
% Compute Hessian
% ---------------
H = sym([]);
for k=1:K
    for l=1:K
        H(k,l) = diff(G(k,1), B(l,1));
    end
end

% Analytical Hessian
% ------------------
z  = Z(:,1);
r  = R(:,1);
HH = X' * ((z*z') .* A  + diag(z .* (A * r))) * X;

assert(~any(any(simplify(H-HH, 100))));
```
:::

Still assuming $\mathbf{S}=\mathbf{I}$, Newton's algorithm consists of
taking the step

$$
\mathbf{b}_n = \mathbf{b}_n - \mathbf{H}_n^{-1}\mathbf{g}_n ~.
$$

However, our Hessian $\mathbf{H}_n$ is not always positive definite, which can
lead to moving in the wrong direction (ascent instead of descent). To solve
this problem, the Gauss-Newton algorithm substitutes $\mathbf{g}_n = \mathbf{0}$
in the Hessian, resulting in the positive definite approximation

$$
\mathbf{H}_n \approx
\mathbf{X}^\mathrm{T}\Big[
    \left(
    \mathbf{z}_n
    \mathbf{z}_n^\mathrm{T}
    \right)
    \odot
    \mathbf{C}^{-1}
\Big]
\mathbf{X}
~.
$$

We found in {cite:p}`balbastre2021model` that the Gauss-Newton step
often overshoots in log-linear model, but that the following positive
definite approximation does not

$$
\mathbf{H}_n \approx
\mathbf{X}^\mathrm{T}\left[
    \left(
    \mathbf{z}_n
    \mathbf{z}_n^\mathrm{T}
    \right)
    \odot
    \mathbf{C}^{-1}
    +
    \operatorname{diag}\left(
        \mathbf{z}_n
        \odot
        \left[
        \mathbf{C}^{-1}
        \left|\mathbf{z}_n - \mathbf{y}_n\right|
        \right]
    \right)
\right]
\mathbf{X}
~,
$$

where $\left|\cdot\right|$ denotes the element-wise absolute value.
The second term may be weighted by a constant $\alpha \in [0, 1]$,
thereby interpolating between the (faster) Gauss-Newton and (more robust)
loaded Hessian.


After convergence, the optimum is reached and the Laplace approximation of
the posterior is

$$
\begin{align*}
p\left(\mathbf{b}_n\mid\mathbf{y}_n\right)
& \approx \mathcal{N}\left(
    \mathbf{b}_n
    ~\middle|~
    \hat{\mathbf{b}}_n,
    \hat{\boldsymbol{\Sigma}}_n
\right)
\\
\hat{\boldsymbol{\Sigma}}_n^{-1}
& =
\mathbf{X}^\mathrm{T}\Big[
    \left(
    \mathbf{z}_n
    \mathbf{z}_n^\mathrm{T}
    \right)
    \odot
    \mathbf{C}^{-1}
\Big]\mathbf{X}
\end{align*}
$$

### Expectation of a LogNormal random variable

:::{admonition} Moments of the LogNormal distribution
Let $\mathbf{x} \sim \mathcal{N}\left(\boldsymbol\mu, \boldsymbol\Sigma\right)$
be a Gaussian variable. The first and second moments of the variable
$\mathbf{y} = \exp.\left(\mathbf{x}\right)$ are (see {cite:p}`halliwell2015lognormal`)

$$
\begin{align*}
\mathbb{E}\left[\mathbf{y}\right] & = \exp.\left(\boldsymbol\mu + \frac{1}{2}\operatorname{diag}\left(\boldsymbol\Sigma\right)\right) ~, \\
\mathbb{E}\left[\mathbf{y}\mathbf{y}^\mathrm{T}\right] & = \left(\exp.\left(\boldsymbol\mu\right)\exp.\left(\boldsymbol\mu\right)^\mathrm{T}\right) \odot \exp.\left(\boldsymbol\Sigma\right) ~.
\end{align*}
$$
:::

We can use this identity to compute the _expected_ model fit
$\mathbb{E}\left[\exp.\left(\mathbf{X}\mathbf{B}\right)\right]$
under the approximate posterior.

## Variational Bayesian optimization

An alternative approach consists of
1. Encode the approximate posterior as a Gaussian

$$
q\left(\mathbf{b}_n\right)
= \mathcal{N}\left(
    \mathbf{b}_n
    ~\middle|~
    \hat{\mathbf{b}}_n,
    \hat{\boldsymbol{\Sigma}}_n
\right)
$$

2. Find the parameters $\hat{\mathbf{b}}_n$ and $\hat{\boldsymbol{\Sigma}}_n$
that makes $q$ closest to the true posterior in terms of the Kullback-Leibler
divergence $KL\left(q ~\middle\|~ p\right)$.

This approach, which replaces _computing_ the posterior with _optimizing_
a _family_ of posteriors, is known as Variational Bayes.

## References
```{bibliography}
:filter: docname in docnames
:style: plain
```
