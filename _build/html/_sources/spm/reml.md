(spm:reml)=
# Heteroscedasticity and ReML

(spm:reml:heteroscedasticity)=
## Heteroscedasticity

In the {ref}`spm:glm` defined in the previous section, the covariance
matrix across rows of $\mathbf{Y}$ must be known. In practice, however,
it is not and is often simply assumed to be a (scaled) identity
$\sigma^2 \mathbf{I}$. In other words, all measurements are supposed to
be independant and with the same variance $\sigma^2$ ({term}`homoscedasticity`).
Note that in this case, the exact value of $\sigma^2$ is not required to compute the
{term}`maximum likelihood` estimate of the model parameters $\mathbf{B}$.

In cases where the homoscedasticity assumption does not hold
({term}`heteroscedasticity`), the covariance matrix must be _estimated_.

Naively, one may first estimate the maximum-likelihood parameter
$\hat{\mathbf{B}}$, assuming homoscedasticity, before estimating the
maximum-likelihood covariance $\hat{\mathbf{C}}$ conditioned on $\hat{\mathbf{B}}$,
reestimate $\hat{\mathbf{B}}$ assuming $\hat{\mathbf{C}}$, and iterate.
This is a _joint_ maximum-likelihood estimation

$$
\hat{\mathbf{B}}, \hat{\mathbf{C}} = \operatorname{arg}\max_{\mathbf{B},\mathbf{C}}
p\left(\mathbf{Y}~\middle|~\mathbf{B},\mathbf{C}\right) ~.
$$

However, the optimization problem is not jointly convex in general and
may therefore not converge.

(spm:reml:reml)=
## Restricted Maximum-Likelihood (ReML)

A more elegant solution consists of only estimating the maximum-likelihood
covariance, and marginalize the model parameters

$$
\hat{\mathbf{C}}
= \operatorname{arg}\max_{\mathbf{C}} p\left(\mathbf{Y}~\middle|~\mathbf{C}\right)
= \operatorname{arg}\max_{\mathbf{C}} \int p\left(\mathbf{Y}~\middle|~\mathbf{B},\mathbf{C}\right)~p\left(\mathbf{B}\right)~\mathrm{d}\mathbf{B} ~.
$$

This requires defining a {term}`prior` over $\mathbf{B}$, which we can
again assume noninformative ($p\left(\mathbf{B}\right)\propto 1$).

This optimization is carried out using an {term}`Expectation-Maximization` (EM)
procedure. Starting with an initial estimate $\mathbf{C}^\ast$, we first compute the
posterior distribution over the latent quantity $\mathbf{B}$

$$
p\left(\mathbf{B} ~\middle|~ \mathbf{Y}, \mathbf{C}^\ast\right)
\propto
\mathcal{N}\left(\mathbf{Y}  ~\middle|~ \mathbf{X}\mathbf{B}, \mathbf{C}^\ast, \mathbf{S}\right) ~.
$$

Following the argument made in section {ref}`spm:glm:inference:bayes`,
we find the this distribution is again Matrix-Gaussian with

$$
q^\ast\left(\mathbf{B}\right)
= p\left(\mathbf{B}\mid\mathbf{Y}, \mathbf{C}^\ast\right)
= \mathcal{N}\left(
    \mathbf{B}
    ~\middle|~
    \underbrace{\left(\mathbf{X}^\mathrm{T}\left.\mathbf{C}^\ast\right.^{-1}\mathbf{X}\right)^{-1}
  \mathbf{X}^\mathrm{T}\left.\mathbf{C}^\ast\right.^{-1}\mathbf{Y}}_{\mathbf{B}^\ast},
    \underbrace{\left(\mathbf{X}^\mathrm{T}\left.\mathbf{C}^\ast\right.^{-1}\mathbf{X}\right)^{-1}}_{\boldsymbol{\Sigma}^\ast},
    \mathbf{S}
\right)
$$

(**E-step**) We then compute a functional in $\mathbf{C}$ that bounds
the true log-likelihood $\ln p\left(\mathbf{Y}~\middle|~\mathbf{C}\right)$ from below:

$$
\begin{align*}
Q\left(\mathbf{C} ~\middle|~ \mathbf{C}^\ast\right)
& = \mathbb{E}_{q^\ast}\left[\ln p\left(\mathbf{Y},\mathbf{B} ~\middle|~ \mathbf{C}\right)\right] \\
& = -\frac{1}{2}\mathbb{E}_{q^\ast}\left[\operatorname{Tr}\left(\mathbf{C}^{-1}\left(\mathbf{X}\mathbf{B}-\mathbf{Y}\right)\mathbf{S}^{-1}\left(\mathbf{X}\mathbf{B}-\mathbf{Y}\right)^\mathrm{T}\right)
+ N\ln\left|\mathbf{C}\right|
+ M\ln\left|\mathbf{S}\right|
+ NM\ln 2\pi\right]
~.
\end{align*}
$$

Keeping only terms that depend on $\mathbf{C}$ and computing the expectation
accordingly, we find

$$
\begin{align*}
-2Q\left(\mathbf{C} ~\middle|~ \mathbf{C}^\ast\right)
=
\operatorname{Tr}\left(
    \mathbf{C}^{-1}\left(\mathbf{X}\mathbf{B}^\ast-\mathbf{Y}\right)
    \mathbf{S}^{-1}\left(\mathbf{X}\mathbf{B}^\ast-\mathbf{Y}\right)^\mathrm{T}
\right)
+
N \operatorname{Tr}\left(
    \boldsymbol{\Sigma}^\ast
    \left(\mathbf{X}^\mathrm{T}\mathbf{C}^{-1}\mathbf{X}\right)
\right)
+
N\ln\left|\mathbf{C}\right| ~.
\end{align*}
$$

(**M-step**) We finally differentiate with respect to the
_{term}`precision matrix`_ $\mathbf{C}^{-1}$ and solve for zero, finding

$$
\mathbf{C} =
\frac{1}{N}
\left(\mathbf{X}\mathbf{B}^\ast-\mathbf{Y}\right)
\mathbf{S}^{-1}
\left(\mathbf{X}\mathbf{B}^\ast-\mathbf{Y}\right)^\mathrm{T}
+ \mathbf{X}\boldsymbol{\Sigma}^\ast\mathbf{X}^\mathrm{T} ~.
$$

This new estimate can now replace the previous $\mathbf{C}^\ast$, and
the EM procedure is repeated until convergence.

## Low-rank parameterization

In practice, the full covariance matrix possesses $\frac{N(N+1)}{2}$ free parameters,
which are often too many to estimate, leading to overfitting and
under-estimated true variance. The alternative, implemented in SPM, consists
in encoding the covariance matrix using a linear combination of "basis functions"

$$
\mathbf{C} = \sum_{j=1}^{J} h_j \mathbf{Q}_j ~,
$$

where the $\left\{\mathbf{Q}_j\right\}_{j=1}^J$ are the user-defined basis
functions and $\mathbf{h} \in \mathbb{R}^J ~\left(\text{with}~J \ll \frac{N(N+1)}{2}\right)$
contains the free-parameters to optimize.

Contrary to the full covariance case, there is no simple closed-form solution
for $\mathbf{h}$. Instead, we replace the closed-form M-step with an
iterative M-step using Gauss-Newton optimization.

Let us write $\mathcal{L} = -2NQ$. We differentiate once with respect to
$\mathbf{C}^{-1}$ and find

$$
\frac{\partial \mathcal{L}}{\partial\mathbf{C}^{-1}}
= \underbrace{\frac{1}{N}\left(\mathbf{X}\mathbf{B}^\ast-\mathbf{Y}\right)
\mathbf{S}^{-1}
\left(\mathbf{X}\mathbf{B}^\ast-\mathbf{Y}\right)^\mathrm{T}}_{\mathbf{R}}
+ \underbrace{\left(\mathbf{X}\boldsymbol{\Sigma}^\ast\mathbf{X}^\mathrm{T}\right)}_{\mathbf{U}}
- \mathbf{C}
$$

Let us note that the first two terms ($\mathbf{R}+\mathbf{U}$)
do not depend on $\mathbf{C}$. We may then apply the chain rule to
obtain the derivative of the functional with respect to one of the parameters

$$
\begin{align*}
\frac{\partial \mathcal{L}}{\partial h_j}
& {}= \operatorname{Tr}\left(
    \frac{\partial \mathbf{C}^{-1}}{\partial h_j}^\mathrm{T}
    \frac{\partial \mathcal{L}}{\partial\mathbf{C}^{-1}}
\right)
{}= \operatorname{Tr}\left(
    \mathbf{C}^{-1}\frac{\partial \mathbf{C}}{\partial h_j}\mathbf{C}^{-1}
    \frac{\partial \mathcal{L}}{\partial\mathbf{C}^{-1}}
\right)
{}= \operatorname{Tr}\left(
    \mathbf{C}^{-1}\mathbf{Q}_j\mathbf{C}^{-1}
    \frac{\partial \mathcal{L}}{\partial\mathbf{C}^{-1}}
\right)
\\
&{}= \operatorname{Tr}\left(
    \mathbf{Q}_j\left(
    \mathbf{C}^{-1}\left(\mathbf{R}+\mathbf{U}\right)\mathbf{C}^{-1}
    -
    \mathbf{C}^{-1}
    \right)
\right)
\end{align*}
$$

We can differentiate one more time (and take into accout the fact that
the covariance matrix, and thefore its gradient, is symmetric), to obtain

$$
\frac{\partial^2 \mathcal{L}}{\partial h_i\partial h_j}
=
\operatorname{Tr}\left(
    \left(
    \mathbf{C}^{-1}\mathbf{Q}_i\mathbf{C}^{-1}\mathbf{Q}_j
    \right)
    \left(
        2\mathbf{C}^{-1}\left(\mathbf{R}+\mathbf{U}\right) - \mathbf{I}
    \right)
\right)
$$

The resulting Hessian matrix is not always positive-definite. We further
simplify it by peforming Fisher's scoring, for which we assume that the
current estimate of $\mathbf{C}$ is optimal, and therefore substitute the
identity $\frac{\partial \mathcal{L}}{\partial\mathbf{C}^{-1}}=\mathbf{0}$
(which leads to $\mathbf{R}+\mathbf{U}=\mathbf{C}$) in the expression of the Hessian.
This results in

$$
\frac{\partial^2 \mathcal{L}}{\partial h_i\partial h_j}
\approx
\operatorname{Tr}\left(
    \mathbf{C}^{-1}\mathbf{Q}_i\mathbf{C}^{-1}\mathbf{Q}_j
\right)
$$

Let us write $\mathbf{g} = \frac{\partial \mathcal{L}}{\partial \mathbf{h}}$
and $\mathbf{H} = \frac{\partial \mathcal{L}}{\partial \mathbf{h}\partial \mathbf{h}^\mathrm{T}}$,
we can now take a _generalized_ M-step

$$
\mathbf{h}^\ast = \mathbf{h}^\ast - \mathbf{H}^{-1}\mathbf{g} ~,
$$

from which we build an updated covariance

$$
\mathbf{C}^\ast = \sum_{j=1}^{J} h^\ast_j \mathbf{Q}_j ~,
$$

before repeating the procedure until convergence.


### Computational efficiency

One may note that we are computing the gradient and Hessian at the previous
estimate $\mathbf{b}^\ast$. Therefore, $\mathbf{C}^\ast$ can be substituted
for $\mathbf{C}$ in the expressions of the gradient and Hessian, which leads
to some simplifications.

Let's first simplify the term
$\left.\mathbf{C}^\ast\right.^{-1}\mathbf{R}\left.\mathbf{C}^\ast\right.^{-1}$:

$$
\begin{align*}
N\left.\mathbf{C}^\ast\right.^{-1}\mathbf{R}\left.\mathbf{C}^\ast\right.^{-1}
& =
\left.\mathbf{C}^\ast\right.^{-1}
\left(\mathbf{X}\mathbf{B}^\ast-\mathbf{Y}\right)
\mathbf{S}^{-1}
\left(\mathbf{X}\mathbf{B}^\ast-\mathbf{Y}\right)^\mathrm{T}
\left.\mathbf{C}^\ast\right.^{-1}
\\
& =
\left.\mathbf{C}^\ast\right.^{-1}
\left(\mathbf{X}\left(\mathbf{X}^\mathrm{T}\left.\mathbf{C}^\ast\right.^{-1}\mathbf{X}\right)^{-1}
  \mathbf{X}^\mathrm{T}\left.\mathbf{C}^\ast\right.^{-1}\mathbf{Y}-\mathbf{Y}\right)
\mathbf{S}^{-1}
\left(\mathbf{X}\left(\mathbf{X}^\mathrm{T}\left.\mathbf{C}^\ast\right.^{-1}\mathbf{X}\right)^{-1}
  \mathbf{X}^\mathrm{T}\left.\mathbf{C}^\ast\right.^{-1}\mathbf{Y}-\mathbf{Y}\right)^\mathrm{T}
\left.\mathbf{C}^\ast\right.^{-1}
\\
& =
\left(\left.\mathbf{C}^\ast\right.^{-1}\mathbf{X}\left(\mathbf{X}^\mathrm{T}\left.\mathbf{C}^\ast\right.^{-1}\mathbf{X}\right)^{-1}
  \mathbf{X}^\mathrm{T}\left.\mathbf{C}^\ast\right.^{-1}-\left.\mathbf{C}^\ast\right.^{-1}\right)
\left(\mathbf{Y}\mathbf{S}^{-1}\mathbf{Y}^\mathrm{T}\right)
\left(\left.\mathbf{C}^\ast\right.^{-1}\mathbf{X}\left(\mathbf{X}^\mathrm{T}\left.\mathbf{C}^\ast\right.^{-1}\mathbf{X}\right)^{-1}
  \mathbf{X}^\mathrm{T}\left.\mathbf{C}^\ast\right.^{-1}-\left.\mathbf{C}^\ast\right.^{-1}\right)^\mathrm{T}
\end{align*}
$$

We note

$$
\mathbf{P} =
\left.\mathbf{C}^\ast\right.^{-1}\mathbf{X}\left(\mathbf{X}^\mathrm{T}\left.\mathbf{C}^\ast\right.^{-1}\mathbf{X}\right)^{-1}
  \mathbf{X}^\mathrm{T}\left.\mathbf{C}^\ast\right.^{-1}-\left.\mathbf{C}^\ast\right.^{-1}
$$

which gives us

$$
N\left.\mathbf{C}^\ast\right.^{-1}\mathbf{R}\left.\mathbf{C}^\ast\right.^{-1}
= \mathbf{P}\left(\mathbf{Y}\mathbf{S}^{-1}\mathbf{Y}^\mathrm{T}\right)\mathbf{P}^\mathrm{T}
~.
$$

Let us now simplify the term
$\left.\mathbf{C}^\ast\right.^{-1}\mathbf{U}\left.\mathbf{C}^\ast\right.^{-1}$:

$$
\begin{align*}
\left.\mathbf{C}^\ast\right.^{-1}\mathbf{U}\left.\mathbf{C}^\ast\right.^{-1}
& =
\left.\mathbf{C}^\ast\right.^{-1}
\mathbf{X}
\left(\mathbf{X}^\mathrm{T}\left.\mathbf{C}^\ast\right.^{-1}\mathbf{X}\right)^{-1}
\mathbf{X}^\mathrm{T}
\left.\mathbf{C}^\ast\right.^{-1}
\\
\Rightarrow~~
\left.\mathbf{C}^\ast\right.^{-1}\mathbf{U}\left.\mathbf{C}^\ast\right.^{-1} -
\left.\mathbf{C}^\ast\right.^{-1}
& =
\left.\mathbf{C}^\ast\right.^{-1}
\mathbf{X}
\left(\mathbf{X}^\mathrm{T}\left.\mathbf{C}^\ast\right.^{-1}\mathbf{X}\right)^{-1}
\mathbf{X}^\mathrm{T}
\left.\mathbf{C}^\ast\right.^{-1} -
\left.\mathbf{C}^\ast\right.^{-1}
= \mathbf{P}
\end{align*}
$$

Therefore

$$
\frac{\partial \mathcal{L}}{\partial h_j}
= \operatorname{Tr}\left(
    \mathbf{Q}_j\left(
    \mathbf{C}^{-1}\left(\mathbf{R}+\mathbf{U}\right)\mathbf{C}^{-1}
    -
    \mathbf{C}^{-1}
    \right)
\right)
= \operatorname{Tr}\left(
    \mathbf{Q}_j\left(
    \mathbf{P}\left(\mathbf{Y}\mathbf{S}^{-1}\mathbf{Y}^\mathrm{T}\right)\mathbf{P}^\mathrm{T}
    + \mathbf{P}
    \right)
\right)
$$

If we further write

$$
\mathbf{V} = \mathbf{P}\left(\mathbf{Y}\mathbf{S}^{-1}\mathbf{Y}^\mathrm{T}\right) + \mathbf{I}
$$

we have

$$
\frac{\partial \mathcal{L}}{\partial h_j}
=
\operatorname{Tr}\left(
    \mathbf{P}\mathbf{Q}_j\mathbf{V}
\right)
$$
