# Glossary of Probability Theory

```{glossary}
Random Variable
    TODO

Probability Distribution
    TODO

PDF
Probability Density Function
    A probability density function over some space is a function that takes
    value in $\mathbb{R}^+$ and whose integral over that space is one:

    $$
        p: \Omega \rightarrow \mathbb{R}^+
        ~~\text{s.t.}~~ \int_\Omega p(x) \mathrm{d}x = 1
    $$

CDF
Cumulative Distribution Function
    Only defined for univariate distributions.
    It is the probability that a random variable takes any value below
    a given one:

    $$
        F: \left\{\begin{array}{rcl}
            \Omega & \rightarrow & [0, 1] \\
            x & \mapsto & \int_{-\infty}^x p(z) \mathrm{d}z
        \end{array}\right.
    $$

Kernel
    TODO

Partition
Partition Function
    TODO

Score
Score Function
    TODO

Entropy
    TODO

Divergence
    TODO

Kullback-Leibler Divergence
    TODO

Likelihood
    In Bayesian theory, a {term}`probability distribution`, $p(\mathbf{X}\mid\mathbf{Z})$,
    defined over some {term}`measurement` $\mathbf{X}$, conditioned on
    some {term}`unobserved` quantity $\mathbf{Z}$`.

    Can also more generally relate to any analytical {term}`probability density`.
    For example, the _joint likelihood_ of a model is the probability density
    defined over the union of all its variables.

Log-Likelihood
    Logarithm of the {term}`likelihood`.

NLL
Negative Log-Likelihood
    Negative {term}`log-likelihood`.

ML
Maximum Likelihood
Maximum-Likelihood
    TODO

Prior
    In Bayesian theory, a {term}`probability distribution`, $p(\mathbf{Z})$,
    defined over some {term}`unobserved` quantity $\mathbf{Z}$, that does not
    depend on any {term}`measurement`.


Posterior
    In Bayesian theory, a {term}`probability distribution`, $p(\mathbf{Z}\mid\mathbf{X})$,
    defined over some {term}`unobserved` quantity $\mathbf{Z}$, conditioned on a
    {term}`measurement` $\mathbf{X}$.

ReML
Restricted Maximum Likelihood
Restricted Maximum-Likelihood
    A {term}`maximum-likelihood` approach where some of the variables
    are {term}`marginalized`.

    See section _{ref}`spm:reml`_.

Bayes' Rule
    TODO

Bayesian Update
    TODO

GLM
General Linear Model
    See section _{ref}`spm:glm`_.

Homoscedastic
Homoscedasticity
    TODO

Heteroscedastic
Heteroscedasticity
    TODO

```
