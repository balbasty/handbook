# Multivariate L1 distribution

## Generalised multivariate Gaussian

In the univariate case, both the Gaussian and Laplace distributions can
be seen as special cases of the generalised Gaussian family defined by
the kernel

$$
f^\ast(\mathbf{x})
= \exp\left(-\left(\frac{\left|x - \mu\right|}{\alpha}\right)^\beta\right)
= \exp\left(-\left(\lambda\left(x - \mu\right)^2\right)^q\right)
~,
$$

where the Gaussian distribution uses $q=1, \lambda=\frac{1}{2\sigma^2}$
and the Laplace distribution uses $q=\frac{1}{2}, \lambda=\frac{1}{b^2}$.

Let's generalise the multivariate Gaussian in a similar way. The kernel
of the multivariate Gaussian is

$$
f^\ast(\mathbf{x}) = \exp\left(
    -\left(\mathbf{x}-\boldsymbol{\mu}\right)^T
    \boldsymbol{\Lambda}
    \left(\mathbf{x}-\boldsymbol{\mu}\right)
\right) ~.
$$

We define the generalised kernel

$$
f^\ast(\mathbf{x}) = \exp\left(
    -\left(\left(\mathbf{x}-\boldsymbol{\mu}\right)^T
    \boldsymbol{\Lambda}
    \left(\mathbf{x}-\boldsymbol{\mu}\right)\right)^q
\right) ~.
$$

## Multivariate Laplace

We find the multivariate Laplace[^mv1] kernel by setting $q=\frac{1}{2}$

$$
f^\ast(\mathbf{x}) = \exp\left(
    -\sqrt{\left(\mathbf{x}-\boldsymbol{\mu}\right)^T
    \boldsymbol{\Lambda}
    \left(\mathbf{x}-\boldsymbol{\mu}\right)}
\right) ~.
$$

The normalization constant $C$ of the probability density function
$p(\mathbf{x}) = Cf^\ast(\mathbf{x})$ is the
integral of the kernel over its support:

$$
C(\boldsymbol{\mu}, \boldsymbol{\Lambda})
= \int f^\ast(\mathbf{x}) \text{d}\mathbf{x}
= \int \exp\left(
    -\sqrt{\left(\mathbf{x}-\boldsymbol{\mu}\right)^T
    \boldsymbol{\Lambda}
    \left(\mathbf{x}-\boldsymbol{\mu}\right)}
\right) \text{d}\mathbf{x}
~.
$$

We perform the change of variable
$\mathbf{z} = \boldsymbol{\Lambda}^\frac{1}{2}\left(\mathbf{x}-\boldsymbol{\mu}\right)$,
whose Jacobian determinant is $\left|\boldsymbol{\Lambda}^\frac{1}{2}\right|$,
and integrate by substitution, leading to

$$
C(\boldsymbol{\mu}, \boldsymbol{\Lambda})
= \left|\boldsymbol{\Lambda}\right|^{-1} \int \exp\left(-\lVert\mathbf{z}\rVert_2\right) \text{d}\mathbf{z}
= 2 \left|\boldsymbol{\Lambda}\right|^{-1} \pi^\frac{n}{2} \frac{\Gamma\left(n\right)}{\Gamma\left(\frac{n}{2}\right)}
~
$$

where $n$ is the number of dimensions. The derivation of the final integral
is available [here](https://math.stackexchange.com/questions/2377072).
When $n=1$, we find the normalization constant of the univariate Laplace
distribution, as expected.

## Bivariate case

Let us focus on the bivariate case, and assume that $\boldsymbol{\mu} = \mathbf{0}$
and $\boldsymbol{\Lambda} = \mathbf{I}$, which we can always obtain by
whitening $\mathbf{x}$, _i.e._,
$\mathbf{z} = \boldsymbol{\Lambda}^\frac{1}{2}\left(\mathbf{x}-\boldsymbol{\mu}\right)$.
Let's further write $\mathbf{z} = [x, y]$.  We obtain the much simpler PDF

$$
p(x, y) = \frac{1}{\pi}
\exp\left(-\sqrt{x^2 + y^2}\right) ~.
$$

We would like to compute the marginal distribution $p(x) = \int p(x, y) \text{d}y$
and the conditional distribution $p(x | y) = \frac{p(x, y)}{\int p(x, y) \text{d}x}$.
Both are obtained by integrating the PDF over one of the two variables. We will
work with the kernel for clarity:

$$
\int f^\ast(x | y) \text{d}x
= \int  \exp\left(-\sqrt{x^2 + y^2}\right)\text{d}x
= 2y K_1\left(y\right)
~,
$$

where $K_1$ is the
[modified Bessel function of the second kind](http://mathworld.wolfram.com/ModifiedBesselFunctionoftheSecondKind.html).
The derivation of the final integral
is available [here](https://math.stackexchange.com/questions/2341193).
Therefore

$$
\begin{align}
p(x | y) & = \frac{\exp\left(-\sqrt{x^2 + y^2}\right)}{2y K_1\left(y\right)} ~, \\
p(x) & = \frac{2}{\pi} x K_1\left(x\right) ~.
\end{align}
$$

The longer form with arbitrary $\boldsymbol{\mu}$
and $\boldsymbol{\Lambda}$ is obtained by a change of variable.

[^mv1]: Note that the name "multivariate Laplace" defines a [different
distribution](https://en.wikipedia.org/wiki/Multivariate_Laplace_distribution).
In the current page, the name refers to out proposed generalised definition
with $q=\frac{1}{2}$.
