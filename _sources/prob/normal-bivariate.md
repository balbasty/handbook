# Bivariate Normal distribution

```{include} ../latex/preambule.md
```

This is a special case of the multivariate Normal distribution in two
dimensions. However, since this distribution is quite used in image
registration, we provide specific derivations that may be useful.

Let us write the covariance matrix as

$$
\mat{\Sigma} =
\begin{pmatrix}
\sigma_x^2 & \alpha\sigma_x\sigma_y \\
\alpha\sigma_x\sigma_y  & \sigma_y^2
\end{pmatrix}
=
\begin{pmatrix}
\sigma_x & 0 \\
0  & \sigma_y
\end{pmatrix}
\times
\underset{\text{correlation matrix}}{\begin{pmatrix}
1 & \alpha \\
\alpha  & 1
\end{pmatrix}}
\times
\begin{pmatrix}
\sigma_x & 0 \\
0  & \sigma_y
\end{pmatrix}
~,
$$

where $\sigma_x \in \mathbb{R}_\ast^+$ and $\sigma_y \in \mathbb{R}_\ast^+$ are
the standard deviations along the first and second dimensions, and
$\alpha \in [0, 1]$ is the correlation coefficient. Inverting a 2x2 matrix
can be done in closed-form and we find the corresponding precision matrix

$$
\begin{align}
\mat{\Lambda} = \mat{\Sigma}^{-1}
& {}= \frac{1}{\sigma_x^2\sigma_y^2\left(1 - \alpha^2\right)}
\begin{pmatrix}
\sigma_y^2 & -\alpha\sigma_x\sigma_y \\
-\alpha\sigma_x\sigma_y  & \sigma_x^2
\end{pmatrix}
\\
\\
& {}= \frac{1}{\alpha^2 - 1}
\begin{pmatrix}
\lambda_x & 0 \\
0  & \lambda_y
\end{pmatrix}
\times
\begin{pmatrix}
1 & \alpha \\
\alpha  & 1
\end{pmatrix}
\times
\begin{pmatrix}
\lambda_x & 0 \\
0  & \lambda_y
\end{pmatrix}
~ ,
\end{align}
$$

where $\lambda_x = \sigma_x^{-2}$ and $\lambda_y = \sigma_y^{-2}$.


:::::{admonition} {term}`Probability Density Function`
````{tab-set}
```{tab-item} Standard
$$
\NCond{x, y}{\alpha}
= \frac{1}{\sqrt{2\pi\left(1 - \alpha^2\right)}}
\exp\left(-\frac{x^2 + y^2 + 2\alpha x y}{2}\right)
$$
```
```{tab-item} Mean and variance
$$
\NCond{x, y}{\mu_x, \mu_y, \sigma_x, \sigma_y, \alpha}
= \frac{1}{\sqrt{2\pi\sigma_x^2\sigma_y^2\left(1 - \alpha^2\right)}}
\exp\left(-\frac{\hat{x}^2 + \hat{y}^2 + 2\alpha \hat{x} \hat{y}}{2}\right)
$$
where $\hat{x} = \frac{x - \mu_x}{\sigma_x}$ and
$\hat{y} = \frac{y - \mu_y}{\sigma_y}$.
```
```{tab-item} Mean and precision
$$
\pcond{x, y}{\mu_x, \mu_y, \lambda_x, \lambda_y, \alpha}
= \NCond{x, y}{\mu_x, \mu_y, \lambda_x^{-1}, \lambda_y^{-1}, \alpha}
= \sqrt{\frac{\lambda_x\lambda_y}{2\pi\left(1 - \alpha^2\right)}}
\exp\left(-\frac{\hat{x}^2 + \hat{y}^2 + 2\alpha \hat{x} \hat{y}}{2}\right)
$$
where $\hat{x} = \lambda_x\left(x - \mu_x\right)$ and
$\hat{y} = \lambda_y\left(x - \mu_y\right)$.
```
:::::

:::::{admonition} Log {term}`Probability Density Function`
````{tab-set}
```{tab-item} Standard
$$
\ln\NCond{x, y}{\alpha}
= -\frac{1}{2}\left(
    x^2 + y^2 + 2\alpha x y
    + \ln\left(1 - \alpha^2\right) + \ln 2\pi
\right)
$$
```
```{tab-item} Mean and variance
$$
\ln\NCond{x, y}{\mu_x, \mu_y, \sigma_x, \sigma_y, \alpha}
= -\frac{1}{2}\left(
    \hat{x}^2 + \hat{y}^2 + 2\alpha \hat{x} \hat{y}
    + \ln\sigma_x^2 + \ln\sigma_y^2
    + \ln\left(1 - \alpha^2\right) + \ln 2\pi
\right)
$$
where $\hat{x} = \frac{x - \mu_x}{\sigma_x}$ and
$\hat{y} = \frac{y - \mu_y}{\sigma_y}$.
```
```{tab-item} Mean and precision
$$
\ln\pcond{x, y}{\mu_x, \mu_y, \lambda_x, \lambda_y, \alpha}
= -\frac{1}{2}\left(
    \hat{x}^2 + \hat{y}^2 + 2\alpha \hat{x} \hat{y}
    - \ln\lambda_x - \ln\lambda_y
    + \ln\left(1 - \alpha^2\right) + \ln 2\pi
\right)
$$
where $\hat{x} = \lambda_x\left(x - \mu_x\right)$ and
$\hat{y} = \lambda_y\left(x - \mu_y\right)$.
```
:::::

:::::{admonition} {term}`Kernel`
````{tab-set}
```{tab-item} Standard
$$
\DistCond[\mathcal{N}]{f^\ast}{x, y}{\alpha}
= \exp\left(-\frac{x^2 + y^2 + 2\alpha x y}{2}\right)
$$
```
```{tab-item} Mean and variance
$$
\DistCond[\mathcal{N}]{f^\ast}{x, y}{\mu_x, \mu_y, \sigma_x^2, \sigma_y^2, \alpha}
= \exp\left(-\frac{\hat{x}^2 + \hat{y}^2 + 2\alpha \hat{x} \hat{y}}{2}\right)
$$
where $\hat{x} = \frac{x - \mu_x}{\sigma_x}$ and
$\hat{y} = \frac{y - \mu_y}{\sigma_y}$.
```
```{tab-item} Mean and precision
$$
\DistCond{f^\ast}{x}{\mu_x, \mu_y, \lambda_x, \lambda_y, \alpha}
= \exp\left(-\frac{\hat{x}^2 + \hat{y}^2 + 2\alpha \hat{x} \hat{y}}{2}\right)
$$
where $\hat{x} = \lambda_x\left(x - \mu_x\right)$ and
$\hat{y} = \lambda_y\left(x - \mu_y\right)$.
```
````
:::::

:::::{admonition} Log {term}`Kernel`
````{tab-set}
```{tab-item} Standard
$$
\ln\DistCond[\mathcal{N}]{f^\ast}{x, y}{\alpha}
= -\frac{x^2 + y^2 + 2\alpha x y}{2}
$$
```
```{tab-item} Mean and variance
$$
\ln\DistCond[\mathcal{N}]{f^\ast}{x, y}{\mu_x, \mu_y, \sigma_x^2, \sigma_y^2, \alpha}
= -\frac{\hat{x}^2 + \hat{y}^2 + 2\alpha \hat{x} \hat{y}}{2}
$$
where $\hat{x} = \frac{x - \mu_x}{\sigma_x}$ and
$\hat{y} = \frac{y - \mu_y}{\sigma_y}$.
```
```{tab-item} Mean and precision
$$
\ln\DistCond{f^\ast}{x}{\mu_x, \mu_y, \lambda_x, \lambda_y, \alpha}
= -\frac{\hat{x}^2 + \hat{y}^2 + 2\alpha \hat{x} \hat{y}}{2}
$$
where $\hat{x} = \lambda_x\left(x - \mu_x\right)$ and
$\hat{y} = \lambda_y\left(x - \mu_y\right)$.
```
````
:::::

:::::{admonition} {term}`Partition Function`
````{tab-set}
```{tab-item} Standard
$$
\Dist[\mathcal{N}]{Z}{\alpha}
= \sqrt{2\pi\left(1-\alpha^2\right)}
$$
```
```{tab-item} Variance
$$
\Dist[\mathcal{N}]{Z}{\sigma_x^2, \sigma_y^2, \alpha}
= \sqrt{2\pi\sigma_x^2\sigma_y^2\left(1-\alpha^2\right)}
$$
```
```{tab-item} Precision
$$
\Dist{Z}{\lambda_x, \lambda_y, \alpha}
= \sqrt{\frac{2\pi\left(1-\alpha^2\right)}{\lambda_x\lambda_y}}
$$
```
````
:::::

:::::{admonition} Log {term}`Partition Function`
````{tab-set}
```{tab-item} Standard
$$
\ln\Dist[\mathcal{N}]{Z}{\alpha}
= \frac{1}{2}\left(\ln 2\pi + \ln\left(1-\alpha^2\right) \right)
$$
```
```{tab-item} Mean and variance
$$
\ln\Dist[\mathcal{N}]{Z}{\sigma_x^2, \sigma_y^2, \alpha}
= \frac{1}{2}\left(\ln 2\pi + \ln\sigma_x^2 + \ln\sigma_y^2 + \ln\left(1-\alpha^2\right) \right)
$$
```
```{tab-item} Mean and precision
$$
\ln\Dist{Z}{\lambda_x, \lambda_y, \alpha}
= \frac{1}{2}\left(\ln 2\pi - \ln\lambda_x - \ln\lambda_y + \ln\left(1-\alpha^2\right) \right)
$$
```
````
:::::

:::::{admonition} {term}`Score Function`
````{tab-set}
```{tab-item} Standard
$$
s\left(\alpha ~;~ x\right)
= \frac{\alpha}{1-\alpha^2} - xy
$$
```
```{tab-item} Mean and variance
$$
\begin{align}
s\left(\mu ~;~ \sigma^2, x\right)
& = \frac{x - \mu}{\sigma^2}
\\
s\left(\sigma^2 ~;~ \mu, x\right)
& = \frac{1}{2\sigma^2}\left(\frac{\left(x - \mu\right)^2}{\sigma^2} - 1\right)
\end{align}
$$
```
```{tab-item} Mean and precision
$$
\begin{align}
s\left(\mu ~;~ \lambda, x\right)
& = \lambda\left(x - \mu\right)
\\
s\left(\lambda ~;~ \mu, x\right)
& = \frac{1}{2}\left(\frac{1}{\lambda} - \left(x - \mu\right)^2\right)
\end{align}
$$
```
````
:::::
