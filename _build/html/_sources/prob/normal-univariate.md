# Univariate Normal distribution

```{include} ../latex/preambule.md
```

:::::{admonition} {term}`Probability Density Function`
````{tab-set}
```{tab-item} Standard
$$
\Dist{\varphi}{x}
= \frac{1}{\sqrt{2\pi}}\exp\left(-\frac{x^2}{2}\right)
$$
```
```{tab-item} Mean and variance
$$
\NCond{x}{\mu, \sigma^2}
= \frac{1}{\sigma}\Dist{\varphi}{\frac{x - \mu}{\sigma}}
= \frac{1}{\sqrt{2\pi\sigma^2}}\exp\left(-\frac{\left(x - \mu\right)^2}{2\sigma^2}\right)
$$
```
```{tab-item} Mean and precision
$$
\pcond{x}{\mu, \lambda}
= \NCond{x}{\mu, \frac{1}{\lambda}}
= \sqrt{\frac{\lambda}{2\pi}}\exp\left(-\frac{\lambda\left(x - \mu\right)^2}{2}\right)
$$
```
```{tab-item} Normal conjugate
$$
\pcond{\mu}{\mu_0, n_0, \lambda}
= \NCond{\mu}{\mu_0, \frac{1}{n_0\lambda}}
= \sqrt{\frac{n_0\lambda}{2\pi}}\exp\left(-\frac{n_0\lambda\left(\mu - \mu_0\right)^2}{2}\right)
$$
```
````
:::::

:::::{admonition} Log {term}`Probability Density Function`
````{tab-set}
```{tab-item} Standard
$$
\ln\Dist{\varphi}{x}
= -\frac{1}{2}\left(x^2 + \ln 2\pi \right)
$$
```
```{tab-item} Mean and variance
$$
\ln\NCond{x}{\mu, \sigma^2}
= -\frac{1}{2}\left(\frac{\left(x - \mu\right)^2}{\sigma^2} + \ln\sigma^2 + \ln 2\pi \right)
$$
```
```{tab-item} Mean and precision
$$
\ln\pcond{x}{\mu, \lambda}
= -\frac{1}{2}\left(\lambda\left(x - \mu\right)^2 - \ln\lambda + \ln 2\pi \right)
$$
```
```{tab-item} Normal conjugate
$$
\ln\pcond{\mu}{\mu_0, n_0, \lambda}
= \ln\NCond{\mu}{\mu_0, \frac{1}{n_0\lambda}}
= -\frac{1}{2}\left(n_0\lambda\left(\mu - \mu_0\right)^2 - \ln n_0 - \ln\lambda + \ln 2\pi \right)
$$
```
````
:::::

:::::{admonition} {term}`Cumulative Distribution Function`
````{tab-set}
```{tab-item} Standard
$$
\Dist{\Phi}{x}
= \frac{1}{2}\left[1 + \operatorname{erf}\left(\frac{x}{\sqrt{2}}\right)\right]
$$
```
```{tab-item} Mean and variance
$$
\DistCond[\mathcal{N}]{F}{x}{\mu, \sigma^2}
= \Dist{\Phi}{\frac{x-\mu}{\sigma}}
= \frac{1}{2}\left[1 + \operatorname{erf}\left(\frac{x-\mu}{\sqrt{2\sigma^2}}\right)\right]
$$
```
```{tab-item} Mean and precision
$$
\DistCond{F}{x}{\mu, \lambda}
= \frac{1}{2}\left[1 + \operatorname{erf}\left(\left(x-\mu\right)\sqrt{\frac{\lambda}{2}}\right)\right]
$$
```
````
:::::

:::::{admonition} {term}`Kernel`
````{tab-set}
```{tab-item} Standard
$$
\Dist{\varphi^\ast}{x}
= \exp\left(-\frac{x^2}{2}\right)
$$
```
```{tab-item} Mean and variance
$$
\DistCond[\mathcal{N}]{f^\ast}{x}{\mu, \sigma^2}
= \exp\left(-\frac{\left(x - \mu\right)^2}{2\sigma^2}\right)
$$
```
```{tab-item} Mean and precision
$$
\DistCond{f^\ast}{x}{\mu, \lambda}
= \exp\left(-\frac{\lambda\left(x - \mu\right)^2}{2}\right)
$$
```
````
:::::

:::::{admonition} Log {term}`Kernel`
````{tab-set}
```{tab-item} Standard
$$
\ln\Dist{\varphi^\ast}{x}
= -\frac{x^2}{2}
$$
```
```{tab-item} Mean and variance
$$
\ln\DistCond[\mathcal{N}]{f^\ast}{x}{\mu, \sigma^2}
= -\frac{1}{2}\left(\frac{\left(x - \mu\right)^2}{\sigma^2}\right)
$$
```
```{tab-item} Mean and precision
$$
\ln\DistCond{f^\ast}{x}{\mu, \lambda}
= -\frac{1}{2}\left(\lambda\left(x - \mu\right)^2\right)
$$
```
````
:::::

:::::{admonition} {term}`Partition Function`
````{tab-set}
```{tab-item} Standard
$$
Z_\varphi
= \sqrt{2\pi}
$$
```
```{tab-item} Mean and variance
$$
\Dist[\mathcal{N}]{Z}{\mu, \sigma^2}
= \sqrt{2\pi\sigma^2}
$$
```
```{tab-item} Mean and precision
$$
\Dist{Z}{\mu, \lambda}
= \sqrt{\frac{2\pi}{\lambda}}
$$
```
````
:::::

:::::{admonition} Log {term}`Partition Function`
````{tab-set}
```{tab-item} Standard
$$
\ln Z_\varphi
= \frac{1}{2}\ln 2\pi
$$
```
```{tab-item} Mean and variance
$$
\ln\Dist[\mathcal{N}]{Z}{\mu, \sigma^2}
= \frac{1}{2}\left(\ln 2\pi + \ln\sigma^2\right)
$$
```
```{tab-item} Mean and precision
$$
\ln\Dist{Z}{\mu, \lambda}
= \frac{1}{2}\left(\ln 2\pi - \ln\lambda\right)
$$
```
````
:::::

:::::{admonition} {term}`Score Function`
````{tab-set}
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

:::::{admonition} {term}`Kullback-Leibler Divergence`
````{tab-set}
```{tab-item} Mean and variance
$$
\KL[\mathcal{N}]{\mu_1, \sigma_1^2}{\mu_0, \sigma_0^2}
= \frac{1}{2\sigma_0^2}\left(\mu_1 - \mu_0\right)^2
+ \frac{1}{2}\left(\frac{\sigma_1^2}{\sigma_0^2} - \ln\left(\frac{\sigma_1^2}{\sigma_0^2}\right) - 1\right)
$$
```
```{tab-item} Mean and precision
$$
\KL{\p{\mu_1, \lambda_1}}{\p{\mu_0, \lambda_0}}
= \frac{\lambda_0}{2}\left(\mu_1 - \mu_0\right)^2
+ \frac{1}{2}\left(\frac{\lambda_0}{\lambda_1} - \ln\left(\frac{\lambda_0}{\lambda_1}\right) - 1\right)
$$
```
```{tab-item} Normal conjugate
$$
\KL{\p{\mu_1, n_1, \lambda}}{\p{\mu_0, n_0, \lambda}}
= \frac{n_0\lambda}{2}\left(\mu_1 - \mu_0\right)^2
+ \frac{1}{2}\left(\frac{n_0}{n_1} - \ln\left(\frac{n_0}{n_1}\right) - 1\right)
$$
```
````
:::::

:::::{admonition} {term}`Entropy`
````{tab-set}
```{tab-item} Standard
$$
H_\varphi
= \frac{1}{2} \left(\ln 2\pi + 1\right)
$$
```
```{tab-item} Mean and variance
$$
\Dist[\mathcal{N}]{H}{\mu, \sigma^2}
= \frac{1}{2} \left(\ln 2\pi + \ln\sigma^2 + 1\right)
$$
```
```{tab-item} Mean and precision
$$
\Dist{H}{\p{\mu, \lambda}}
= \frac{1}{2} \left(\ln 2\pi - \ln\lambda + 1\right)
$$
```
````
:::::


:::::{admonition} {term}`Maximum-Likelihood`
````{tab-set}
```{tab-item} Mean and variance
$$
\begin{align}
\left(\hat{\mu}, \hat{\sigma}^2 \right)
& {}= \argmax_{\mu, \sigma^2}~\NCond{\vec{x}}{\mu, \sigma^2}
\\
& {}= \left(
    \frac{1}{N}\sum_n x_n,
    \left[\frac{1}{N}\sum_{n}x_n\right]^2 - \frac{1}{N} \sum_n x_n^2
\right)
\\
& {}= \left( \overline{x}, \overline{x}^2 - \overline{x^2} \right)
\end{align}
$$
```
```{tab-item} Mean (fixed variance)
$$
\hat{\mu}
= \argmax_{\mu}~\NCond{\vec{x}}{\mu, \sigma^2}
= \frac{1}{N}\sum_n x_n
= \overline{x}
$$
```
```{tab-item} Variance (fixed mean)
$$
\hat{\sigma}^2
= \argmax_{\sigma^2}~\NCond{\vec{x}}{\mu, \sigma^2}
= \frac{1}{N}\sum_n (x_n - \mu)^2
= \overline{x-\mu}
$$
```
````
:::::


:::::{admonition} {term}`Bayesian Update`
````{tab-set}
```{tab-item} Mean
$$
\begin{align}
\mu & \sim \N{\mu_0, \frac{1}{n_0\lambda}}
\\
x \mid \mu & \sim \N{\mu, \lambda^{-1}}
\\
\Rightarrow~ \mu \mid \vec{x}
& \sim \N{\frac{n_0 \mu_0 + n \hat{\mu}}{n_0 + n}, \frac{1}{\left(n_0+n\right)\lambda}}
\end{align}
$$
```
```{tab-item} Variance
$$
\begin{align}
\sigma^2 & \sim \NamedDist{Inv\text{-}Gamma}{\frac{n_0}{2}, \frac{n_0\sigma_0^{2}}{2}}
\\
x \mid \sigma^2 & \sim \N{\mu, \sigma^2}
\\
\Rightarrow~ \sigma^2 \mid \vec{x}
& \sim \NamedDist{Inv\text{-}Gamma}{\frac{n_0 + n}{2}, \frac{n_0\sigma_0^{2} + n\hat{\sigma}^2}{2}}
\end{align}
$$
```
```{tab-item} Precision
$$
\begin{align}
\lambda & \sim \NamedDist{Gamma}{\frac{n_0}{2}, \frac{n_0\lambda_0^{-1}}{2}}
\\
x \mid \lambda & \sim \N{\mu, \lambda^{-1}}
\\
\Rightarrow~ \lambda \mid \vec{x}
& \sim \NamedDist{Gamma}{\frac{n_0 + n}{2}, \frac{n_0\lambda_0^{-1} + n\hat{\lambda}^{-1}}{2}}
\end{align}
$$
```
```{tab-item} Mean and Variance
$$
\begin{alignat}{4}
&& \sigma^2 & \sim \NamedDist{Inv\text{-}Gamma}{\frac{n_0}{2}, \frac{n_0\sigma_0^{2}}{2}}
\\
&& \mu \mid \sigma^2 & \sim \N{\mu_0, \frac{\sigma^2}{m_0}}
\\
&& x \mid \mu, \sigma^2 & \sim \N{\mu, \sigma^2}
\\
\Rightarrow~
&& \sigma^2 \mid \vec{x}
& \sim \NamedDist{Inv\text{-}Gamma}{\frac{n_0 + n}{2}, \frac{n_0\sigma_0^{2} + n\hat{\sigma}^2}{2}}
\\
&& \mu \mid \vec{x}, \sigma^2
& \sim \N{\frac{m_0 \mu_0 + n \hat{\mu}}{m_0 + n}, \frac{\sigma^2}{m_0+n}}
\end{alignat}
$$
```
```{tab-item} Mean and Precision
$$
\begin{alignat}{4}
&& \sigma^2 & \sim \NamedDist{Gamma}{\frac{n_0}{2}, \frac{n_0}{2\lambda_0}}
\\
&& \mu \mid \lambda & \sim \N{\mu_0, \frac{1}{m_0\lambda}}
\\
&& x \mid \mu, \lambda & \sim \N{\mu, \lambda^{-1}}
\\
\Rightarrow~
&& \lambda \mid \vec{x}
& \sim \NamedDist{Gamma}{\frac{n_0 + n}{2}, \frac{n_0\lambda_0^{-1} + n\hat{\lambda}^{-1}}{2}}
\\
&& \mu \mid \vec{x}, \lambda
& \sim \N{\frac{m_0 \mu_0 + n \hat{\mu}}{m_0 + n}, \frac{1}{\lambda\left(m_0+n\right)}}
\end{alignat}
$$
```
````
:::::
