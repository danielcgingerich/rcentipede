# CENTIPEDE (R implementation)
### Introduction
This is an R implementation of the CENTIPEDE model first proposed by R Pique Regi et al. 2011, Accurate inference of transcription factor binding from DNA sequence and chromatin accessibility data. Though the model was originally implemented for DNAse-seq datasets, we extend its applications to single nuclei ATAC-seq data. 

### Model overview
Model notation – we have matrix $X \in \mathbb{R}^{N \times 200} $, where each row corresponds to the transcription factor binding site (TFBS) in question, and each column corresponds to the base pair (bp) position flanking the TFBS, from $-100$ (upstream) to $100$ (downstream). The values of $X$ indicate the number of cutsites detected at position $j$. We denote $X_i$ as row $i$ of $X$. We denote $R_i=\sum_j X_{ij}$. Additionally, we introduce a latent variable, $Z_i$, to denote whether TFBS $i$ is bound ($Z_i = 1$) or not ($Z_i = 0$).

Marginalizing out $Z$, we have: 

$$ P(\text{Data}) = \pi P(\text{Data} | Z =1 ) + (1-\pi) P(\text{Data} | Z = 0 )$$ 

Where $Data$ is our Tn5 cutsite data around each motif site, and $\pi \in [0, 1]$ is the mixture proportion. The Tn5 distribution around a motif is parameterized as a multinomial distribution conditional on the number of Tn5 cutsites observed within the window and the binding state (bound/unbound), and the total number of cutsites is parametrized as a negative binomial:

$$ P(X_i | Z_i=1, R_i, p, \mu_1, \phi_1) = \text{Multinomial} (X_i | R_i, p) \cdot  \text{NegBin} (R_i | \mu_1, \phi_1) $$

$$ P(X_i | Z_i=0, R_i, \lambda, \mu_0, \phi_0) = \text{Multinomial} (X_i | R_i, \lambda) \cdot  \text{NegBin} (R_i | \mu_0, \phi_0) $$

Where $p$ is a vector of multinomial probabilities, $\lambda = (\frac{1}{200}, ..., \frac{1}{200})$ is a vector of equal probabilities parameterizing a uniform multinomial distribution, $\mu_1$ and $\phi_1$ are the mean and dispersion for the bound TF negative binomial, and $\mu_0$ and $\phi_0$ are the mean and dispersion of the unbound negative binomial.

The complete likelihood (if $Z$ were known) is: 

$$
\begin{aligned}
P(X | Z, R) = & \prod_{i=1}^N ( \pi \text{Multinomial} (X_i | R_i , p) \cdot \text{NegBin} (R_i | \mu_1, \phi_1 ) )^{Z_i}  \\
& \cdot ( (1-\pi) \text{Multinomial}(X_i | R_i, \lambda ) \cdot \text{NegBin} (R_i | \mu_0, \phi_0))^{1-Z_i}
\end{aligned}
$$

And the corresponding log-likelihood: 

$$
\begin{aligned}
    p(X | Z, R, p, \lambda, \mu_1, \mu_0, \phi_1, \phi_0) = & \sum_{i=1}^N Z_i \log \pi + Z_i X_i \log p + Z_i \big( \log \Gamma (R_i + \phi_1) \\
    & - \log \Gamma( \phi_1 ) - \log \Gamma (R_i + 1 ) + \phi_1 \log (\phi_1) \\
    & - \phi_1 \log (\phi_1 + \mu_1 ) + R_i \log ( \mu_1 ) - R_i \log (\phi_1 + \mu_1 ) \big) \\
    & + (1 - Z_i) \log (1 - \pi ) + (1-Z_i) X_i \log (\lambda) \\
    & + (1 - Z_i) \big( \log \Gamma( R_i + \phi_0 ) - \log \Gamma (\phi_0 ) \\
    & - \log \Gamma (R_i + 1) + \phi_0 \log ( \phi_0 ) - \phi_0 \log (\phi_0 + \mu_0) \\ 
    &+ R_i \log (\mu_0) - R_i \log (\phi_0 + \mu_0) \big) 
\end{aligned}
$$

### Estimating the model parameters
As $Z_i$ is an unobserved latent variable, CENTIPEDE uses expectation maximization to estimate the parameters. 

*E-step* &mdash; First, we calculate the expectation of $Z$ given the initial parameter estimates: 

$\mathbb{E}(Z_i | X_i, R_i, p, \mu_1, \mu_0, \phi_1, \phi_0 )= $

$$ \frac{
\pi \text{Multinomial}(X_i | p, R_i) \cdot \text{NegBin} (R_i | \mu_1, \phi_1) 
}{
\pi \text{Multinomial}(X_i | p, R_i) \cdot \text{NegBin} (R_i | \mu_1, \phi_1 + (1-\pi) \text{Multinomial}(X_i | R_i, \lambda) \cdot \text{NegBin} (R_i | \mu_0, \phi_0)
}$$

*M-step* &mdash; Replacing $Z_i$ from the complete log likelihood with its expectation, the parameter $p$, representing the probabilities of the multinomial distribution, is updated as: 

$$p^{(t+1)} = \frac{
\sum_i X_{ij} \mathbb{E}(Z_i^{(t)} | X_i, \theta^{(t)})
}{
\sum_i \sum_j X_{ij} \mathbb{E} (Z_i^{(t)} | X_i, \theta^{(t)} )
}$$

Where $t$ is the current iteration and $\theta$ is shorthand for the current parameters. 

Unfortunately, the parameters of the bound and unbound negative binomial distributions, $\mu_1$, $\phi_1$, $\mu_0$, and $\phi_0$ have no closed form solutions. However, optimizing the associated log likelihood components is equivalent to a weighted generalized linear model, and are thus updated using the function ```glm.nb``` from the R package MASS, with the weights set to $\mathbb{E}(Z_i | X_i, \theta^{(t)} )$ for the bound case, and $1 - \mathbb{E}(Z_i | X_i, \theta^{(t)} )$ for the unbound case.

The quantity $\pi$ is parametrized by a sigmoid function, and incorporates prior information: 

$$ \pi (\beta) = \frac{1}{1 + e^{-W \beta}}$$

Where W is a $\mathbb{R}^{N \times m}$ column matrix containing the motif PWM scores and any additional prior information, such as the motif conservation scores or output from predictive machine learning models. $\beta$ is optimized by the R function optim. We iterate the CENTIPEDE EM algorithm until the negative log likelihood decreases by less than 1e-3. As TF-footprints have been known to have TF-specific signatures, we run the CENTIPEDE algorithm on each TF individually within each cell subtype. For downstream analysis, we select TFs with $P(\text{bound} | \text{Data} ) = \mathbb{E} (Z_i | \theta, X, R) \geq 0.95 $.
