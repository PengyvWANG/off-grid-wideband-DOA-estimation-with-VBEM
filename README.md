# off-grid-wideband-DOA-estimation-with-VBEM
Implementation of journal paper entitled 'An off-grid wideband DOA estimation method with the variational Bayes expectation-maximization framework' for DOA estimation [Signal Processing]

**Notation:** Codes for comparison methods are not available. Please remove them from `main.m` manually.

## Citation
```
@article{WANG2022108423,
title = {An off-grid wideband DOA estimation method with the variational Bayes expectation-maximization framework},
journal = {Signal Processing},
volume = {193},
pages = {108423},
year = {2022},
issn = {0165-1684},
doi = {https://doi.org/10.1016/j.sigpro.2021.108423},
url = {https://www.sciencedirect.com/science/article/pii/S0165168421004606},
author = {Pengyu Wang and Huichao Yang and Zhongfu Ye},
keywords = {DOA Estimation, Wideband applications, Sparse Bayesian learning, Expectation-maximization algorithm, Variational approach},
abstract = {Wideband direction-of-arrival (DOA) estimation has important applications in many fields. To improve the performance of wideband DOA estimation, an off-grid approach with the variational Bayes Expectation-Maximization (VBEM) framework is proposed. In this work, a hierarchical Bernoulli-Gaussian model with several latent variables and parameters is formulated to describe the array observation. A latent indicator vector is assigned to describe the existence of the potential sources as well as promote sparsity in the angle domain. The VBEM framework is utilized to solve the latent variables and parameters. By modeling the existence and non-existence of potential sources separately and utilizing the wideband data jointly, our method performs better than the existing methods. The proposed method does not require any prior knowledge of the number of sources nor the initial DOAs, and it can solve both coherent and incoherent wideband sources. Simulation results show that the proposed method is superior to the existing wideband DOA estimation methods with few snapshots in frequency.}
}
```
