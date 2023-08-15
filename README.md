# MaxWiK: Maxima Weighted Isolation Kernel Mapping Method

[![License](https://img.shields.io/badge/License-AGPL-orange.svg)](https://github.com/tugHall/MaxWiK/blob/main/LICENSE.md)


**MaxWiK** _(**Max**ima **W**eighted-**i**solation **K**ernel mapping method)_ is a machine learning method of meta-sampling based on Isolation Kernel and Kernel mean embedding.
For more details of the method, please, be kind to read the paper:

[Iurii S. Nagornov, Approximate Bayesian Computation Based on Maxima Weighted Isolation Kernel Mapping, arXiv.2201.12745, 2022](https://doi.org/10.48550/arXiv.2201.12745)

Authors and contributor list:
---
_**Iurii (Yuri) Nagornov**_ (Maintainer, Author)

All questions and requests can be sent to iurii.nagornov@aist.go.jp or nagornov.yuri@gmail.com 

Short description
---

**Motivation:** A branching processes model yields an unevenly stochastically distributed dataset that consists of sparse and dense regions. This work addresses the problem of precisely evaluating parameters for such a model. Applying a branching processes model to an area such as cancer cell evolution faces a number of obstacles, including high dimensionality and the rare appearance of a result of interest. We take on the ambitious task of obtaining the coefficients of a model that reflects the relationship of driver gene mutations and cancer hallmarks on the basis of personal data regarding variant allele frequencies. 

**Method:** An approximate Bayesian computation method based on Isolation Kernel is developed. The method involves the transformation of row data to a Hilbert space (mapping) and the measurement of the similarity between simulated points and maxima weighted Isolation Kernel mapping related to the observation point. We also design a meta-sampling algorithm for parameter estimation that requires no gradient calculation and is dimension independent. The advantages of the proposed machine learning method are more clearly can be illustrated using multidimensional data as well as a specific branching processes model like cancer cell evolution.
 
**Package:** This software is a package named **MaxWiK** contains Approximate Bayesian Computation methods to choose a single parameter for a single observation point. 

To install, please, use the archive file 'MaxWiK_1.0.0.tar.gz':

    utils::install.packages("./MaxWiK_1.0.0.tar.gz", repos = NULL, type = "source")

To see how it works, please, be kind use the templates. To get templates, please, use command:

    copy_templates( dir = './' )   
    # dir can be any working folder where template will be copied


## Cite package MaxWiK

For publication, please, be kind to use next references related to MaxWiK software:

- [Iurii S. Nagornov, Approximate Bayesian Computation Based on Maxima Weighted Isolation Kernel Mapping, arXiv.2201.12745, 2022](https://doi.org/10.48550/arXiv.2201.12745)

- [Iurii S. Nagornov, Overfitting Problem in the Approximate Bayesian Computation Method Based on Maxima Weighted Isolation Kernel, The 36th Annual Conference of the Japanese Society for Artificial Intelligence, 2S5-IS-2c-05, 2022](https://confit.atlas.jp/guide/event/jsai2022/subject/2S5-IS-2c-05/tables?cryptoId=)

- [Package MaxWiK: https://github.com/tugHall/MaxWiK](https://github.com/tugHall/MaxWiK.git)


