GPA
===

GPA (**G**enetic analysis incorporating **P**leiotropy and **A**nnotation) is a statistical approach to prioritizing GWAS results by integrating pleiotropy information and annotation data, proposed in Chung et al. (2014). 'GPA' package provides computationally efficient and user friendly interface to fit the GPA models and implement the hypothesis testing for the pleiotropy and the enrichment of annotation for the associated SNPs. [The 'GPA' vignette](https://github.com/dongjunchung/GPA/blob/master/inst/doc/GPA-example.pdf?raw=true) provides a good start point for the step-by-step data analysis using 'GPA' package. Please check [our GPA Google Group](https://groups.google.com/d/forum/gpa-user-group) for discussions and questions regarding genetic data analysis using 'GPA' package. The following two help pages provide a good start point for the genetic analysis using the 'GPA' package, including the overview of 'GPA' package and the example command lines:

```
library(GPA)
package?GPA
class?GPA
```

ShinyGPA
========

ShinyGPA is an interactive and flexible visualization framework to investigate the pleiotropic architecture using GWAS results, proposed in Kortemeier et al. (2017). The following help page provides the overview of ShinyGPA and the example command lines:

```
library(GPA)
?shinyGPA
```

Installation
============ 

The stable versions of 'GPA' package can be obtained from the following URLs:

Package source: [https://github.com/dongjunchung/GPA_binary/blob/master/GPA_1.1-0.tar.gz?raw=true](https://github.com/dongjunchung/GPA_binary/blob/master/GPA_1.1-0.tar.gz?raw=true)

Windows binary: [https://github.com/dongjunchung/GPA_binary/blob/master/GPA_1.1-0.zip?raw=true](https://github.com/dongjunchung/GPA_binary/blob/master/GPA_1.1-0.zip?raw=true)

Mac OS/X binary: [https://github.com/dongjunchung/GPA_binary/blob/master/GPA_1.1-0.tgz?raw=true](https://github.com/dongjunchung/GPA_binary/blob/master/GPA_1.1-0.tgz?raw=true)

To install the developmental versions of 'GPA' package, it's easiest to use the 'devtools' package. Note that the 'GPA' package depends on the 'Rcpp' package, which also requires appropriate setting of *Rtools* and *Xcode* for Windows and Mac OS/X, respectively.

```
#install.packages("devtools")
library(devtools)
install_github("dongjunchung/GPA")
```

References
==========

Chung D\*, Yang C\*, Li C, Gelernter J, and Zhao H (2014), "GPA: A statistical approach to prioritizing GWAS results by integrating pleiotropy information and annotation data," *PLoS Genetics*, 10: e1004787. (\* joint first authors)

Kortemeier E, Ramos PS, Hunt KJ, Kim HJ, Hardiman G, and Chung D (2017), "ShinyGPA: An interactive and dynamic visualization toolkit for genetic studies."
