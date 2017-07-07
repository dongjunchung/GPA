GPA
===

GPA is a statistical approach to prioritizing GWAS results by integrating pleiotropy information and annotation data. 'GPA' package provides computationally efficient and user friendly interface to fit GPA models and implement the hypothesis testing for the pleiotropy and the enrichment of annotation for the associated SNPs.

The stable versions of GPA package can be obtained from the following URLs:

Package source: https://github.com/dongjunchung/GPA_binary/blob/master/GPA_0.9-3.tar.gz?raw=true

Windows binary: https://github.com/dongjunchung/GPA_binary/blob/master/GPA_0.9-3.zip?raw=true

Mac OS/X binary: https://github.com/dongjunchung/GPA_binary/blob/master/GPA_0.9-3.tgz?raw=true

GPA vignette provides a good start point for the step-by-step data analysis using GPA package and it can be found at Package source: https://github.com/dongjunchung/GPA/blob/master/inst/doc/GPA-example.pdf?raw=true. Please check https://groups.google.com/d/forum/gpa-user-group for discussions and questions regarding genetic data analysis using GPA package. You can track development of GPA package at http://github.com/dongjunchung/GPA.

Usage
===========

The following two help pages provide a good start point for the genetic analysis using GPA package, including the overview of GPA package and the example command lines:

```
library(GPA)
package?GPA
class?GPA
```

Development
===========

To install the development version of GPA, it's easiest to use the 'devtools' package. Note that GPA depends on the 'Rcpp' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.

```
#install.packages("devtools")
library(devtools)
install_github("dongjunchung/GPA")
```

References
==========

Chung D\*, Yang C\*, Li C, Gelernter J, and Zhao H (2014), "GPA: A statistical approach to prioritizing GWAS results by integrating pleiotropy information and annotation data." PLoS Genetics, 10: e1004787. (\* joint first authors)
