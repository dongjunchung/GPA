GPA
===

GPA is a statistical approach to prioritizing GWAS results by integrating pleiotropy information and annotation data. 'GPA' package provides computationally efficient and user friendly interface to fit GPA models and implement the hypothesis testing for the pleiotropy and the enrichment of annotation for the associated SNPs.

The stable versions of GPA package can be obtained from the following URLs:

Package source: https://github.com/dongjunchung/GPA_binary/blob/master/GPA_0.9-0.tar.gz?raw=true

Windows binary: https://github.com/dongjunchung/GPA_binary/blob/master/GPA_0.9-0.zip?raw=true

Mac OS/X binary: https://github.com/dongjunchung/GPA_binary/blob/master/GPA_0.9-0.tgz?raw=true


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

Chung D\*, Yang C\*, Li C, Gelernter J, and Zhao H (2014), "GPA: A statistical approach to prioritizing GWAS results by integrating pleiotropy information and annotation data." (\* joint first authors; submitted)
