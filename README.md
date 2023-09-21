# Genes SKATO Analysis 
##### Code repository used for the publication of HIV type 1 IFN genes association with HIV Viral Load
Author of the study:
- Sara Bohnstedt Moerup (sara.bohnstedt.moerup@regionh.dk)

Initial Code Development by:
- Cameron MacPherson (cameron@biostone.consulting)

Modification of code and dockerisation by:
- Preston Leung (preston.yui.sum.leung@regionh.dk)
 
Additional contacts:
- Daniel Murray (daniel.dawson.murray@regionh.dk)
- Joanne Reekie (joanne.reekie@regionh.dk)

## Citation:
Please cite this paper when using this tool (To be updated).

## Tool Description
Originally developed for testing SNPs within a single gene and its association with a specific, this tool has been modified to include SNPs from different genes, allowing gene sets to be tested. Genes SKATO Analysis are a set of RScripts put into a docker such that it can be downloaded and executed on the user's computer without the need to install required packages.  

## Downloading docker image
The image that runs the analysis can be retrieved via:
```
docker pull pdawgzgg/gene_analysis_skato
```


## Tool parameters and description
### Run command
```
docker run \
--rm \
-v /TRUE_PATH/TO/DATA:/dockerWorkDir/TORUN \
pdawgzgg/gene_analysis_skato:v1.0 \
-p /dockerWorkDir/TORUN/Analysis.parameters.csv \
-o /dockerWorkDir/TORUN/OUTPUT.csv
```
**Note:** The memory requirement could be quite intensive when dealing with large SNP files. Could potentially use upward of **36GB** of memory.

| Docker Parameters | Description |
|----|----|
| --rm | Automatically remove the container when it exits |
| -v | Bind mount a volume |


| Required Parameters | Description |
|----|----|
| -o | Filename to write the output to. |
| -p | Filename of the input file. |

**Note:** In the above example command, the data required to run the tool is assumed be in `/dockerWorkDir/TORUN/` folder or whatever name used to in place of `TORUN`. The output can be accessed in `/TRUE_PATH/TO/DATA` when it is completed.

### Expected file format of the input file
In the input file used for the `-p` parameter, the `.csv` (comma-separated) is expected to store a dataframe columns about the data being used and the location of data.

An example of the parameter file would look like this:
|entrez|window|pidfile|confounder.1|confounder.2|confounder.3|cohort|outcome|ignore|affy_file|pheno_file|geno_file
|----|----|----|----|----|----|----|----|----|----|----|----|
|entrez ID| +/- window size from gene region|/PATH/TO/PIDFILE|Covariate 1|Covariate 2|Covariate 3|COHORT-ID|Outcome Variable|1 to ignore else 0|/PATH/TO/AFFY|/PATH/TO/PHENO|/PATH/TO/GENO|

**Note:** If more covariates are to be needed for your dataset, simply add a column `confounder.N` where `N` is an integer.


|File Types | Description |
|----|----|
|pidfile| This is a list (`.txt` is fine) of subject/patient IDs that are to be included in the run. IDs are matched against IDs in **pheno_file** and **geno_file**. |
|affy_file| A dataframe-like file storing SNP information. Contains 4 columns: AffyID, rsID, chromosome and position stored in `.RDS` format. The AffyID column stores the SNP probe id name. This is the name that details the SNP count for each patient/subject. It **may or may not** be identical to rsID. |
|pheno_file| A dataframe-like file storing the clinical data and/or phenotypes (columns) of each subject/patient (rows) stored in `.RDS` format. |
|geno_file| A dataframe-like file storing subject/patient id (rows) by SNPs (columns). The number 1 or 2 represents the copies of SNP present and zero for absence for each subject/patient. Stored in `.RDS` format. The column names for the SNPs should be same as the AffyIDs in **affy_file**. |

`geno_file` and `affy_file` can be extracted from plink files using `plink` tool. Briefly, `geno_file` is a modified form of the `.raw` file. In the `.raw` file, **IID** is used distinguish subjects and is implemented as the rownames of the dataframe in the `.RDS` file. The columns **FID, PAT, MAT, SEX and PHENOTYPE** are not used. The following columns should be the SNP ids used by the SNP array (aka Affymetrix SNP IDs or AffyID). For `affy_file`, this can be extracted from the `.bim` file. The modification just adds an extra column for rsIDs, where the rsID is associated to the matching AffyID. For how to generate `.RDS` check out [this blog](https://www.r-bloggers.com/2016/12/remember-to-use-the-rds-format/).

[![CHIP](https://chip.dk/Portals/0/CHIP_new.png?ver=2020-10-01-104734-463)](https://chip.dk)

[![INSIGHT](https://chip.dk/portals/0/files/INSIGHT/INSIGHT-logo.png?ver=2020-06-22-123834-000)](http://insight.ccbr.umn.edu)

----

##### R Packages used:
  - Arora S, Morgan M, Carlson M, Pagès H (2022). _GenomeInfoDb: Utilities for manipulating chromosome names, including modifying them to follow a particular naming style_. R package version 1.30.1, <URL: https://bioconductor.org/packages/GenomeInfoDb>.
  - Bache S, Wickham H (2022). _magrittr: A Forward-Pipe Operator for R_. R package version 2.0.3, <URL: https://CRAN.R-project.org/package=magrittr>.
  - Bates D, Maechler M, Jagan M (2023). _Matrix: Sparse and Dense Matrix Classes and Methods_. R package version 1.5-4, <URL: https://CRAN.R-project.org/package=Matrix>.
  - Bengtsson H (2003). “The R.oo package - Object-Oriented Programming with References Using Standard R Code.” In Hornik K, Leisch F, Zeileis A (eds.), _Proceedings of the 3rd International Workshop on Distributed Statistical Computing (DSC 2003)_. <URL: https://www.r-project.org/conferences/DSC-2003/Proceedings/Bengtsson.pdf>.
  - Bengtsson H (2003). “The R.oo package - Object-Oriented Programming with References Using Standard R Code.” In Hornik K, Leisch F, Zeileis A (eds.), _Proceedings of the 3rd International Workshop on Distributed Statistical Computing (DSC 2003)_. <URL: https://www.r-project.org/conferences/DSC-2003/Proceedings/Bengtsson.pdf>.
  - Bengtsson H (2022). _R.utils: Various Programming Utilities_. R package version 2.12.2, <URL: https://CRAN.R-project.org/package=R.utils>.
  - Dey R, Lee S (2020). _SPAtest: Score Test and Meta-Analysis Based on Saddlepoint Approximation_. R package version 3.1.2, <URL: https://CRAN.R-project.org/package=SPAtest>.
  - Huber W, Carey VJ, Gentleman R, Anders S, Carlson M, Carvalho BS, Bravo HC, Davis S, Gatto L, Girke T, Gottardo R, Hahne F, Hansen KD, Irizarry RA, Lawrence M, Love MI, MacDonald J, Obenchain V, Ole's AK, Pag`es H, Reyes A, Shannon P, Smyth GK, Tenenbaum D, Waldron L, Morgan M (2015). “Orchestrating high-throughput genomic analysis with Bioconductor.” _Nature Methods_, *12*(2), 115-121. <URL: http://www.nature.com/nmeth/journal/v12/n2/full/nmeth.3252.html>.
  - Huber, W., Carey, J. V, Gentleman, R., Anders, S., Carlson, M., Carvalho, S. B, Bravo, C. H, Davis, S., Gatto, L., Girke, T., Gottardo, R., Hahne, F., Hansen, D. K, Irizarry, A. R, Lawrence, M., Love, I. M, MacDonald, J., Obenchain, V., Ole's, K. A, Pag`es, H., Reyes, A., Shannon, P., Smyth, K. G, Tenenbaum, D., Waldron, L., Morgan, M. (2015). “Orchestrating high-throughput genomic analysis with Bioconductor.” _Nature Methods_, *12*(2), 115-121. <URL: http://www.nature.com/nmeth/journal/v12/n2/full/nmeth.3252.html>.
  - Lawrence M, Huber W, Pagès H, Aboyoun P, Carlson M, Gentleman R, Morgan M, Carey V (2013). “Software for Computing and Annotating Genomic Ranges.” _PLoS Computational Biology_, *9*. doi: 10.1371/journal.pcbi.1003118 (URL: https://doi.org/10.1371/journal.pcbi.1003118), <URL: http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003118>.
  - Lawrence M, Huber W, Pagès H, Aboyoun P, Carlson M, Gentleman R, Morgan M, Carey V (2013). “Software for Computing and Annotating Genomic Ranges.” _PLoS Computational Biology_, *9*. doi: 10.1371/journal.pcbi.1003118 (URL: https://doi.org/10.1371/journal.pcbi.1003118), <URL: http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003118>.
  - Lawrence M, Huber W, Pagès H, Aboyoun P, Carlson M, Gentleman R, Morgan M, Carey V (2013). “Software for Computing and Annotating Genomic Ranges.” _PLoS Computational Biology_, *9*. doi: 10.1371/journal.pcbi.1003118 (URL: https://doi.org/10.1371/journal.pcbi.1003118), <URL: http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003118>.
  - Lee S, Zhao Z, Miropolsky wcfL, Wu M (2023). _SKAT: SNP-Set (Sequence) Kernel Association Test_. R package version 2.2.5, <URL: https://CRAN.R-project.org/package=SKAT>.
  - Makowski D, Lüdecke D, Patil I, Thériault R, Ben-Shachar M, Wiernik B (2023). “Automated Results Reporting as a Practical Tool to Improve Reproducibility and Methodological Best Practices Adoption.” _CRAN_. <URL: https://easystats.github.io/report/>.
  - Mark A, Thompson R, Afrasiabi C, Wu C (2021). _mygene: Access MyGene.Info_ services_. R package version 1.30.0.
  - Morgan M (2023). _BiocManager: Access the Bioconductor Project Package Repository_. R package version 1.30.20, <URL: https://CRAN.R-project.org/package=BiocManager>.
  - Müller K, Wickham H (2023). _tibble: Simple Data Frames_. R package version 3.2.1, <URL: https://CRAN.R-project.org/package=tibble>.
  - Pagès H, Carlson M, Falcon S, Li N (2021). _AnnotationDbi: Manipulation of SQLite-based annotations in Bioconductor_. R package version 1.56.2, <URL: https://bioconductor.org/packages/AnnotationDbi>.
  - Pagès H, Lawrence M, Aboyoun P (2022). _S4Vectors: Foundation of vector-like and list-like containers in Bioconductor_. R package version 0.32.4, <URL: https://bioconductor.org/packages/S4Vectors>.
  - Qiu Y, Mei J (2022). _RSpectra: Solvers for Large-Scale Eigenvalue and SVD Problems_. R package version 0.16-1, <URL: https://CRAN.R-project.org/package=RSpectra>.
  - R Core Team (2021). _R: A Language and Environment for Statistical Computing_. R Foundation for Statistical Computing, Vienna, Austria. <URL: https://www.R-project.org/>.
  - Robinson D, Hayes A, Couch S (2023). _broom: Convert Statistical Objects into Tidy Tibbles_. R package version 1.0.4, <URL: https://CRAN.R-project.org/package=broom>.
  - Venables WN, Ripley BD (2002). _Modern Applied Statistics with S_, Fourth edition. Springer, New York. ISBN 0-387-95457-0, <URL: https://www.stats.ox.ac.uk/pub/MASS4/>.
  - von Jouanne-Diedrich H (2017). _OneR: One Rule Machine Learning Classification Algorithm with Enhancements_. R package version 2.2, <URL: https://CRAN.R-project.org/package=OneR>.
  - Wickham H, Bryan J, Barrett M (2022). _usethis: Automate Package and Project Setup_. R package version 2.1.6, <URL: https://CRAN.R-project.org/package=usethis>.
  - Wickham H, François R, Henry L, Müller K, Vaughan D (2023). _dplyr: A Grammar of Data Manipulation_. R package version 1.1.2, <URL: https://CRAN.R-project.org/package=dplyr>.
  - Wickham H, Henry L (2023). _purrr: Functional Programming Tools_. R package version 1.0.1, <URL: https://CRAN.R-project.org/package=purrr>.
  - Wickham H, Hester J, Chang W, Bryan J (2022). _devtools: Tools to Make Developing R Packages Easier_. R package version 2.4.5, <URL: https://CRAN.R-project.org/package=devtools>.


----
##### Special Remarks
_I would like to give some special acknowledgement to Sara Bohnstedt Moerup for her effort and perseverance in learning and adapting to running analysis programs. As a clinician-scientist, it is not easy to jump into bioinformatics, especially when needing to deal with backend coding. Recognition of these efforts when there is a clear need to bridge the skill gap is often understated. So here, I wish express my appreciation and recognition of Sara's efforts to understand the messy side of bioinformatics in this project._ 

-Preston

----
This README file is written using [Dillinger](dillinger.io).
