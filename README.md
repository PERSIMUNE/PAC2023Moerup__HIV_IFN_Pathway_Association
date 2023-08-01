# Genes SKATO Analysis 
##### Code repository used for the publication of HIV type 1 IFN genes association with HIV Viral Load

Initial Development by:
- Cameron MacPherson (cameron@biostone.consulting)

Modification of code and dockerisation by:
- Preston Leung (preston.yui.sum.leung@regionh.dk)

Author of the study:
- Sara Bohnstedt Moerup (sara.bohnstedt.moerup@regionh.dk)

## Citation:
Please cite this paper when using this tool.

## Tool Description
Originally developed for testing SNPs within a single gene and its association with a specific, this tool has been modified to include SNPs from different genes, allowing gene sets to be tested. Genes SKATO Analysis are a set of RScripts put into a docker such that it can be downloaded and executed on the user's computer without the need to install required packages.  

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
| Docker Parameters | Description |
|----|----|
| --rm | Automatically remove the container when it exits |
| -v | Bind mount a volume |


| Required Parameters | Description |
|----|----|
| -o | Filename to write the output to. |
| -p | Filename of the input file. |

**Note: ** In the above example command, the data required to run the tool is assumed be in `/dockerWorkDir/TORUN/` folder or whatever name used to in place of `TORUN`. The output can be accessed in `/TRUE_PATH/TO/DATA` when it is completed.

### Expected file format of the input file
In the input file used for the `-p` parameter, the `.csv` (comma-separated) is expected to store a dataframe columns about the data being used and the location of data.

An example of the parameter file would look like this:
|entrez|window|pidfile|confounder.1|confounder.2|confounder.3|cohort|outcome|ignore|affy_file|pheno_file|geno_file
|----|----|----|----|----|----|----|----|----|----|----|----|
|entrez ID| +/- window size from gene region|/PATH/TO/PIDFILE|Covariate 1|Covariate 2|Covariate 3|COHORT-ID|Outcome Variable|1 to ignore else 0|/PATH/TO/AFFY|/PATH/TO/PHENO|/PATH/TO/GENO|

*Note:* If more covariates are to be needed for your dataset, simply add a column `confounder.N` where `N` is an integer.


|File Types | Description |
|----|----|
|pidfile| This is a list (`.txt` is fine) of subject/patient IDs that are to be included in the run. |
|affy_file| |
|pheno_file| |
|geno_file| |


[![CHIP](https://chip.dk/Portals/0/CHIP_new.png?ver=2020-10-01-104734-463)](https://chip.dk)
[![INSIGHT](https://chip.dk/portals/0/files/INSIGHT/INSIGHT-logo.png?ver=2020-06-22-123834-000)](http://insight.ccbr.umn.edu)
