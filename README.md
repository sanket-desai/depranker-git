## DepRanker: a gene impact score calculator (prioritization method) for RNAi / CRISPR screen results

Developed by Dutt lab

Version v0.1.0

### Features
**DepRanker works downstream to the shRNA analysis pipelines like edgeR and accomodates information from gene prioritization tools like ROAST and DEMETER. It integrates gene expression, copy number variation and methylation data to rank / prioritize genes from pooled screen data.**

DepRanker accepts following input files:
-	Toptags file from edgeR shRNA analysis
-	ROAST output file
-	DEMETER gene score file
-	Expression file
-	Copy number variation file
-	Methylation file 


### Availability and Implementation
The open source DepRanker package is implemented in Python and the source code is available from [https://github.com/sanket-desai/depranker-git](https://github.com/sanket-desai/depranker-git).

## Installation
DepRanker package depends on <i>[Numpy](https://numpy.org/)</i> python package may be installed using <i>[pip](https://pypi.org/project/pip/)</i> as follows:

```
python -m pip install --user numpy
```
Installation of DepRanker package can be performed by following standard installation process. Go to the package directory (where setup.py file is placed) and run the following command:
```
sudo python setup.py
```
User will be indicated upon the successfull installation of the package on the terminal.
## Usage
```
usage: depranker.py [-h] [-roast ROAST_RESULT_FILE] [-demeter DEMETER_SCORE_FILE]
                    [-toptags EDGER_TOPTAGS_FILE] [-exprs EXPRESSION_FILE]
                    [-cnv COPY_NUMBER_VARIATION_FILE] [-meth METHYLATION_BETA_FILES] [-out OUTPUT_FILE]

DepRanker: A gene impact score calculator (prioritization method) for RNAi / CRISPR screen results
Developed by Dutt lab
Version v0.1.0

optional arguments:
  -h, --help            show this help message and exit
  -roast ROAST_RESULT_FILE
                        Roast result file
  -demeter DEMETER_SCORE_FILE
                        Demeter score file
  -toptags EDGER_TOPTAGS_FILE
                        EdgeR toptags result file
  -exprs EXPRESSION_FILE
                        Gene expression file
  -cnv COPY_NUMBER_VARIATION_FILE
                        Copy number variation file
  -meth METHYLATION_BETA_FILE
                        Methylation Beta values for individual gene CpG regions
  -out OUTPUT_FILE      Ouput file

```

### License
This project is licensed under the MIT license.
