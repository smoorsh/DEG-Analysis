# DEG-Analysis
An R script that uses the DEseq2 library to find differential gene expression between normal and cancerous tissues and creates a volcano plot.

## Installation and Preparation Instructions
Obtain normal and cancer gene expression data from GTEX and the Genome Data Commons TCGA data. Once uploaded to your working directory, run the following commands and use MergeGEM to merge the files for preprocessing.
Be sure to replace 'CANCERTYPE' and 'cancer' with your cancer type from your file name.

```
git clone https://github.com/feltus/mergegem.git
python ./mergegem/mergegem.py CANCERTYPE_GTEX CANCERTYPE_TCGA cancer-gtex-tcga.txt
cat cancer-gtex-tcga.txt | sed -e 's/\.[0-9]*//g' -e 's/ *$//' > cancer-gtex-tcga-integer.txt
cat cancer-gtex-tcga-integer.txt | awk '!a[$1]++' > cancer-gtex-tcga-integer-unique.txt
cat cancer-gtex-tcga-integer-unique.txt| sed 's/-/_/g' > cancer-gtex-tcga-clean.txt
cat cancer-gtex-tcga-clean.txt | sed 's/\s/,/g' >  cancer-gtex-tcga-clean.csv
cat cancer-gtex-tcga.labels.txt| sed 's/-/_/g' >  cancer-gtex-tcga-dash.labels.txt
cat cancer-gtex-tcga-dash.labels.txt | sed 's/\s/,/g' | sed 's/$/,GTEX_CANCERTYPE_TCGA/' > cancer-gtex-tcga.comparison.tmp
cat cancer-gtex-tcga.comparison.tmp | sed 's/sample,label,GTEX_CANCERTYPE_TCGA/Sample,Group,Comparison/' > cancer-gtex-tcga.comparison.csv
```
You should now have a count file and an annotation file, cancer-gtex-tcgs-clean.csv and cancer-gtex-tcga.comparison.csv, respectively.

Install BiocManager if it is not installed on the Rstudio Console:
```
if (!require("BiocManager", quietly = TRUE))
```
Install the DESeq2 and any other necessary libraries if they are not already installed. Replace 'library' with the name of the library to be installed.
```
BiocManager::install("library")
```
## Clone the repository
Run the following command line to clone DEG-Analysis
```
git clone https://github.com/smoorsh/DEG-Analysis.git
```
## Run the Program
Replace each argument with the title of your input or desired output file, and replace the plot-header with the header for your volcano plot.
```
Rscript ./DEG-analysis/DEG_analysis.R countfile-name.csv annotationfile-name.csv DEGoutfile-name.csv plot-header plot-filename.png
```
