# Personalized-Disease-Diagnosis-Assistants
SCADDx and LOADDx can produce a short personalized ranked list of the most likely diseases for each patient at a requested time-point, from a Knowledge Base with thousands of diseases. They enable personalized disease diagnosis by combining patients’ temporal gene expression data with a Knowledge Base. The proposed algorithms can predict the most likely diseases for each patient along with their most affected genes, which can support healthcare professionals in their decision-making.

## How to Install and Run the Software
LOADDx and SCADDx are written in R programming language and can run on Windows, MacOS, or Linux. However, the R programming language is required to be installed before you run the software/code.

### Required Programming Language:
R version 3.6.2 or above

You can download the latest version of R from here:
* [Install R](https://www.r-project.org/)

### Required packages:
Install all the packages listed in the requirements.txt file

### Steps to run the code:
1. Download the provided R code, gene expression datasets, and knowledge base, and keep them in the same folder. 
2. Open the terminal.
3. Go to the folder where you downloaded all the codes and datasets. You can use cd command for that.
4. Run SCADDx using this command: 
```
R CMD BATCH SCADDx.R
```
5. Run LOADDx using this command: 
```
R CMD BATCH LOADDx.R
```
These commands will create a .Rout (output) file in the same folder. This .Rout file will contain all the results. 
