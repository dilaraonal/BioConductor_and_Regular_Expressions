# BioConductor & Regular Expressions
## Bioconductor
Bioconductor is an open-source project that provides R packages specifically designed for bioinformatics and computational biology. It offers powerful tools for analyzing genomics, transcriptomics, proteomics, and metabolomics data.
**Install Bioconductor Manager**

```r
install.packages("BiocManager")
```

**Install Bioconductor Packages that are necessarily**

```r
BiocManager::install("GEOquery") 
BiocManager::install("biomaRt") 
BiocManager::install("GenomicFeatures") 
BiocManager::install("TxDb.Mmusculus.UCSC.mm9.knownGene") # The genome data of the living species used was downloaded
```
## Regular Expressions
Both tidyr and dplyr are part of the tidyverse package collection in R. They are designed for data manipulation and cleaning, making it easier to work with data frames. Regular expressions (regex) and the tidyr package are both powerful tools for data cleaning and transformation. They are often used together in text manipulation tasks when preparing datasets for analysis.

## Necessarily Library
```r
library(tidyr)
library(dplyr)
```

| Function          | Description |
|------------------|------------|
| `pivot_longer()` | 	Converts wide data to long format |
| `pivot_wider()` | Converts long data to wide format |
| `separate()` | Splits one column into multiple columns |
| `unite()` | 	Combines multiple columns into one |


## Applications 

Low carbohydrate ketogenic diet (LCKD) is used to treat epilepsy and obesity. Okuda (2019) collected transcriptomic samples from the livers and brains of regular-fed (chow diet) mice and LCKD-fed mice. The data is available in GEO (Gene Expression Omnibus) with the ID of GSE115342. There are four categories in the dataset, and three replicate samples in each category, making the total number of samples in the dataset 12. The categories are, chow-diet cortex, chow-diet liver, LCKD cortex, and LCKD liver.

**Okuda, T. (2019). A low-carbohydrate ketogenic diet promotes ganglioside synthesis via the transcriptional regulation of ganglioside metabolism-related genes. Scientific reports, 9(1), 7627.**

**a)** Use Geoquery package to automatically download the dataset. What is the dimension of the expression data in the dataset? 

**b)** Check the names of the samples from phenoData (the “title” field). From phenotype data, use regular expression to get the indices of samples for each of the four categories in the data. Your regular expressions must include at least one wildcard character. Print the number of samples for each category. 

**c)** Use t.test to calculate p-values of the genes for comparing mRNA levels of regular(chow)-diet and LCKD cases for each tissue (cortex or liver) separately. [Hint: you should put t.test into a for loop to do this. Also, the output of t-test is a list. You should record $p-value element of the list]. Plot a scatter plot to show relation between log10 of p-values of cortex and log10 of p-values of liver using ggplot2. Compare/discuss shortly the plots. 

**d)** Choose five genes with lowest p-values in the cortex data, and their associated gene symbols and ENSEMBL gene IDs using featureData. [Hint: order() function is very practical to identify indices of the lowest values in the p-value vector. After getting the order, you can get corresponding gene symbols, gene Ensembl ID, chromosome information etc]. In which chromosomes are the five genes located? Use BioMarts package to get Entrez ID’s of these genes. Then, use GenomicFeatures package to identify the genomic intervals of those genes (strand information, chromosome name and genomic range). [Hint: you can use metadata column in TxDb mouse database “TxDb.Mmusculus.UCSC.mm9.knownGene” and %in% operator to identify the matches in a single line) 

**e)** Convert the expression data into a dataframe. Take the transpose of the dataframe so that rows show samples, and columns show genes. Make sure that the output is also a dataframe by using as.data.frame. Create a new column in the transposed dataframe by using the following four category names: brain_lckd, liver_lckd, brain_chow, liver_chow. To do so, start from the phenoData “title” field. Use regular expressions with at least one wild character to convert the sample names in the phenoData “title” field to those new category names. 

**f)** It was found that LCKD diet alters expression of Gm2a gene, a gene with a role in ganglioside degradation. Create a boxplot of Gm2a gene using ggplot2 package, where each sample category should be shown with a different boxplot. Remember to enhance the visual quality of your plot by using what you learned in the class. 

**g)** For each gene in the data, calculate average mRNA expression in each sample category, and save the results into a new dataframe. The new dataframe will have 4 rows, each corresponding to a category, and about 50.000 columns, each corresponding to a gene. Replace the row names with the new category names defined in (e) and remove the category column from the dataframe. Take the transpose of this dataframe such that the rows show genes, the columns show categories. 

**h)** Make a scatter plot of mRNA levels of all genes in the brain_lckd category versus liver_lckd category. Use ggplot2 library. Remember to enhance the visual quality of your plot by using what you learned in the class. Make another plot, this time between brain_chow category and the liver_chow category. Is there any correlation between gene expression levels from two different tissues when the diet type is the same? 

**i)** Now, you want to compare expression levels of all the genes in the four categories using boxplot. For this, you will first “tidy” the data using tidyR package. Rather than having 4 columns for each gene, we want to have two columns: one listing category names, the other one listing the corresponding expression value. Use the tidy version to filter this new dataframe such that only rows having values higher than 1000 remains. Then, plot a boxplot to show the distribution of mRNA levels for each of four categories. Remember to enhance the visual quality of your plot by using jitter plots and using opacity options. Discuss the boxplot in terms of different behaviours in different categories.

