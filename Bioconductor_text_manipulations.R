###A###
# The BiocManager package was installed and the necessary Bioconductor packages were downloaded through this package.
install.packages("BiocManager")
BiocManager::install("GEOquery") 
library(GEOquery) # Geoquery package was used to download the data

gse_data <- getGEO("GSE115342", AnnotGPL=TRUE) # The annotation file that came with the relevant dataset was also downloaded.
gse_data  <- gse_data [[1]] # The first element of the list object was selected to directly contain the GEO dataset object
dimension_data <- dim(gse_data) 
###B###

titles <- pData(gse_data)$title # Titles have been reached for Pheno data users
# Using wild characters, indexes of titles belonging to different categories were found with the grep function
Cortex_chow<- grep(".*_cortex_chow_.*_[1-3]$", titles) 
Cortex_LCKD <- grep(".*_cortex_LCKD_.*_[1-3]$", titles)
Liver_chow <-grep(".*_liver_chow_.*_[1-3]$", titles)
Liver_LCKD <- grep(".*_liver_LCKD_.*_[1-3]$", titles)
#The number of samples in each category was specified with the obtained inex information
print(paste("Category 1: Chow-diet cortex, Number of samples:",length(Cortex_chow)))
print(paste("Category 2: LCKD cortex, Number of samples:",length(Cortex_LCKD)))
print(paste("Category 3: Chow-diet liver, Number of samples:",length(Liver_chow)))
print(paste("Category 4: LCKD liver,  Number of samples:",length(Liver_LCKD)))

###C###
data_exprs<-exprs(gse_data) # A data set containing mRNA levels information in the data was drawn (experiment data)

# t.test was applied to compare LCKD and chow diet applied in cortex tissue and p values were reached.
cortex_p_values <- NULL # An empty variable has been created to transfer p values to
for (i in 1:nrow(data_exprs)) {
  cortex_p_values[i] <- t.test(data_exprs[i, Cortex_chow], data_exprs[i,Cortex_LCKD])$p.value
}
# t.tes was applied to compare LCKD and chow diet applied to liver tissue and p values were reached.
liver_p_values <- NULL
for (i in 1:nrow(data_exprs)) {
  liver_p_values[i] <- t.test(data_exprs[i, Liver_chow], data_exprs[i, Liver_LCKD])$p.value
}

library(ggplot2) 
# The reason for taking the -log10 of p-values is to convert small p-values into more easily interpretable values.
p_values_data <- data.frame(
  Cortex = -log10(cortex_p_values),
  Liver = -log10(liver_p_values)
)

# The relationship between cortex and liver tissues was shown on the plot
ggplot(p_values_data)+
  aes(x = Cortex, y = Liver) +
  geom_jitter(pch = ".", alpha = 0.3, color = "red") + # The reason for using geom_jitter is to observe overplotting in the geom_point function
  labs(title = " Relation Between Cortex and Liver P Values",
       x = "P values of Cortex", y = "P values of Liver") + 
  theme(plot.title = element_text(hjust = 0.5)) 


###D###

#BioMart was installed to access Entrez ID information, and GenomicFeature package was installed to learn genomic features.
BiocManager::install("biomaRt") 
BiocManager::install("GenomicFeatures") 
BiocManager::install("TxDb.Mmusculus.UCSC.mm9.knownGene") # The genome data of the living species used was downloaded

library(biomaRt)
library(GenomicFeatures)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)

feature_data <- fData(gse_data) #  gen bilgilerinin yer ald?????? feature veri k??msei kullan??ld??
sorted_genes <- order(cortex_p_values)
lowest_pval_genes <-feature_data[sorted_genes[1:5], ] #5 gene information with the smallest p value

gene_symbols <- lowest_pval_genes$GENE_SYMBOL
ensembl_ids <- lowest_pval_genes$ENSEMBL_ID
chromosome_inf <- lowest_pval_genes$CHROMOSOMAL_LOCATION
# Information about the 5 genes with the smallest p values obtained from the feature data was obtained.
for (i in 1:5) {
    cat("\nGene Symbol:", gene_symbols[i],
      "\nEnsembl ID",ensembl_ids[i],
       "\nChromosome Information:", chromosome_inf[i])}


ensembl <- useMart("ENSEMBL_MART_ENSEMBL") # # The required data set to start using the database was obtained with the useMart function
ensemble_data <- useDataset("mmusculus_gene_ensembl", ensembl) # dataset from selected database was used
entrez_ids <- getBM(attributes = c("ensembl_transcript_id", "entrezgene_id"), # The desired output is the Entrez gene information in the selected data.
                    filters = "ensembl_transcript_id", # Ensembl id information accessed with feature data added as input
                    values = ensembl_ids, #Obtained ensemble id values
                    mart = ensemble_data)

txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
indices_gene_ids <- which(as.numeric(genes(txdb)$gene_id) %in% entrez_ids$entrezgene_id) #indexes of the obtained entrez ids were found
# The desired gene information was obtained by using index information in the txdb dataset
genes_txdb <- genes(txdb)
for (i in indices_gene_ids ) {
  print(genes_txdb[i])
}

###E###

df_exprs_data <- as.data.frame(data_exprs) # expression dataset converted to data frame
t_df_exprs_data <- as.data.frame(t(df_exprs_data)) # The transposon of the data set was retrieved

pattern <- ".*_(cortex|liver)_(chow|LCKD)_.*_[1-3]$" # A pattern covering all titles was created

categories <- gsub(pattern, "\\1_\\2",titles) # Required words of the titles in the titles variable have been replaced
#titles have been adjusted to the desired final state
new_titles <- gsub("cortex_chow", "brain_chow", categories)
Sample_Titles <- gsub("cortex_LCKD", "brain_LCKD", new_titles)

t_df_exprs_data <- cbind(Sample_Titles , t_df_exprs_data) # added titles information to data frame

###F###

gm2a_gene <- feature_data[feature_data$GENE_SYMBOL == "Gm2a", ] # Gm2a gene information was extracted from the feature dataset
gm2a_exprs_levels <- t_df_exprs_data[,c(1,gm2a_gene$ID)] # The mRNA levels in the expression data set were reached with the index information.
colnames(gm2a_exprs_levels)[2] <- "Gm2a_mRNA_levels" #corrected the column name of the created data frame

ggplot(gm2a_exprs_levels) + 
  aes(x = Sample_Titles , y = Gm2a_mRNA_levels ) +
  geom_boxplot(alpha = 0)+
  geom_jitter(alpha = 0.3, #Created a point cloud
              color = "blue")  +
  labs(title = "Gm2a gene mRNAs levels",
     x = "Sample categories", y = "mRNA levels") + 
  theme(plot.title = element_text(hjust = 0.5)) 

###G###
library(dplyr)

average_mRNAs <- aggregate(. ~ Sample_Titles, data = t_df_exprs_data, FUN = mean) # They were grouped according to sample titles and averages were taken.
rownames(average_mRNAs) <- average_mRNAs$Sample_Titles # row names replaced with title names
average_mRNAs <- subset(average_mRNAs, select = -Sample_Titles) # column containing title names has been removed
average_mRNAs_t <- t(average_mRNAs) # The transposon was taken to reach the desired format

###H###
average_mRNAs_t_df <- as.data.frame(average_mRNAs_t) # The transposon received data was converted into a data frame

# A plot was drawn to observe the relationship of different tissues in the same diet type

ggplot(average_mRNAs_t_df) +
  aes(x = brain_LCKD, y = liver_LCKD) + 
  geom_jitter(pch = ".", alpha = 0.3, color = "red") +
  labs(title = "Brain versus Liver in LCKD",
       x = "Brain LCKD", y = "Liver LCKD") + 
  theme(plot.title = element_text(hjust = 0.5))


ggplot(average_mRNAs_t_df) +
  aes(x = brain_chow, y =liver_chow) + 
  geom_jitter(pch = ".", alpha = 0.3, color = "red") +
  labs(title = "Brain versus Liver in Chow-diet",
       x = "Brain Chow-diet", y = "Liver Chow-diet") + 
  theme(plot.title = element_text(hjust = 0.5)) 

###I###

library(tidyr)
library(dplyr)
# pivot_longer function was used to convert the data to the new desired format
new_format <-average_mRNAs_t_df %>% pivot_longer(cols = c("brain_LCKD","liver_LCKD","brain_chow","liver_chow"),
                                                names_to = "Sample_Category",
                                                values_to = "mRNA_Levels") %>% filter(mRNA_Levels > 1000) # Those with expression levels higher than 1000 were selected

# Plot was created to show how the mRNAs averages were distributed for each category
ggplot(new_format) + 
  aes(x = Sample_Category , y = mRNA_Levels ) +
  geom_boxplot(alpha = 0)+
  geom_jitter(alpha = 0.3, 
              color = "blue")  +
  labs(title = "Distribution of mRNA Levels",
       x = "Sample categories", y = "mRNA levels") + 
  theme(plot.title = element_text(hjust = 0.5)) 
 
