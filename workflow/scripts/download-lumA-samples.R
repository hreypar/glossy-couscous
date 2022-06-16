########################################################################
# Download data from TCGA
#
# Breast cancer and healthy controls
# Expression counts matrices
########################################################################
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(magrittr)

############################ read in samples table #####################
my.lumA.samples = read.csv("TCGA-BRCA-LumA-samples.csv")

############################### make GDC query #########################
tcga.brca.lumA.expression.query <- GDCquery(project = "TCGA-BRCA",
					                                           data.category = "Transcriptome Profiling",
										                                          data.type = "Gene Expression Quantification",
										                                          workflow.type = "STAR - Counts", 
															                                         legacy = FALSE,
															                                         barcode=my.lumA.samples$barcode)

# wow this is actually the download
GDCdownload(tcga.brca.lumA.expression.query)

# A SummarizedExperiment object has three main matrices 
# that can be accessed using the SummarizedExperiment package)
tcga.brca.lumA.expression = GDCprepare(tcga.brca.lumA.expression.query, summarizedExperiment=F) #fails if T

###############################  Edits
tcga.brca.lumA.expression$gene_id %>% unique() %>% length()
#[1] 60664
tcga.brca.lumA.expression$gene_name %>% unique() %>% length()
#[1] 59428

# select only protein coding genes
tcga.brca.lumA.expression$gene_type %>% table() %>% sort()

tcga.brca.lumA.expression <- tcga.brca.lumA.expression[tcga.brca.lumA.expression$gene_type == "protein_coding", ]

# prepare features
table(grepl(";", tcga.brca.lumA.expression$gene_id))
#FALSE 
#19962

table(grepl(";", tcga.brca.lumA.expression$gene_name))
#FALSE 
#19962

features = paste(tcga.brca.lumA.expression$gene_id, tcga.brca.lumA.expression$gene_name, sep = ";")

# select only TPM
table(grepl("tpm_unstranded", colnames(tcga.brca.lumA.expression)))
#FALSE  TRUE 
#3178   635 

tpm_cols <- colnames(tcga.brca.lumA.expression)[grepl("tpm_unstranded", colnames(tcga.brca.lumA.expression))]

tcga.brca.lumA.expression <- tcga.brca.lumA.expression[, ..tpm_cols]
rm(tpm_cols)

# separate normal and cancer
healthy.bc <- paste0("tpm_unstranded_", my.lumA.samples$barcode[my.lumA.samples$sample_type == "Solid Tissue Normal"])
cancer.bc <- paste0("tpm_unstranded_", my.lumA.samples$barcode[my.lumA.samples$sample_type == "Primary Tumor"])

healthy.matrix <- cbind(features, tcga.brca.lumA.expression[, ..healthy.bc])
colnames(healthy.matrix) <- gsub(pattern = "tpm_unstranded_", "", colnames(healthy.matrix))

cancer.matrix <- cbind(features, tcga.brca.lumA.expression[, ..cancer.bc])
colnames(cancer.matrix) <- gsub(pattern = "tpm_unstranded_", "", colnames(cancer.matrix))

############################### write out tables for ARACNE
write.table(x = cancer.matrix, file = "brca/tcga-brca-lumA.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(x = features, file = "brca/tcga-brca-lumA_features.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(x = healthy.matrix, file = "healthy/tcga-brca-lumA-healthy.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(x = features, file = "healthy/tcga-brca-lumA-healthy_features.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

