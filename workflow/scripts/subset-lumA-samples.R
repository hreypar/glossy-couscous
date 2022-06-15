########################################################################
# Download samples data from TCGA
#
# Breast cancer and healthy controls
# Luminal A subtype
#
#################### import libraries and set options ##################
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(magrittr)

########################  obtain RNA-seq ###############################
TCGAbiolinks:::getProjectSummary("TCGA-BRCA")
# Breast
# Ductal and lobular neoplasms
# Female

tcga.brca.expression.query <- GDCquery(project = "TCGA-BRCA",
				                          data.category = "Transcriptome Profiling",
							                     data.type = "Gene Expression Quantification",
							                     workflow.type = "STAR - Counts", 
									                        legacy = FALSE)

tcga.brca.expression = getResults(tcga.brca.expression.query)

################## build my samples data frame #########################
# A TCGA barcode is composed of a collection of metadata identifiers. 
# Each barcode specifically identifies a TCGA data element.
# https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/

# cases is the barcode
tcga.brca.expression$cases %>% length()
#[1] 1226
tcga.brca.expression$cases %>% unique() %>% length()
#[1] 1226

# characters 1 to 12 of the barcode means up to participant
tcga.brca.expression$cases %>% substr(1,12) %>% unique() %>% length()
#[1] 1095

# one participant can have more than one tissue sample
# characters 1 to 15 of the barcode means up to tissue sample
tcga.brca.expression$cases %>% substr(1,15) %>% unique() %>% length()
#[1] 1215

# one tissue sample can have more than one vial
# characters 1 to 16 of the barcode means up to tissue sample vial
tcga.brca.expression$cases %>% substr(1,16) %>% unique() %>% length()
#[1] 1221

# one tissue sample vial can have more than one portion
# this means a portion is a biological sample
# characters 1 to 19 of the barcode means up to portion
tcga.brca.expression$cases %>% substr(1,19) %>% unique() %>% length()
#[1] 1221

################## build my data frame
my.brca.samples <- data.frame(barcode = tcga.brca.expression$cases,
			                                    participant = substr(tcga.brca.expression$cases,1,12),
							                                  tissue_sample = substr(tcga.brca.expression$cases,1,15),
							                                  portion = substr(tcga.brca.expression$cases,1,19),
											                                sample_type = tcga.brca.expression$sample_type)

# check out which sample types we have
my.brca.samples$sample_type %>% table() %>%
	  barplot(las=1, border=F, main = "TCGA - BRCA type of sample (RNA-seq)", col ="#69b3a2", 
		            names.arg = paste0(names(table(my.brca.samples$sample_type)),
					                                    " (", table(my.brca.samples$sample_type), ")")
			              )

################## download subtypes ################## 
tcga.brca.subtypes = TCGAquery_subtype(tumor = "brca")

dim(tcga.brca.subtypes)
#[1] 1087   24

# see how many unique patients
tcga.brca.subtypes$patient %>% unique() %>% length()
#[1] 1087

# check the subtypes
table(tcga.brca.subtypes$BRCA_Subtype_PAM50)
#Basal   Her2   LumA   LumB     NA Normal 
#192     82    562    209      2     40 

# MCF-7 cell line is the Luminal A subtype 
# expressing estrogen receptor (ER)
# expressing progesterone receptor (PR)

################## keep Luminal A samples
tcga.brca.LumA <- tcga.brca.subtypes[tcga.brca.subtypes$BRCA_Subtype_PAM50 == "LumA", ]

# see how many of my patients are subtype LumA
my.brca.samples$participant %in% tcga.brca.LumA$patient %>% table()
#FALSE  TRUE 
#591   635 

my.brca.samples <- my.brca.samples[my.brca.samples$participant %in% tcga.brca.LumA$patient, ]
my.brca.samples$pam50 = "LumA"

my.brca.samples[my.brca.samples$sample_type == "Solid Tissue Normal", "pam50"] <- "normal_tissue"

my.brca.samples$pathologic_stage = lapply(my.brca.samples$participant, function(x) {
						    tcga.brca.LumA[tcga.brca.LumA$patient == x, "pathologic_stage"]
				      }) %>% unlist() 

################## figure out repeated participants or samples thing
#table(my.brca.samples$participant) %>% sort()
#my.brca.samples[my.brca.samples$participant == "TCGA-A7-A0DB", ]

################## Add clinical data ################## 
tcga.brca.clinical <- GDCquery_clinic("TCGA-BRCA","clinical")

tcga.brca.clinical %<>% 
	  select(c("bcr_patient_barcode","gender", "ajcc_pathologic_stage","race","vital_status"))

length(unique(tcga.brca.clinical$bcr_patient_barcode))
#[1] 1098

lapply(my.brca.samples$participant, function(x) {
		   which(tcga.brca.clinical$bcr_patient_barcode == x)
				      }) %>% unlist() -> indx 

my.brca.samples <- cbind(my.brca.samples, tcga.brca.clinical[indx, ])
rm(indx)

#table(my.brca.samples$participant == my.brca.samples$bcr_patient_barcode)
#TRUE 
#635 

my.brca.samples <- my.brca.samples[ , names(my.brca.samples) != "bcr_patient_barcode"]

################## Write out table ################## 
write.csv(x = my.brca.samples, file = "TCGA-BRCA-LumA-samples.csv", quote = FALSE, row.names = FALSE)

