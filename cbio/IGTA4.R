# Load libraries
  library(survival)
  library(survminer)
  library(plyr)
  library(tibble)
  library(forcats)

# Pull data from cbioportal
  AKR1C3_AMP <- read.delim(file = "alterations_across_samples (2).tsv")# 0=living 1=deceased
  META_METABRIC <- read.delim(file = "brca_metabric_clinical_data (1).tsv")# 0=living 1=deceased
  AKR1C3_MICROARRAY_METABRIC <- read.delim(file = "mRNA expression z-Scores relative to all samples (log microarray).txt")

  # Merge data
    names(AKR1C3_MICROARRAY_METABRIC)[names(AKR1C3_MICROARRAY_METABRIC)=="SAMPLE_ID"] <- "Sample.ID"
    names(AKR1C3_MICROARRAY_METABRIC)[names(AKR1C3_MICROARRAY_METABRIC)=="STUDY_ID"] <- "Study.ID"
    MERGE <- merge(AKR1C3_AMP, AKR1C3_MICROARRAY_METABRIC, by = "Sample.ID")
    MERGE_2 <- merge(MERGE, META_METABRIC, by = "Patient.ID")
    MERGE_2$Overall.Survival.Status <- gsub(":.*","", MERGE_2$Overall.Survival.Status) # remove everything after :

  # Subset merged data
    colnames(MERGE_2)
    subset <- MERGE_2[c(1,5,11,17:19,22:23,27:28,32,36:37)]
    colnames(subset)
  
  # Change sructure
    str(subset)
    subset$Overall.Survival.Status <- as.numeric(as.character(subset$Overall.Survival.Status))
    str(subset)

  # Flip factor levels, WANT: no_alteration/negative = 1 and amp/positive = 2 
    subset$SGK3.x <- fct_rev(subset$SGK3.x)
    str(subset)

  # Remove missing surivival data
    subset<-na.omit(subset, cols = c("Overall.Survival..Months.", "Overall.Survival.Status"))# remove rows without OS data

# cox model for gene expression
  aka_cox <- coxph(formula = Surv(Overall.Survival..Months., Overall.Survival.Status) ~ SGK3.x, data = subset)
  beta_coef_i_want <- coef(summary(aka_cox))[1]
  hazard_ratio_i_want <- coef(summary(aka_cox))[2]
  p_value_i_want <- coef(summary(aka_cox))[5]
  summary(aka_cox)
# Histogram of expression data
  colnames(subset)
  hist(subset[,3])

# cox model
#hazard ratio - positive value indicated greater risk (opposite of linear model)
aka_cox <- coxph(formula = Surv(OS_MONTHS, OS_STATUS) ~ AKR1C3, data = sub_merge_AKR1C3_OS)
coxph(formula = Surv(OS_MONTHS, OS_STATUS) ~ AKR1C3, data = sub_merge_AKR1C3_OS)

beta_coef_i_want <- coef(summary(aka_cox))[1]
hazard_ratio_i_want <- coef(summary(aka_cox))[2]
p_value_i_want <- coef(summary(aka_cox))[5]
summary(aka_cox)
#95% CI 1.58-6.50


