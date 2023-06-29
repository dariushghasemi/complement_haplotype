#!/usr/bin/Rscript

#==================================#
#         Haplotypes Analysis
#         June 28, 2023
#==================================#

# R script to run the Haplo.GLM model available in Haplo.Stats 
# package version 1.8.9 for the LP activation associated variants

#----------#
# importing neccessary libtraries
library(haplo.stats)
library(tidyverse)

# finding version of the used package
packageVersion("haplo.stats")
# or
# sessionInfo()

#----------#
# reading phenotype data from txt files into dataframe
phenotype <- read.delim(
    # directory towards phenotypic data
    "/home/dnoce/ComplementSystem/data/res_date_Plate_LP_invRank.ped",
    # defining columns classes to avoid corropted class for columns
    colClasses = c(rep("character", 5), rep("numeric", 1))) %>%
    # tidying the columns names for the sake of readibility
    rename(AID = X.FAM_ID, LP_activate = ODmblpRATIOinvRankres, Sex = SEX)

str(phenotype)

# check if the participants IDs are the same in both two first columns
table(phenotype$AID == phenotype$IND_ID)


#----------#
# reading genotype data from txt files into dataframe
genotype <- read.table(
    # directory towards genetic data
    "dosage_of_complement_variants.txt",
    # defining columns names
    col.names  = c("AID", "chromosome", "position", "MARKER_ID", "REF", "ALT", "AF", "Dosage"),
    # defining columns classes
    colClasses = c(rep("character", 6), rep("numeric", 2))
    )

str(genotype)
#----------#

# joining phenotype and genotype data
pheno_geno <- genotype  %>%
    # creating informative lables for variants CHR:POS
    mutate(SNPid = str_c("chr", chromosome, ":", position)) %>%
    # reshaing the data into wide format to bring variants into the columns
    pivot_wider(
        id_cols = AID,
        names_from = SNPid, 
        values_from = Dosage
        ) %>%
    # merging with phenotype file
    inner_join(
        .,
        phenotype %>% select(AID, LP_activate),
        by = "AID"
        )

# taking a look at the merged data 
pheno_geno %>% head()

cat("\n --------------------------------------------------------- \n")
cat(" Genotype and phenotype data were successfully merged!         ")
cat("\n --------------------------------------------------------- \n")

#----------#
# Sample size and number of variables 
dim(pheno_geno)

#----------#
# selecting the variants columns from the merged data
loci <- pheno_geno %>% select(starts_with("chr10:"))

#----------#
cat("\n --------------------------------------------------------- \n")
cat(" Variants were picked up for haplotype regression!             ")
cat("\n --------------------------------------------------------- \n")

#----------#
# changing the format in a way to have two columns for each varint
# before selecting, as the dosage are imputed variants,
# the dosage level should be converted to integer.
MBL2 <- geno1to2(round(loci, 0), locus.label = colnames(loci))

#----------#
cat("\n --------------------------------------------------------- \n")
cat(" Variant were successfully splitted into minor/major alleles columns!")
cat("\n --------------------------------------------------------- \n")

#----------#
# Defining genotype and phenotype data for haplo.GLM

# Setup genotype data
haplo_genotype <- setupGeno(
    MBL2,
    miss.val = c(0, NA),
    locus.label = colnames(loci)
    )

#attributes(haplo_genotype)

#----------#
# GLM data (merged phenotype and genotype data with minor/major allele indicators)
haplo_dataset <- data.frame(
    # variants dosages
    haplo_genotype,
    # phenotype of interest
    pheno_geno %>% select(- AID, - starts_with("chr"))
    )

# taking a look at the merged data for haplo.GLM funtion
head(haplo_dataset)

#----------#
cat("\n --------------------------------------------------------- \n")
cat(" Phenotype and genotype data are defined for haplo.GLM!       ")
cat("\n --------------------------------------------------------- \n")



#----------#
# Fitting regression model
# glm fit with haplotypes, additive gender covariate on gaussian response
#----------#

# Defining Haplo.GLM model for iteration via map function
hap_model <- function(df){
  
  # making PCs vector to adjust
  PCs <- paste0("PC", 1:10, collapse = " + ")
  
  # defining model formula
  my_formula <- paste("trait ~ haplo_genotype")

  # Set a common random number seed for both models
  common_iseed <- 777
  
  # Set haplo.em.control parameters
  em_ctrl <- haplo.em.control(
    n.try = 2,
    iseed = common_iseed,
    insert.batch.size = 2,
    max.haps.limit = 4e6,
    min.posterior = 1e-5)
  
  # fiting the model
  model_fit <- haplo.glm(
    my_formula,
    family = gaussian,
    data   = df,
    na.action = "na.geno.keep",
    locus.label = colnames(loci),
    x = TRUE,
    control = haplo.glm.control(haplo.freq.min = .01,
                                em.c = em_ctrl))
  
  return(model_fit)
}

#----------#
# Retrieving haplotypes from fitted models and their frequencies
haplo_extract <- function(model) {
  
  haplo_set <- 
    summary(model)$haplotypes %>%
    rownames_to_column(var = "Haplotype") %>%
    mutate(Haplotype = str_replace(
      Haplotype,
      "(?<=\\.)\\d{1,2}(?!\\d)",
      sprintf("%03d", as.numeric(str_extract(Haplotype, "(?<=\\.)\\d{1,2}(?!\\d)")))),
      Haplotype = str_replace_all(Haplotype,
                                  c("haplo_genotype." = "Haplo_",
                                    "haplo.base"      = "Reference")))
  
  return(haplo_set)
}
#----------#

#-----------------------------------------------------#
#----------------   Haplotype-Trait   ----------------
#-----------------------------------------------------#

#----------#
cat("\n --------------------------------------------------------- \n")
cat(" Running Haplo.GLM model . . .                                 ")
cat("\n --------------------------------------------------------- \n")

#----------#
# Running the model
results <- haplo_dataset %>%
  pivot_longer(cols      = - c(haplo_genotype),
               names_to  = "trait_name",
               values_to = "trait") %>%
  group_by(trait_name) %>%
  nest() %>%
  mutate(model     = data  %>% map(hap_model),
         haplotype = model %>% map(haplo_extract),
         glance    = model %>% map(broom::glance),
         tidy      = model %>% map(broom::tidy))


#----------#
cat("\n --------------------------------------------------------- \n")
cat(" Haplo.GLM model was fitted and show the results!              ")
cat("\n --------------------------------------------------------- \n")

#----------#
# show the structure of the results of the second step  
results %>% unnest(tidy)
#----------#

# saving full model results
saveRDS(results, "29-Jun-23_full_variants.RDS")


# From now on the rest of the nalysis will be done on personal R core machine (Thu, 14:00, 29-Jun-23).
# Dariush