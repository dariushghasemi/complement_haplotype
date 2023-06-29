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
    "data/dosage_of_complement_variants.txt",
    # defining columns names
    col.names  = c("AID", "chromosome", "position", "MARKER_ID", "REF", "ALT", "AF", "Dosage"),
    # defining columns classes
    colClasses = c(rep("character", 6), rep("numeric", 2))
    )

str(genotype)
#----------#

variant_annotation <- data.frame(
    "CHROM" = c(rep(10, each = 7)),
    "POS"   = c(44854402, 54528236, 54531226, 54531242, 54531685, 54533360, 54540783),
    "Gene"  = c("CXCL12", rep("MBL2", each = 6)),
    "VEP_annot" = c("intronic", "missense", "missense", "missense", "upstream gene", "upstream gene", "intergenic")
    )
#----------#

# joining phenotype and genotype data
pheno_geno <- genotype %>%
    # filtering the variant in CXCL12 gene for main analysis
    filter(position != 44854402) %>%
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
cat(" Phenotype and genotype data are defined for haplo.GLM!        ")
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
      "(?<=\\.)\\d{1}(?!\\d)",
      sprintf("%02d", as.numeric(str_extract(Haplotype, "(?<=\\.)\\d{1}(?!\\d)")))),
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
# saving full model results
saveRDS(results, "29-Jun-23_selected_variants.RDS")

# Drop unnecessary results
results_shrinked <- results %>% select(trait_name, haplotype, tidy)

# save the results in excel file
#install.packages("writexl") #format(., digits = 17)

results %>% ungroup() %>% select(tidy) %>% unnest(tidy) %>% mutate(p.value = format(p.value, digits = 18)) %>%
	write.csv("30-Jun-23_selected_variants.csv", row.names = F, quote = F)

#----------#
# show the structure of the results of the second step  
results_shrinked %>% unnest(tidy)

#----------#
cat("\n --------------------------------------------------------- \n")
cat(" Haplo.GLM results was saved!                                  ")
cat("\n --------------------------------------------------------- \n")

#----------#
# Heatmap

results_tidy <-
  results %>%
  select(- c(data, model, haplotype, glance)) %>%
  unnest(tidy) %>%
  ungroup() %>%
  mutate(associated = ifelse(p.value <= 0.05, "Yes", "No"),
         term = str_replace(term,
                            "(?<=\\.)\\d{1}(?!\\d)",
                            sprintf("%02d", as.numeric(str_extract(term, "(?<=\\.)\\d{1}(?!\\d)")))),
         term = str_replace(term, "haplo_genotype.", "Haplo_")) %>%
  filter(!str_detect(term, "(Intercept)|PC|Sex|Age|rare"))

#----------#
# Haplotype plots

results_haplo_plot <-
  results %>% 
  select(trait_name, haplotype) %>%
  unnest(haplotype) %>%
  ungroup() %>%
  pivot_longer(cols = - c(trait_name, Haplotype, hap.freq), #, p.value
               names_to = "SNP",
               values_to = "labeled_Allele") %>%
  mutate(POS = str_replace(SNP, "chr10.", "")) %>%
  inner_join(genotype %>% select(position, REF, ALT, AF) %>% distinct(position, .keep_all = T), by = c("POS" = "position")) %>%
  inner_join(variant_annotation %>% mutate(across(POS, as.character)), by = "POS") %>%
  #mutate(VEP_annot = str_replace_all(VEP_annot, "_variant", "")) %>%
  #select(- c(SNP.x, CHROM, POS, ID)) %>%
  #rename(SNP = SNP.y) %>%
  # In EPACTS wiki referred -> AC: Total Non-reference Allele Count
  # https://genome.sph.umich.edu/wiki/EPACTS
  mutate(#associated = ifelse(p.value <= 0.01, "Yes", "No"),
         Allele = case_when(labeled_Allele == 1 & AF >= 0.5 ~ ALT,
                            labeled_Allele == 1 & AF <  0.5 ~ REF,
                            labeled_Allele == 2 & AF <  0.5 ~ ALT,
                            labeled_Allele == 2 & AF >= 0.5 ~ REF)) %>%
  select(- c(REF, ALT))

#----------#

cat("\n --------------------------------------------------------- \n")
cat(" Drawing haplotypes plots . . .                                ")
cat("\n --------------------------------------------------------- \n")

results_haplo_plot

haplo_plot <- function(trait, df){
  
  my_plt <- df %>%
    ggplot(aes(vep2, term2)) +
    geom_point(aes(color = diallelic), size = 5.3, alpha = .75, show.legend = F) +
    geom_text(aes(label = Allele), color = "grey20", size = 4, vjust = .45) +
    facet_wrap(~ SNP2, scales = "free_x", nrow = 1) +
    geom_hline(yintercept = max(df$N_haplo) - .5, lty = 1, linewidth = .7, color = "grey50") +
    scale_color_manual(values = c("deepskyblue1", "green1", "magenta1", "#FF3434", "gold1", "grey50", "navyblue", "orange2")) +
    #labs(x = "", y = "") + # paste(trait, "= haplotype + Sex + Age + PC1:10")
    theme_classic() +
    theme(panel.background = element_rect(fill = "white"),
          panel.spacing = unit(0, "lines"),
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text.x    = element_text(size = 8, face = "bold", angle = 35, vjust = 1, hjust = 1.0),
          legend.position = "bottom",
          legend.key.size = unit(.2, 'cm'),
          legend.title    = element_text(size = 12, face = "bold"),
          legend.text     = element_text(size = 12),
          axis.text.x     = element_text(size = 8, face = "bold", angle = 25, vjust = 1.2, hjust = 1.1),
          axis.text.y     = element_text(size = 8, face = "bold"),
          axis.title      = element_blank())
  
  return(my_plt)
}
#----------#

haplo_plt_full <- results_tidy %>%
  group_by(trait_name) %>%
  right_join(results_haplo_plot, by = c("term" = "Haplotype", "trait_name")) %>%
  filter(term != "Haplo_rare",
  #       term == "Reference" | p.value < 0.05
  ) %>%
  group_by(trait_name, SNP) %>%
  mutate(N_haplo   = n_distinct(term),
         N_allele  = n_distinct(Allele),
         diallelic = if_else(N_allele == 1, "", Allele),
         SNP       = str_replace(SNP, "chr10.", "10:")) %>%
  # shrinking the plot to colored variants
  #filter(N_allele == 2) %>% 
  ungroup() %>%
  mutate(term2 = paste0(term, " (", round(hap.freq, 2), ", ", round(estimate, 2), ", ", round(p.value, 3), ")"),
         term2 = str_replace_all(term2, "NA, NA", "Beta, P-value"),
         SNP2 = factor(SNP, levels = unique(SNP), ordered = TRUE),
         vep2 = factor(paste0(VEP_annot, "_", Gene), levels = unique(paste0(VEP_annot, "_", Gene))),
         ) %>%
  group_by(trait_name) %>%
  nest() %>%
  mutate(#plot = data %>% map(point_range_automatic),
    plot = map2(trait_name, data, haplo_plot),
    filename = paste0("30-Jun-23_", trait_name, "_reconstructed_haplotypes_selected.png")) %>%
  ungroup() %>%
  select(plot, filename)

#----------#
# saving haplotypes plots in ong format
walk2(haplo_plt_full$plot,
      haplo_plt_full$filename,
      ~ ggsave(plot = .x, filename = .y, width = 7.5, height = 4, dpi = 300, units = "in"))


#----------#
cat("\n --------------------------------------------------------- \n")
cat(" Haplotypes plot is saved successfully!                        ")
cat("\n --------------------------------------------------------- \n")

#----------#

# From now on the rest of the nalysis will be done on personal R core machine (Thu, 14:00, 29-Jun-23).
# Darius