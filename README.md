# complement_haplotype
Haplotype analysis for complement system activation GWAS project

- Have not started yet to work on it (12-Jun-23).
- Still in Florence for summer school in Bayesian Causal Inference (14-Jun-23).
- Still need to take time and focus on the analysis plan.
- Trying to open up a spot to run the analysis in this weekend (23-Jun-23).

### Preparing genotype file and script for haplotype analysis; started on Wed, 19:40, 28-jun-23:

- Extracting dosage information of individulas listed in `sample.list` - a tab-seperated file of two columns: chromosome and position - from GRCh37 vcf file. 

```bash
# Creating variants positions file

data.frame(
    "CHROM" = c(rep(10,each = 6)),
    "POS"   = c(44854402, 54531226, 54531242, 54531685, 54533360, 54540783)
    ) %>%
    write.table(
        "variants.list", 
        col.names = F, 
        quote = F, 
        row.names = F, 
        sep = "\t")

# Checking that the number of columns in each rows are consistent
cat variants.list | awk -F\\t '{print NF}' | sort | uniq -c

# extracting genotype data of variants found in variants.list file
bcftools query -f '[%SAMPLE\t%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\t%DS\n]' /shared/statgen/CHRIS5000/Imputation/HRCv1.1.new/chr10.rsq03.vcf.gz -R variants.list -o dosage_of_complement_variants.txt
```
