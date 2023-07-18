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
    "POS"   = c(44854402, 54528236, 54531226, 54531242, 54531685, 54533360, 54540783)
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
bcftools query -f '[%SAMPLE\t%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\t%R2\t%DS\n]'
	 /shared/statgen/CHRIS5000/Imputation/HRCv1.1.new/chr10.rsq03.vcf.gz 
	 -R data/variants.list 
	 -o data/dosage_of_complement_variants.txt
```

- Haplotype association analysis using the entire 6 variants in MBL2 gene is done (Thu, 18:30, 29-Jun-23).

- Genes names and VEP annotation were also added to the plot. The analysis completed here for now after running the model on the entire 7 variants as a sensitivity analysis (Fri, 00:30, 30-Jun-23).

- The threshold for rare haplotype was lowered to `haplo.freq.min = 0.001` to have more rare haplotypes in the set (Sun, 11:30, 09-Jul-23).

- After running the main and sensitivity analyses, a legacy approach is desired now. To do so, several steps have been retaken. Three other variants added to the list of target variants. The dosage levels as well as imputation accuracy (R2) were extracted from VCF file. The analysis is rerun using the previously developed haplotype analysis (Wed, 19:00, 13-Jul-23).

- To double check the results and avoid any misspecification of alleles to the variants in the haplotypes illustrations, the characteristics of the variants which were extracted earlier from VCF file is exported in text format. Besides, the alleles in major/minor binary format were saved in a text file `13-Jul-23_lagacy_haplotypes.txt` to be compared with  `info_variants.txt` file (Tue, 15:15, 18-Jul-23).

```bash
# define column names
echo -e "AID\tCHROM\tPOS\tID\tREF\tALT\tAF\tR2\tDS" > info_variants.txt

# take the variants' information of a random participant
echo -e "AID\tCHROM\tPOS\tID\tREF\tALT\tAF\tR2\tDS" > info_variants.txt
```

Dariush