library(here)
library(tidyverse)


#Load stats
stats <- read_tsv(here("analysis/stacks/populations.sumstats.tsv"), skip=1) %>% 
  rename(Locus = "# Locus ID")

#Load chromosome lengths
fai <- read_tsv(here("analysis/stacks/populations.loci.fasta.fai"),
                     col_names = c("Locus", "len", "bytes", "len_per_line", "bytes_per_line"))

snps <- read_tsv(here("analysis/stacks/populations.snps.vcf"), comment="##") %>% 
  rename(CHROM = "#CHROM") %>%
  mutate(CHROM = paste0("CLocus_", CHROM))

duplicates <- snps %>% group_by(CHROM, POS) %>%
  filter(n() > 1)

snps <- snps %>% 
  group_by(CHROM) %>%
  filter(n() == 1) %>% # Remove loci containing more than 1 SNP
  ungroup() %>% 
  left_join(fai, by= c("CHROM" = "Locus")) %>% #Add length information
  filter(POS > 30 & len > (POS+30)) #Remove SNPs to close to loci edges

# Add MAF and filter out SNPs not present in all three combinations
snps <- snps %>% 
  rowwise() %>%
  mutate(homozygot_ref = sum(str_count(c_across(contains("trimmed")), "0/0")),
         homozygot_alt = sum(str_count(c_across(contains("trimmed")), "1/1")),
         heterozygote = sum(str_count(c_across(contains("trimmed")), "0/1")),
         Not_present = sum(str_count(c_across(contains("trimmed")), "\\./\\."))
  ) %>% 
  filter(homozygot_ref > 0 &
           homozygot_alt > 0 &
           heterozygote > 0)

snps <- snps %>% 
  mutate(MAF= ((homozygot_ref *2) + heterozygote) / ((homozygot_ref + heterozygote + homozygot_alt)*2)) %>% 
  mutate(MAF=round(ifelse((MAF > 0.5), 1-MAF, MAF), digits=2),
         Prop_samp = round((homozygot_ref + heterozygote + homozygot_alt)/(homozygot_ref + heterozygote + homozygot_alt + Not_present),digits=2),
         AF = as.numeric(str_replace_all(INFO, ".*AF=", "")))

#Now plot number of snps with values of MAF and Prop_samp
ggplot(snps, aes(x = AF, y = Prop_samp)) +
  geom_tile(aes(fill = ..count..), stat = "bin2d") +
  scale_fill_viridis_c() +
  labs(title = "Heatmap of Prop_samples vs MAF",
       x = "MAF",
       y = "Prop_samples",
       fill = "Count")

# I now take every SNP with sufficient MAF

interest <- snps %>% filter(Prop_samp == 1) %>% 
  arrange(desc(AF)) %>% 
  filter(AF > 0.3) %>% 
  mutate(coordinates=paste0(CHROM,":", POS-30, "-", POS+30)) %>% 
  select(-MAF)

write_tsv(interest, here("analysis/interest_SNPs.tsv"))

if (nrow(duplicates > 0)) {
  print("WARNING: There are SNPs in the same locus at the same location. Check them for MAF calculation")
}





