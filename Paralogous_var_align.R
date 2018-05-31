Packages = c("tidyverse", "dplyr", "ggplot2", "ggsignif", "huxtable")
lapply(Packages, library, character.only = TRUE)
# library(tidyverse)
# library(dplyr)
# library(ggplot2)
# library(ggsignif)
ref_data = read.csv(file="/media/nick/Data/Users/N/Documents/PhD/Paralogues/data_files/clinvar_20171029.vcf_onlyVariantslist", sep="\t", header=FALSE, col.names = c("Variant_pos", "ID"))

#Pathogenic

#python processed files
# p.paralog_data = read.csv(file="/media/nick/Data/Users/N/Documents/PhD/Paralogues/data_files/clinvar_20171029_onlyPathogenic.out_paraloc_paralogs2.noQC", sep="\t", header=FALSE, col.names = c("Variant_pos", "ID", paste("paralog", 1:132, sep = "")))
p.paralog_data = read.csv(file="/media/nick/Data/Users/N/Documents/PhD/Paralogues/data_files/clinvar_20171029_onlyPathogenic.out_paraloc_paralogs2.para_con", sep="\t", header=FALSE, col.names = c("Variant_pos", "ID", paste("paralog", 1:132, sep = "")))
# p.paralog_data = read.csv(file="/media/nick/Data/Users/N/Documents/PhD/Paralogues/data_files/clinvar_20171029_onlyPathogenic.out_paraloc_paralogs2", sep="\t", header=FALSE, col.names = c("Variant_pos", "ID", paste("paralog", 1:132, sep = "")))
# p.paralog_data = read.csv(file="/media/nick/Data/Users/N/Documents/PhD/Paralogues/data_files/clinvar_20171029_onlyPathogenic.out_paraloc_paralogs2.1", sep="\t", header=FALSE, col.names = c("Variant_pos", "ID", paste("paralog", 1:132, sep = "")))
#p.paralog_data = read.csv(file="/media/nick/Data/Users/N/Documents/PhD/Paralogues/data_files/clinvar_20171029_onlyPathogenic.out_paraloc_paralogs2.2", sep="\t", header=FALSE, col.names = c("Variant_pos", "ID", paste("paralog", 1:132, sep = "")))

#tableized files
p.tableized_data = read.csv(file="/media/nick/Data/Users/N/Documents/PhD/Paralogues/data_files/clinvar_20171029_onlyPathogenic.out_paraloc_tableized", sep = "\t", header=TRUE, stringsAsFactors=FALSE)
p.tableized_data$Variant_pos = paste(p.tableized_data$CHROM,p.tableized_data$POS, sep = " ")
p.tableized_data$REF_Amino_acids = sapply(p.tableized_data[,"Amino_acids"],strsplit, "/")
p.tableized_data$REF_Amino_acids = sapply(p.tableized_data[,"REF_Amino_acids"],unlist)
p.tableized_data$REF_Amino_acids = sapply(p.tableized_data[,"REF_Amino_acids"],function(x) x[1])
p.tableized_data$ALT_Amino_acids = sapply(p.tableized_data[,"Amino_acids"],strsplit, "/")
p.tableized_data$ALT_Amino_acids = sapply(p.tableized_data[,"ALT_Amino_acids"],unlist)
p.tableized_data$ALT_Amino_acids = sapply(p.tableized_data[,"ALT_Amino_acids"],function(x) x[2])

p.tableized_data = p.tableized_data[,c("Variant_pos","ID","REF","ALT","REF_Amino_acids","ALT_Amino_acids","Codons","Paralogue_Vars")]

# p.tableized_data = p.tableized_data[!is.na(p.tableized_data[,"Paralogue_Vars"]),]
p.tableized_data = p.tableized_data[,c("Variant_pos","ID","REF","ALT","REF_Amino_acids","ALT_Amino_acids","Codons")]

# p.tableized_data$Paralogue_Vars = lapply(p.tableized_data[,"Paralogue_Vars"],strsplit,"&")
# p.tableized_data$Paralogue_Vars = lapply(p.tableized_data[,"Paralogue_Vars"],unlist)
# p.tableized_data$Paralogue_Vars = lapply(p.tableized_data[,"Paralogue_Vars"],function(x) x[-1])
# p.tableized_data$Paralogue_Vars = lapply(p.tableized_data$Paralogue_Vars,function(x) lapply(x,strsplit,"_"))
# p.tableized_data$Paralogue_Vars = lapply(p.tableized_data[,"Paralogue_Vars"],unlist)
# p.tableized_data$Paralogue_Vars = lapply(p.tableized_data[,"Paralogue_Vars"],function(x) x[2])

p.paralog_data = left_join(p.paralog_data,p.tableized_data, by = c("Variant_pos" = "Variant_pos", "ID" = "ID"))

# p.ref_data = read.csv(file="/media/nick/Data/Users/N/Documents/PhD/Paralogues/data_files/clinvar_20171029_onlyPathogenic.vcf_onlyVariantslist", sep="\t", header=FALSE, col.names = c("Variant_pos", "ID"))
# p.ref_data = read.csv(file="/media/nick/Data/Users/N/Documents/PhD/Paralogues/data_files/clinvar_20171029_onlyPathogenic.vcf_tableized", sep = "\t", header=TRUE, stringsAsFactors=FALSE)
# p.ref_data$Variant_pos = paste(p.ref_data$CHROM,p.ref_data$POS, sep = " ")
# p.ref_data = p.ref_data[,c("Variant_pos","ID","REF","ALT")]
p.ref_data = p.tableized_data

p.gathered_paralog_data = filter(gather(p.paralog_data, paralog, paralog_pos, paste("paralog", 1:132, sep = ""), factor_key = TRUE), paralog_pos != "")
ptop.Total_paralog_annotations = left_join(p.gathered_paralog_data,p.ref_data, by = c("paralog_pos" = "Variant_pos"))
ptop.Total_paralog_annotations[is.na(ptop.Total_paralog_annotations)] = 0
ptop.num_of_paralog_anno = sum(ptop.Total_paralog_annotations$ID.y!=0)

#Ref alleles
ptop.Total_paralog_annotations = ptop.Total_paralog_annotations[ptop.Total_paralog_annotations$REF_Amino_acids.x==ptop.Total_paralog_annotations$REF_Amino_acids.y,]
ptop.num_of_paralog_anno = sum(ptop.Total_paralog_annotations$ID.y!=0)

#Alt alleles
ptop.Total_paralog_annotations = ptop.Total_paralog_annotations[ptop.Total_paralog_annotations$ALT_Amino_acids.x==ptop.Total_paralog_annotations$ALT_Amino_acids.y,]
ptop.Total_paralog_annotations = ptop.Total_paralog_annotations[!duplicated(ptop.Total_paralog_annotations$ID.x),]
ptop.Total_paralog_annotations[is.na(ptop.Total_paralog_annotations)] = 0
ptop.num_of_paralog_anno = sum(ptop.Total_paralog_annotations$ID.y!=0)


ptoa.Total_paralog_annotations = left_join(p.gathered_paralog_data,ref_data, by = c("paralog_pos" = "Variant_pos"))
ptoa.num_of_paralog_anno = sum(!is.na(ptoa.Total_paralog_annotations$ID.y))

p.num_interactions_per_annotated_var = ptop.num_of_paralog_anno/length(p.paralog_data$Variant_pos)



#Benign
b.paralog_data = read.csv(file="/media/nick/Data/Users/N/Documents/PhD/Paralogues/data_files/clinvar_20171029_onlyBenign.out_paraloc_paralogs2.noQC", sep="\t", header=FALSE, col.names = c("Variant_pos", "ID", paste("paralog", 1:132, sep = "")))
# b.paralog_data = read.csv(file="/media/nick/Data/Users/N/Documents/PhD/Paralogues/data_files/clinvar_20171029_onlyBenign.out_paraloc_paralogs2", sep="\t", header=FALSE, col.names = c("Variant_pos", "ID", paste("paralog", 1:132, sep = "")))
#b.paralog_data = read.csv(file="/media/nick/Data/Users/N/Documents/PhD/Paralogues/data_files/clinvar_20171029_onlyBenign.out_paraloc_paralogs2.1", sep="\t", header=FALSE, col.names = c("Variant_pos", "ID", paste("paralog", 1:132, sep = "")))
#b.paralog_data = read.csv(file="/media/nick/Data/Users/N/Documents/PhD/Paralogues/data_files/clinvar_20171029_onlyBenign.out_paraloc_paralogs2.2", sep="\t", header=FALSE, col.names = c("Variant_pos", "ID", paste("paralog", 1:132, sep = "")))

b.ref_data = read.csv(file="/media/nick/Data/Users/N/Documents/PhD/Paralogues/data_files/clinvar_20171029_onlyBenign.vcf_onlyVariantslist", sep="\t", header=FALSE, col.names = c("Variant_pos", "ID"))

b.gathered_paralog_data = filter(gather(b.paralog_data, paralog, paralog_pos, paste("paralog", 1:132, sep = ""), factor_key = TRUE), paralog_pos != "")
btob.Total_paralog_annotations = left_join(b.gathered_paralog_data,b.ref_data, by = c("paralog_pos" = "Variant_pos"))
btob.num_of_paralog_anno = sum(!is.na(btob.Total_paralog_annotations$ID.y))

btoa.Total_paralog_annotations = left_join(b.gathered_paralog_data,ref_data, by = c("paralog_pos" = "Variant_pos"))
btoa.num_of_paralog_anno = sum(!is.na(btoa.Total_paralog_annotations$ID.y))

b.num_interactions_per_annotated_var = btob.num_of_paralog_anno/length(b.paralog_data$Variant_pos)

#Pathogenic to Benign
ptob.Total_paralog_annotations = left_join(p.gathered_paralog_data, b.ref_data, by = c("paralog_pos" = "Variant_pos"))
ptob.num_of_paralog_anno = sum(!is.na(ptob.Total_paralog_annotations$ID.y))

#Benign to Pathogenic
btop.Total_paralog_annotations = left_join(b.gathered_paralog_data, p.ref_data, by = c("paralog_pos" = "Variant_pos"))
btop.num_of_paralog_anno = sum(!is.na(btop.Total_paralog_annotations$ID.y))

#contingency table
con_table = matrix(c(ptop.num_of_paralog_anno, ptob.num_of_paralog_anno, btop.num_of_paralog_anno, btob.num_of_paralog_anno),ncol = 2)
fisher.test(con_table)
fisher.test(matrix(c(50425,4819,34843,196), ncol = 2)) #What James did
# fisher.test(matrix(c(50425,308,34843,219), ncol = 2)) #What James did

#for variants aligning to pathogenic variants
con_table2 = matrix(
  c(length(b.paralog_data$ID),
    length(p.paralog_data$ID),
    btop.num_of_paralog_anno,
    ptop.num_of_paralog_anno
  ), ncol = 2
)
colnames(con_table2) = c("Number of variants in genes with paralogs", "Number of variants aligned to other pathogenic variants")
rownames(con_table2) = c("benign", "pathognic")

# con_table2 = as.table(con_table2)

#for variants aligning to benign variants
con_table3 = matrix(
  c(length(b.paralog_data$ID),
    length(p.paralog_data$ID),
    btob.num_of_paralog_anno,
    ptob.num_of_paralog_anno
  ), ncol = 2
)
colnames(con_table3) = c("Number of variants in genes with paralogs", "Number of variants aligned to other benign variants")
rownames(con_table3) = c("benign", "pathognic")

#stats
#for con_table2 (variants aligned to pathogenic)
con_table2_TP = ptop.num_of_paralog_anno
con_table2_FP = btop.num_of_paralog_anno
con_table2_TN = length(b.paralog_data$ID)-btop.num_of_paralog_anno
con_table2_FN = length(p.paralog_data$ID)-ptop.num_of_paralog_anno

con_table2_PPV = con_table2_TP/(con_table2_TP+con_table2_FP)
con_table2_Sensitivty = con_table2_TP/(con_table2_TP+con_table2_FN)
con_table2_Specifity = con_table2_TN/(con_table2_TN+con_table2_FP)

# con_table2_Prevalence = length(p.paralog_data$ID)/(length(p.paralog_data$ID)+length(b.paralog_data$ID))
con_table2_Prevalence = length(p.ref_data$ID)/(length(p.ref_data$ID)+length(b.ref_data$ID))
minus_con_table2_Specifity = 1 - as.numeric(con_table2_Specifity)
minus_con_table2_Prevalence = 1 - as.numeric(con_table2_Prevalence)

con_table2_PPV_prevalence = (con_table2_Sensitivty * con_table2_Prevalence) / ((con_table2_Sensitivty * con_table2_Prevalence) + (minus_con_table2_Specifity * minus_con_table2_Prevalence))

con_table2_FPR = con_table2_FP/(con_table2_FP+con_table2_TN)

#for con_table3 (variants aligned to benign)
con_table3_TP = btob.num_of_paralog_anno
con_table3_FP = ptob.num_of_paralog_anno
con_table3_TN = length(p.paralog_data$ID)-ptob.num_of_paralog_anno
con_table3_FN = length(b.paralog_data$ID)-btob.num_of_paralog_anno

con_table3_PPV = con_table3_TP/(con_table3_TP+con_table3_FP)
con_table3_Sensitivty = con_table3_TP/(con_table3_TP+con_table3_FN)
con_table3_FPR = con_table3_FP/(con_table3_FP+con_table3_TN)
##graphs
#barplot
# x_vars = c(
#   "total \n pathogenic \n variant",
#   "total \n benign \n variant",
#   "pathogenic \n variants only \n in genes \n w/ paralogs",
#   "benign \n variants only \n in genes \n w/ paralogs",
#   "pathogenic \n variants \n align to pathogenic",
#   "benign \n variants align \n to pathogenic",
#   "pathogenic \n variants \n align to benign",
#   "benign \n variants align \n to benign"
#   )
x_vars = c(
  "total variants",
  "total variants",
  "variants only \n in genes \n w/ paralogs",
  "variants only \n in genes \n w/ paralogs",
  "variants \n aligned to \n pathogenic",
  "variants \n aligned to \n pathogenic",
  "variants \n aligned to \n benign",
  "variants \n aligned to \n benign"
)
clinical_significance = c(
  "pathogenic","benign","pathogenic","benign","pathogenic","benign","pathogenic","benign"
  )
y_vars = c(
  length(p.ref_data$ID),
  length(b.ref_data$ID),
  length(p.paralog_data$ID),
  length(b.paralog_data$ID),
  ptop.num_of_paralog_anno,
  btop.num_of_paralog_anno,
  ptob.num_of_paralog_anno,
  btob.num_of_paralog_anno
           )
dt = data.frame(
  x_vars,
  y_vars,
  clinical_significance
  )
dt$x_vars <- factor(dt$x_vars, levels = unique(dt$x_vars))
leg = c("pathogenic", "benign")

g = ggplot(dt, aes(x=x_vars, y=y_vars, group = clinical_significance, fill=clinical_significance))
g + geom_bar(stat="identity", position = "dodge") + 
  # geom_signif(
  #   y_position = c(50000,47000,44000), xmin = c(5,5,5), xmax = c(8,7,6), annotations = c("*","*","*"))
  # geom_signif(
  #   y_position = c(10000), xmin = c(3), xmax = c(3.45), annotations = c("*"), tip_length = 0) +
  labs(x = "", y = "number of variants") +
  geom_text(aes(label=y_vars), position=position_dodge(width=0.9), vjust=-0.25)

#All together chart
x_vars = c(
  "total \n variants",
  "total \n variants",
  "variants only \n in genes \n w/ paralogs",
  "variants only \n in genes \n w/ paralogs",
  "variants \n aligned to \n pathogenic \n no QC",
  "variants \n aligned to \n pathogenic \n no QC",
  "variants \n aligned to \n benign \n no QC",
  "variants \n aligned to \n benign \n no QC",
  "variants \n aligned to \n pathogenic \n QC1",
  "variants \n aligned to \n pathogenic \n QC1",
  "variants \n aligned to \n benign \n QC1",
  "variants \n aligned to \n benign \n QC1",
  "variants \n aligned to \n pathogenic \n QC2",
  "variants \n aligned to \n pathogenic \n QC2",
  "variants \n aligned to \n benign \n QC2",
  "variants \n aligned to \n benign \n QC2"
)
clinical_significance = c(
  "pathogenic","benign","pathogenic","benign","pathogenic","benign","pathogenic","benign","pathogenic","benign","pathogenic","benign","pathogenic","benign","pathogenic","benign"
)
y_vars = c(
  50425,
  34843,
  13611,
  5300,
  4819,
  196,
  308,
  219,
  4437,
  44,
  205,
  92,
  2338,
  12,
  104,
  32
)
dt = data.frame(
  x_vars,
  y_vars,
  clinical_significance
)
dt$x_vars <- factor(dt$x_vars, levels = unique(dt$x_vars))
leg = c("pathogenic", "benign")

g = ggplot(dt, aes(x=x_vars, y=y_vars, group = clinical_significance, fill=clinical_significance))
g + geom_bar(stat="identity", position = "dodge") + 
  # geom_signif(
  #   y_position = c(50000,47000,44000), xmin = c(5,5,5), xmax = c(8,7,6), annotations = c("*","*","*"))
  # geom_signif(
  #   y_position = c(10000), xmin = c(3), xmax = c(3.45), annotations = c("*"), tip_length = 0) +
  labs(x = "", y = "number of variants") +
  geom_text(aes(label=y_vars), position=position_dodge(width=0.9), vjust=-0.25)





#GeneOntology
all_clinvar_genes = read.csv(
  file = "/media/nick/Data/Users/N/Documents/PhD/Paralogues/all_clinvar_genes_protein_class_Chart.txt",
  sep="\t",
  header = FALSE,
  colClasses = c("NULL",NA,NA,"NULL","NULL"),
  col.names = c("","Protein_class", "num_genes_all_clinvar", "", "")
)

# all_clinvar_genes$percent_of_tot = all_clinvar_genes$num_genes/sum(all_clinvar_genes$num_genes)

patho_paralog_snp_genes = read.csv(
  file = "/media/nick/Data/Users/N/Documents/PhD/Paralogues/patho_paralog_snps_protein_class_Chart.txt",
  sep="\t",
  header = FALSE,
  colClasses = c("NULL",NA,NA,"NULL","NULL"),
  col.names = c("","Protein_class", "num_genes_patho", "", "")
)

# patho_paralog_snp_genes$percent_of_tot = patho_paralog_snp_genes$num_genes/sum(patho_paralog_snp_genes$num_genes)
# 


benign_paralog_snp_genes = read.csv(
  file = "/media/nick/Data/Users/N/Documents/PhD/Paralogues/benign_paralog_snps_protein_class_Chart.txt",
  sep="\t",
  header = FALSE,
  colClasses = c("NULL",NA,NA,"NULL","NULL"),
  col.names = c("","Protein_class", "num_genes_benign", "", "")
)

# benign_paralog_snp_genes$percent_of_tot = benign_paralog_snp_genes$num_genes/sum(benign_paralog_snp_genes$num_genes)

all_genes = full_join(benign_paralog_snp_genes, patho_paralog_snp_genes, by="Protein_class")
all_genes = full_join(all_genes, all_clinvar_genes, by="Protein_class")
all_genes[13,"num_genes_patho"] = 0

all_clinvar_genes_bar = ggplot(all_genes, aes(x=reorder(Protein_class, num_genes_all_clinvar), y=num_genes_all_clinvar, fill=reorder(Protein_class, num_genes_all_clinvar)))
all_clinvar_genes_bar + geom_bar(stat="identity", position = "dodge") +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  geom_text(aes(label=num_genes_all_clinvar), position=position_dodge(width=0.9), vjust=0.5) +
  labs(x = "Protein_class", y = "hits") +
  guides(fill = FALSE) + coord_flip()

patho_paralog_snp_genes_bar = ggplot(all_genes, aes(x=reorder(Protein_class, num_genes_patho), y=num_genes_patho, fill=reorder(Protein_class, num_genes_all_clinvar)))
patho_paralog_snp_genes_bar + geom_bar(stat="identity", position = "dodge") +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  geom_text(aes(label=num_genes_patho), position=position_dodge(width=0.9), vjust=0.5) +
  labs(x = "Protein_class", y = "hits") +
  guides(fill = FALSE) + coord_flip()

benign_paralog_snp_genes_bar = ggplot(all_genes, aes(x=reorder(Protein_class, num_genes_benign), y=num_genes_benign, fill=reorder(Protein_class, num_genes_all_clinvar)))
benign_paralog_snp_genes_bar + geom_bar(stat="identity", position = "dodge") +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  geom_text(aes(label=num_genes_benign), position=position_dodge(width=0.9), vjust=0.5) +
  labs(x = "Protein_class", y = "hits") +
  guides(fill = FALSE) + coord_flip()

# pie = ggplot(
#   patho_paralog_snp_genes,
#   aes(x=factor(1),
#       y=percent_of_tot,
#       fill=factor(Molecular_Function)
#       )
#   )
# pie = pie + geom_bar(width = 1, stat="identity")
# pie = pie + coord_polar(theta="y")


# pie = ggplot(patho_paralog_snp_genes, aes(x = factor(0), fill = factor(patho_paralog_snp_genes$percent_of_tot))) +
#   geom_bar(width = 1)
# pie + coord_polar(theta = "y")
