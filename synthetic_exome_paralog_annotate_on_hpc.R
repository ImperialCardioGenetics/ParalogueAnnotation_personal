Packages = c("tidyverse", "plyr", "dplyr")
new.packages = Packages[!(Packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "https://cran.ma.imperial.ac.uk/")
# lapply(Packages, library, character.only = TRUE)
library("plyr")
library("dplyr")
library("tidyverse")

Paralogous_var_align = function(paralogs2_file, 
                                paralog_tableized_file, 
                                joining_tableized_file=paralog_tableized_file, 
                                Subset=c("NA","NA"), 
                                Overlap = "NA",
                                paraz_cutoff = "NA"){ #Function for joining together variant paralogous locations
        # Subset can either be "Gene" or "ID" with a list of respective gene symbols or id values; leave as "NA" to not subset anything
        # Overlap can either be 0 to ignore variants in overlapping genes or 1 to only look at variants in overlapping genes; leave as "NA" to analyse everything
        if (endsWith(paralogs2_file, ".gz")){
                paralog_data = gzfile(paralogs2_file)
                max_no_col = (max(count.fields(paralog_data, sep = "\t"))-5) #-6 for "Variant_pos", "ID", "Gene", "Ref", "Alt", and "\n"; -5 for above python script
                paralog_data = read.csv(file=gzfile(paralogs2_file), sep="\t", header=FALSE, col.names = c("Variant_pos", "ID", "Gene", "REF", "ALT", paste("paralog", 1:max_no_col, sep = "")))
                paralog_tableized_data = read.csv(file=gzfile(paralog_tableized_file), sep = "\t", header=TRUE, stringsAsFactors=FALSE)
        } else {
                paralog_data = file(paralogs2_file)
                max_no_col = (max(count.fields(paralog_data, sep = "\t"))-5) #-6 for "Variant_pos", "ID", "Gene", "Ref", "Alt", and "\n"; -5 for above python script
                paralog_data = read.csv(file=paralogs2_file, sep="\t", header=FALSE, col.names = c("Variant_pos", "ID", "Gene", "REF", "ALT", paste("paralog", 1:max_no_col, sep = "")))
                paralog_tableized_data = read.csv(file=paralog_tableized_file, sep = "\t", header=TRUE, stringsAsFactors=FALSE)
        }
        
        if (Subset[1] == "Gene"){
                paralog_data = dplyr::filter(paralog_data, Gene %in% Subset[-1])
        } else if (Subset[1] == "ID"){
                paralog_data = dplyr::filter(paralog_data, ID %in% Subset[-1])
        } 
        
        paralog_tableized_data = paralog_tableized_data[!is.na(paralog_tableized_data$Amino_acids),]
        
        paralog_tableized_data$Variant_pos = paste(paralog_tableized_data$CHROM,paralog_tableized_data$POS, sep = " ")
        paralog_tableized_data$REF_Amino_acids = sapply(paralog_tableized_data[,"Amino_acids"],strsplit, "/")
        paralog_tableized_data$REF_Amino_acids = sapply(paralog_tableized_data[,"REF_Amino_acids"],function(x) x[1])
        paralog_tableized_data$ALT_Amino_acids = sapply(paralog_tableized_data[,"Amino_acids"],strsplit, "/")
        paralog_tableized_data$ALT_Amino_acids = sapply(paralog_tableized_data[,"ALT_Amino_acids"],function(x) x[2])
        paralog_data = left_join(paralog_data,paralog_tableized_data, by =  c("Variant_pos", "ID", "REF", "ALT", "Gene" = "SYMBOL"))
        
        paralog_data = distinct(paralog_data)
        paralog_data = paralog_data[!is.na(paralog_data$Protein_position) & !(grepl("-",paralog_data$Protein_position)),]
        
        if (Overlap == 0){#CAN UPDATE THIS TO WORK ON rsID or "." by first pasting together variant pos and gene symbol and then n_occur on that new coloumn
                n_occur = data.frame(table(paralog_data$ID))
                
                # overlapp_genes = n_occur[n_occur$Freq>1,]
                # overlapp_genes = filter(paralog_data, ID %in% overlapp_genes$Var1)
                
                n_occur = n_occur[n_occur$Freq==1,]
                paralog_data = filter(paralog_data, ID %in% n_occur$Var1)
        } else if (Overlap == 1){
                n_occur = data.frame(table(paralog_data$ID))
                n_occur = n_occur[n_occur$Freq>1,]
                paralog_data = filter(paralog_data, ID %in% n_occur$Var1)
        } else {
                n_occur = NA
        }
        
        input_variants = distinct(paralog_data[c("Variant_pos", "ID", "REF", "ALT")])
        
        if (!(paraz_cutoff == "NA")){
                paralog_data = paralog_data[!is.na(paralog_data$Para_Z_score) & paralog_data$Para_Z_score >= paraz_cutoff,]
        } else {
                paralog_data = paralog_data
        }
        
        gathered_paralog_data = filter(gather(paralog_data, paralog, paralog_pos, paste("paralog", 1:max_no_col, sep = ""), factor_key = TRUE), paralog_pos != "")
        
        joining_tableized_data = read.csv(file=joining_tableized_file, sep = "\t", header=TRUE, stringsAsFactors=FALSE)
        joining_tableized_data = joining_tableized_data[!is.na(joining_tableized_data$Amino_acids),]
        
        joining_tableized_data$Variant_pos = paste(joining_tableized_data$CHROM,joining_tableized_data$POS, sep = " ")
        joining_tableized_data$REF_Amino_acids = sapply(joining_tableized_data[,"Amino_acids"],strsplit, "/")
        joining_tableized_data$REF_Amino_acids = sapply(joining_tableized_data[,"REF_Amino_acids"],function(x) x[1])
        joining_tableized_data$ALT_Amino_acids = sapply(joining_tableized_data[,"Amino_acids"],strsplit, "/")
        joining_tableized_data$ALT_Amino_acids = sapply(joining_tableized_data[,"ALT_Amino_acids"],function(x) x[2])
        joining_tableized_data = joining_tableized_data[!is.na(joining_tableized_data$Protein_position) & 
                                                                !(grepl("-", joining_tableized_data$Protein_position)) &
                                                                !(grepl("-", joining_tableized_data$Amino_acids)) &
                                                                !(grepl("\\*", joining_tableized_data$Amino_acids)) &
                                                                grepl("/", joining_tableized_data$Amino_acids) &
                                                                str_count(joining_tableized_data$REF) == str_count(joining_tableized_data$ALT) &
                                                                str_count(joining_tableized_data$Amino_acids) == 3,]
        ref_data = joining_tableized_data
        
        Left_joined_gathered_paralog_data = left_join(gathered_paralog_data, ref_data, by = c("paralog_pos" = "Variant_pos"))
        Total_paralog_annotations = distinct(Left_joined_gathered_paralog_data)
        Total_paralog_annotations = Total_paralog_annotations[Total_paralog_annotations$ID.x != Total_paralog_annotations$ID.y,]# & #remove self annotations that arise from overlapping genes that are paralogues to itself
        # Total_paralog_annotations$Gene != Total_paralog_annotations$SYMBOL,] #&
        # !is.na(Total_paralog_annotations$Protein_position.x) & #remove variants that are not protein coding
        # !(grepl("-",p.paralogous_var_align$Total_paralog_annotations$Protein_position.x)),] #remove variants that are not snv
        
        # input_variants = distinct(Total_paralog_annotations[c("Variant_pos","ID.x","REF.x","ALT.x")])
        
        Unique_variant_annotations = Total_paralog_annotations[!is.na(Total_paralog_annotations$ID.y),]
        Unique_variant_gene_annotations = distinct(Unique_variant_annotations[c("Variant_pos", "ID.x", "Gene", "REF.x", "ALT.x")])
        Unique_variant_annotations = distinct(Unique_variant_annotations[c("Variant_pos", "ID.x", "REF.x", "ALT.x")])
        
        Unique_variant_IDs_not_annotated = setdiff(input_variants$ID,Unique_variant_annotations$ID.x)
        Unique_variants_not_annotated = distinct(Total_paralog_annotations[Total_paralog_annotations$ID.x %in% Unique_variant_IDs_not_annotated,c("Variant_pos", "ID.x", "REF.x", "ALT.x")])
        Unique_variant_genes_not_annotated = distinct(Total_paralog_annotations[Total_paralog_annotations$ID.x %in% Unique_variant_IDs_not_annotated,c("Variant_pos", "ID.x", "Gene", "REF.x", "ALT.x")])
        
        true_num_of_paralog_anno = length(Unique_variant_annotations$Variant_pos)
        num_of_paralog_anno = sum(!is.na(Total_paralog_annotations$ID.y))
        return(list("paralog_data" = paralog_data, 
                    "input_variants" = input_variants, 
                    "Unique_variant_annotations" = Unique_variant_annotations,
                    "Unique_variant_gene_annotations" = Unique_variant_gene_annotations,
                    "Unique_variant_IDs_not_annotated" = Unique_variant_IDs_not_annotated,
                    "Unique_variant_genes_not_annotated" = Unique_variant_genes_not_annotated,
                    "Unique_variants_not_annotated" = Unique_variants_not_annotated,
                    "gathered_paralog_data" = gathered_paralog_data, 
                    "Left_joined_gathered_paralog_data" = Left_joined_gathered_paralog_data,
                    "Total_paralog_annotations" = Total_paralog_annotations, 
                    "num_of_paralog_anno" = num_of_paralog_anno, 
                    "true_num_of_paralog_anno" = true_num_of_paralog_anno, 
                    "ref_data" = ref_data, 
                    "max_no_col" = max_no_col,
                    "Subset" = Subset,
                    "n_occur" = n_occur))
}

setwd("/work/nyl112/ParalogueAnnotation_personal")
# JUST COMMENT OUT THE BITS THAT DONT NEED RUNNING
# #chrom1
# for (i in 0:9){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_1/synthetic.vep.cov.table_chrom1_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_1/synthetic.vep.cov.table_chrom1_wIDs_proper_split00",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_1/synthetic.vep.cov.table_chrom1_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_1/synthetic.vep.cov.table_chrom1_wIDs_proper_split00",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_1/synthetic_chrom1_para_con_split00",i,".paralogous_var_align.RDS"))
# }
# for (i in 10:99){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_1/synthetic.vep.cov.table_chrom1_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_1/synthetic.vep.cov.table_chrom1_wIDs_proper_split0",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_1/synthetic.vep.cov.table_chrom1_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"), 
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_1/synthetic.vep.cov.table_chrom1_wIDs_proper_split0",i,".out_paraloc_tableized"), 
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_1/synthetic_chrom1_para_con_split0",i,".paralogous_var_align.RDS"))
# }
# for (i in 100:360){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_1/synthetic.vep.cov.table_chrom1_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_1/synthetic.vep.cov.table_chrom1_wIDs_proper_split",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_1/synthetic.vep.cov.table_chrom1_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"), 
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_1/synthetic.vep.cov.table_chrom1_wIDs_proper_split",i,".out_paraloc_tableized"), 
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_1/synthetic_chrom1_para_con_split",i,".paralogous_var_align.RDS"))
# }

# #chrom2
# for (i in 0:9){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_2/synthetic.vep.cov.table_chrom2_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_2/synthetic.vep.cov.table_chrom2_wIDs_proper_split00",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_2/synthetic.vep.cov.table_chrom2_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_2/synthetic.vep.cov.table_chrom2_wIDs_proper_split00",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_2/synthetic_chrom2_para_con_split00",i,".paralogous_var_align.RDS"))
# }
# for (i in 10:99){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_2/synthetic.vep.cov.table_chrom2_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_2/synthetic.vep.cov.table_chrom2_wIDs_proper_split0",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_2/synthetic.vep.cov.table_chrom2_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"), 
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_2/synthetic.vep.cov.table_chrom2_wIDs_proper_split0",i,".out_paraloc_tableized"), 
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_2/synthetic_chrom2_para_con_split0",i,".paralogous_var_align.RDS"))
# }
# for (i in 100:263){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_2/synthetic.vep.cov.table_chrom2_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_2/synthetic.vep.cov.table_chrom2_wIDs_proper_split",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_2/synthetic.vep.cov.table_chrom2_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"), 
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_2/synthetic.vep.cov.table_chrom2_wIDs_proper_split",i,".out_paraloc_tableized"), 
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_2/synthetic_chrom2_para_con_split",i,".paralogous_var_align.RDS"))
# }

# #chrom3
# for (i in 0:9){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_3/synthetic.vep.cov.table_chrom3_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_3/synthetic.vep.cov.table_chrom3_wIDs_proper_split00",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_3/synthetic.vep.cov.table_chrom3_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_3/synthetic.vep.cov.table_chrom3_wIDs_proper_split00",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_3/synthetic_chrom3_para_con_split00",i,".paralogous_var_align.RDS"))
# }
# for (i in 10:99){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_3/synthetic.vep.cov.table_chrom3_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_3/synthetic.vep.cov.table_chrom3_wIDs_proper_split0",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_3/synthetic.vep.cov.table_chrom3_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_3/synthetic.vep.cov.table_chrom3_wIDs_proper_split0",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_3/synthetic_chrom3_para_con_split0",i,".paralogous_var_align.RDS"))
# }
# for (i in 100:204){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_3/synthetic.vep.cov.table_chrom3_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_3/synthetic.vep.cov.table_chrom3_wIDs_proper_split",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_3/synthetic.vep.cov.table_chrom3_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_3/synthetic.vep.cov.table_chrom3_wIDs_proper_split",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_3/synthetic_chrom3_para_con_split",i,".paralogous_var_align.RDS"))
# }
# 
# #chrom4
# for (i in 0:9){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_4/synthetic.vep.cov.table_chrom4_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_4/synthetic.vep.cov.table_chrom4_wIDs_proper_split00",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_4/synthetic.vep.cov.table_chrom4_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_4/synthetic.vep.cov.table_chrom4_wIDs_proper_split00",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_4/synthetic_chrom4_para_con_split00",i,".paralogous_var_align.RDS"))
# }
# for (i in 10:99){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_4/synthetic.vep.cov.table_chrom4_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_4/synthetic.vep.cov.table_chrom4_wIDs_proper_split0",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_4/synthetic.vep.cov.table_chrom4_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_4/synthetic.vep.cov.table_chrom4_wIDs_proper_split0",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_4/synthetic_chrom4_para_con_split0",i,".paralogous_var_align.RDS"))
# }
# for (i in 100:139){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_4/synthetic.vep.cov.table_chrom4_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_4/synthetic.vep.cov.table_chrom4_wIDs_proper_split",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_4/synthetic.vep.cov.table_chrom4_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_4/synthetic.vep.cov.table_chrom4_wIDs_proper_split",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_4/synthetic_chrom4_para_con_split",i,".paralogous_var_align.RDS"))
# }
# 
# #chrom5
# for (i in 0:9){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_5/synthetic.vep.cov.table_chrom5_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_5/synthetic.vep.cov.table_chrom5_wIDs_proper_split00",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_5/synthetic.vep.cov.table_chrom5_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_5/synthetic.vep.cov.table_chrom5_wIDs_proper_split00",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_5/synthetic_chrom5_para_con_split00",i,".paralogous_var_align.RDS"))
# }
# for (i in 10:99){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_5/synthetic.vep.cov.table_chrom5_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_5/synthetic.vep.cov.table_chrom5_wIDs_proper_split0",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_5/synthetic.vep.cov.table_chrom5_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_5/synthetic.vep.cov.table_chrom5_wIDs_proper_split0",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_5/synthetic_chrom5_para_con_split0",i,".paralogous_var_align.RDS"))
# }
# for (i in 100:159){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_5/synthetic.vep.cov.table_chrom5_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_5/synthetic.vep.cov.table_chrom5_wIDs_proper_split",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_5/synthetic.vep.cov.table_chrom5_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_5/synthetic.vep.cov.table_chrom5_wIDs_proper_split",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_5/synthetic_chrom5_para_con_split",i,".paralogous_var_align.RDS"))
# }
# 
# #chrom6
# for (i in 0:9){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_6/synthetic.vep.cov.table_chrom6_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_6/synthetic.vep.cov.table_chrom6_wIDs_proper_split00",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_6/synthetic.vep.cov.table_chrom6_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_6/synthetic.vep.cov.table_chrom6_wIDs_proper_split00",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_6/synthetic_chrom6_para_con_split00",i,".paralogous_var_align.RDS"))
# }
# for (i in 10:99){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_6/synthetic.vep.cov.table_chrom6_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_6/synthetic.vep.cov.table_chrom6_wIDs_proper_split0",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_6/synthetic.vep.cov.table_chrom6_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_6/synthetic.vep.cov.table_chrom6_wIDs_proper_split0",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_6/synthetic_chrom6_para_con_split0",i,".paralogous_var_align.RDS"))
# }
# for (i in 100:179){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_6/synthetic.vep.cov.table_chrom6_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_6/synthetic.vep.cov.table_chrom6_wIDs_proper_split",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_6/synthetic.vep.cov.table_chrom6_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_6/synthetic.vep.cov.table_chrom6_wIDs_proper_split",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_6/synthetic_chrom6_para_con_split",i,".paralogous_var_align.RDS"))
# }
# 
# #chrom7
# for (i in 0:9){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_7/synthetic.vep.cov.table_chrom7_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_7/synthetic.vep.cov.table_chrom7_wIDs_proper_split00",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_7/synthetic.vep.cov.table_chrom7_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_7/synthetic.vep.cov.table_chrom7_wIDs_proper_split00",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_7/synthetic_chrom7_para_con_split00",i,".paralogous_var_align.RDS"))
# }
# for (i in 10:99){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_7/synthetic.vep.cov.table_chrom7_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_7/synthetic.vep.cov.table_chrom7_wIDs_proper_split0",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_7/synthetic.vep.cov.table_chrom7_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_7/synthetic.vep.cov.table_chrom7_wIDs_proper_split0",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_7/synthetic_chrom7_para_con_split0",i,".paralogous_var_align.RDS"))
# }
# for (i in 100:172){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_7/synthetic.vep.cov.table_chrom7_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_7/synthetic.vep.cov.table_chrom7_wIDs_proper_split",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_7/synthetic.vep.cov.table_chrom7_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_7/synthetic.vep.cov.table_chrom7_wIDs_proper_split",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_7/synthetic_chrom7_para_con_split",i,".paralogous_var_align.RDS"))
# }
# 
# #chrom8
# for (i in 0:9){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_8/synthetic.vep.cov.table_chrom8_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_8/synthetic.vep.cov.table_chrom8_wIDs_proper_split00",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_8/synthetic.vep.cov.table_chrom8_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_8/synthetic.vep.cov.table_chrom8_wIDs_proper_split00",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_8/synthetic_chrom8_para_con_split00",i,".paralogous_var_align.RDS"))
# }
# for (i in 10:99){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_8/synthetic.vep.cov.table_chrom8_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_8/synthetic.vep.cov.table_chrom8_wIDs_proper_split0",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_8/synthetic.vep.cov.table_chrom8_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_8/synthetic.vep.cov.table_chrom8_wIDs_proper_split0",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_8/synthetic_chrom8_para_con_split0",i,".paralogous_var_align.RDS"))
# }
# for (i in 100:123){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_8/synthetic.vep.cov.table_chrom8_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_8/synthetic.vep.cov.table_chrom8_wIDs_proper_split",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_8/synthetic.vep.cov.table_chrom8_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_8/synthetic.vep.cov.table_chrom8_wIDs_proper_split",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_8/synthetic_chrom8_para_con_split",i,".paralogous_var_align.RDS"))
# }
# 
# #chrom9
# for (i in 0:9){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_9/synthetic.vep.cov.table_chrom9_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_9/synthetic.vep.cov.table_chrom9_wIDs_proper_split00",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_9/synthetic.vep.cov.table_chrom9_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_9/synthetic.vep.cov.table_chrom9_wIDs_proper_split00",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_9/synthetic_chrom9_para_con_split00",i,".paralogous_var_align.RDS"))
# }
# for (i in 10:99){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_9/synthetic.vep.cov.table_chrom9_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_9/synthetic.vep.cov.table_chrom9_wIDs_proper_split0",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_9/synthetic.vep.cov.table_chrom9_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_9/synthetic.vep.cov.table_chrom9_wIDs_proper_split0",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_9/synthetic_chrom9_para_con_split0",i,".paralogous_var_align.RDS"))
# }
# for (i in 100:144){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_9/synthetic.vep.cov.table_chrom9_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_9/synthetic.vep.cov.table_chrom9_wIDs_proper_split",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_9/synthetic.vep.cov.table_chrom9_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_9/synthetic.vep.cov.table_chrom9_wIDs_proper_split",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_9/synthetic_chrom9_para_con_split",i,".paralogous_var_align.RDS"))
# }
# 
# #chrom10
# for (i in 0:9){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_10/synthetic.vep.cov.table_chrom10_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_10/synthetic.vep.cov.table_chrom10_wIDs_proper_split00",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_10/synthetic.vep.cov.table_chrom10_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_10/synthetic.vep.cov.table_chrom10_wIDs_proper_split00",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_10/synthetic_chrom10_para_con_split00",i,".paralogous_var_align.RDS"))
# }
# for (i in 10:99){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_10/synthetic.vep.cov.table_chrom10_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_10/synthetic.vep.cov.table_chrom10_wIDs_proper_split0",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_10/synthetic.vep.cov.table_chrom10_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_10/synthetic.vep.cov.table_chrom10_wIDs_proper_split0",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_10/synthetic_chrom10_para_con_split0",i,".paralogous_var_align.RDS"))
# }
# for (i in 100:145){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_10/synthetic.vep.cov.table_chrom10_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_10/synthetic.vep.cov.table_chrom10_wIDs_proper_split",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_10/synthetic.vep.cov.table_chrom10_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_10/synthetic.vep.cov.table_chrom10_wIDs_proper_split",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_10/synthetic_chrom10_para_con_split",i,".paralogous_var_align.RDS"))
# }
# 
# #chrom11
# for (i in 0:9){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_11/synthetic.vep.cov.table_chrom11_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_11/synthetic.vep.cov.table_chrom11_wIDs_proper_split00",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_11/synthetic.vep.cov.table_chrom11_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_11/synthetic.vep.cov.table_chrom11_wIDs_proper_split00",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_11/synthetic_chrom11_para_con_split00",i,".paralogous_var_align.RDS"))
# }
# for (i in 10:99){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_11/synthetic.vep.cov.table_chrom11_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_11/synthetic.vep.cov.table_chrom11_wIDs_proper_split0",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_11/synthetic.vep.cov.table_chrom11_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_11/synthetic.vep.cov.table_chrom11_wIDs_proper_split0",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_11/synthetic_chrom11_para_con_split0",i,".paralogous_var_align.RDS"))
# }
# for (i in 100:208){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_11/synthetic.vep.cov.table_chrom11_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_11/synthetic.vep.cov.table_chrom11_wIDs_proper_split",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_11/synthetic.vep.cov.table_chrom11_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_11/synthetic.vep.cov.table_chrom11_wIDs_proper_split",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_11/synthetic_chrom11_para_con_split",i,".paralogous_var_align.RDS"))
# }
# 
# #chrom12
# for (i in 0:9){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_12/synthetic.vep.cov.table_chrom12_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_12/synthetic.vep.cov.table_chrom12_wIDs_proper_split00",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_12/synthetic.vep.cov.table_chrom12_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_12/synthetic.vep.cov.table_chrom12_wIDs_proper_split00",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_12/synthetic_chrom12_para_con_split00",i,".paralogous_var_align.RDS"))
# }
# for (i in 10:99){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_12/synthetic.vep.cov.table_chrom12_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_12/synthetic.vep.cov.table_chrom12_wIDs_proper_split0",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_12/synthetic.vep.cov.table_chrom12_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_12/synthetic.vep.cov.table_chrom12_wIDs_proper_split0",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_12/synthetic_chrom12_para_con_split0",i,".paralogous_var_align.RDS"))
# }
# for (i in 100:191){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_12/synthetic.vep.cov.table_chrom12_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_12/synthetic.vep.cov.table_chrom12_wIDs_proper_split",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_12/synthetic.vep.cov.table_chrom12_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_12/synthetic.vep.cov.table_chrom12_wIDs_proper_split",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_12/synthetic_chrom12_para_con_split",i,".paralogous_var_align.RDS"))
# }
# 
# #chrom13
# for (i in 0:9){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_13/synthetic.vep.cov.table_chrom13_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_13/synthetic.vep.cov.table_chrom13_wIDs_proper_split00",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_13/synthetic.vep.cov.table_chrom13_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_13/synthetic.vep.cov.table_chrom13_wIDs_proper_split00",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_13/synthetic_chrom13_para_con_split00",i,".paralogous_var_align.RDS"))
# }
# for (i in 10:63){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_13/synthetic.vep.cov.table_chrom13_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_13/synthetic.vep.cov.table_chrom13_wIDs_proper_split0",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_13/synthetic.vep.cov.table_chrom13_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_13/synthetic.vep.cov.table_chrom13_wIDs_proper_split0",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_13/synthetic_chrom13_para_con_split0",i,".paralogous_var_align.RDS"))
# }
# 
# #chrom14
# for (i in 0:9){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_14/synthetic.vep.cov.table_chrom14_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_14/synthetic.vep.cov.table_chrom14_wIDs_proper_split00",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_14/synthetic.vep.cov.table_chrom14_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_14/synthetic.vep.cov.table_chrom14_wIDs_proper_split00",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_14/synthetic_chrom14_para_con_split00",i,".paralogous_var_align.RDS"))
# }
# for (i in 10:99){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_14/synthetic.vep.cov.table_chrom14_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_14/synthetic.vep.cov.table_chrom14_wIDs_proper_split0",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_14/synthetic.vep.cov.table_chrom14_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_14/synthetic.vep.cov.table_chrom14_wIDs_proper_split0",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_14/synthetic_chrom14_para_con_split0",i,".paralogous_var_align.RDS"))
# }
# for (i in 100:117){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_14/synthetic.vep.cov.table_chrom14_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_14/synthetic.vep.cov.table_chrom14_wIDs_proper_split",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_14/synthetic.vep.cov.table_chrom14_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_14/synthetic.vep.cov.table_chrom14_wIDs_proper_split",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_14/synthetic_chrom14_para_con_split",i,".paralogous_var_align.RDS"))
# }
# 
# #chrom15
# for (i in 0:9){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_15/synthetic.vep.cov.table_chrom15_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_15/synthetic.vep.cov.table_chrom15_wIDs_proper_split00",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_15/synthetic.vep.cov.table_chrom15_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_15/synthetic.vep.cov.table_chrom15_wIDs_proper_split00",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_15/synthetic_chrom15_para_con_split00",i,".paralogous_var_align.RDS"))
# }
# for (i in 10:93){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_15/synthetic.vep.cov.table_chrom15_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_15/synthetic.vep.cov.table_chrom15_wIDs_proper_split0",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_15/synthetic.vep.cov.table_chrom15_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_15/synthetic.vep.cov.table_chrom15_wIDs_proper_split0",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_15/synthetic_chrom15_para_con_split0",i,".paralogous_var_align.RDS"))
# }
# #94 1st half
# info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_15/synthetic.vep.cov.table_chrom15_wIDs_proper_split0","94_half1",".out_paraloc_paralogs2.para_con"))
# if (info$size == 0) next
# tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_15/synthetic.vep.cov.table_chrom15_wIDs_proper_split0","94_half1",".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
# if (dim(tableized_file)[1] == 0) next
# p.normal_PA = Paralogous_var_align(
#   paste0("./data/all_possible_mutations/synthetic_exome/chrom_15/synthetic.vep.cov.table_chrom15_wIDs_proper_split0","94_half1",".out_paraloc_paralogs2.para_con"),
#   paste0("./data/all_possible_mutations/synthetic_exome/chrom_15/synthetic.vep.cov.table_chrom15_wIDs_proper_split0","94_half1",".out_paraloc_tableized"),
#   "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
# save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_15/synthetic_chrom15_para_con_split0","94_half1",".paralogous_var_align.RDS"))
# #94 2nd half
# info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_15/synthetic.vep.cov.table_chrom15_wIDs_proper_split0","94_half2",".out_paraloc_paralogs2.para_con"))
# if (info$size == 0) next
# tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_15/synthetic.vep.cov.table_chrom15_wIDs_proper_split0","94_half2",".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
# if (dim(tableized_file)[1] == 0) next
# p.normal_PA = Paralogous_var_align(
#   paste0("./data/all_possible_mutations/synthetic_exome/chrom_15/synthetic.vep.cov.table_chrom15_wIDs_proper_split0","94_half2",".out_paraloc_paralogs2.para_con"),
#   paste0("./data/all_possible_mutations/synthetic_exome/chrom_15/synthetic.vep.cov.table_chrom15_wIDs_proper_split0","94_half2",".out_paraloc_tableized"),
#   "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
# save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_15/synthetic_chrom15_para_con_split0","94_half2",".paralogous_var_align.RDS"))
# for (i in 95:99){
#   info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_15/synthetic.vep.cov.table_chrom15_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
#   if (info$size == 0) next
#   tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_15/synthetic.vep.cov.table_chrom15_wIDs_proper_split0",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#   if (dim(tableized_file)[1] == 0) next
#   p.normal_PA = Paralogous_var_align(
#     paste0("./data/all_possible_mutations/synthetic_exome/chrom_15/synthetic.vep.cov.table_chrom15_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"),
#     paste0("./data/all_possible_mutations/synthetic_exome/chrom_15/synthetic.vep.cov.table_chrom15_wIDs_proper_split0",i,".out_paraloc_tableized"),
#     "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#   save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_15/synthetic_chrom15_para_con_split0",i,".paralogous_var_align.RDS"))
# }
# for (i in 100:129){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_15/synthetic.vep.cov.table_chrom15_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_15/synthetic.vep.cov.table_chrom15_wIDs_proper_split",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_15/synthetic.vep.cov.table_chrom15_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_15/synthetic.vep.cov.table_chrom15_wIDs_proper_split",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_15/synthetic_chrom15_para_con_split",i,".paralogous_var_align.RDS"))
# }

# #chrom16
# for (i in 0:9){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_16/synthetic.vep.cov.table_chrom16_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_16/synthetic.vep.cov.table_chrom16_wIDs_proper_split00",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_16/synthetic.vep.cov.table_chrom16_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_16/synthetic.vep.cov.table_chrom16_wIDs_proper_split00",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_16/synthetic_chrom16_para_con_split00",i,".paralogous_var_align.RDS"))
# }
# for (i in 10:99){
#         info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_16/synthetic.vep.cov.table_chrom16_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
#         if (info$size == 0) next
#         tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_16/synthetic.vep.cov.table_chrom16_wIDs_proper_split0",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
#         if (dim(tableized_file)[1] == 0) next
#         p.normal_PA = Paralogous_var_align(
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_16/synthetic.vep.cov.table_chrom16_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"),
#                 paste0("./data/all_possible_mutations/synthetic_exome/chrom_16/synthetic.vep.cov.table_chrom16_wIDs_proper_split0",i,".out_paraloc_tableized"),
#                 "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
#         save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_16/synthetic_chrom16_para_con_split0",i,".paralogous_var_align.RDS"))
# }
for (i in 100:157){
        info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_16/synthetic.vep.cov.table_chrom16_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"))
        if (info$size == 0) next
        info = file(paste0("./data/all_possible_mutations/synthetic_exome/chrom_16/synthetic.vep.cov.table_chrom16_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"))
        max_no_col = (max(count.fields(info, sep = "\t"))-5)
        if (max_no_col == 0) next
        tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_16/synthetic.vep.cov.table_chrom16_wIDs_proper_split",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
        if (dim(tableized_file)[1] == 0) next
        p.normal_PA = Paralogous_var_align(
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_16/synthetic.vep.cov.table_chrom16_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"),
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_16/synthetic.vep.cov.table_chrom16_wIDs_proper_split",i,".out_paraloc_tableized"),
                "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
        save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_16/synthetic_chrom16_para_con_split",i,".paralogous_var_align.RDS"))
}

#chrom17
for (i in 0:9){
        info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_17/synthetic.vep.cov.table_chrom17_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
        if (info$size == 0) next
        info = file(paste0("./data/all_possible_mutations/synthetic_exome/chrom_17/synthetic.vep.cov.table_chrom17_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
        max_no_col = (max(count.fields(info, sep = "\t"))-5)
        if (max_no_col == 0) next
        tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_17/synthetic.vep.cov.table_chrom17_wIDs_proper_split00",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
        if (dim(tableized_file)[1] == 0) next
        p.normal_PA = Paralogous_var_align(
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_17/synthetic.vep.cov.table_chrom17_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"),
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_17/synthetic.vep.cov.table_chrom17_wIDs_proper_split00",i,".out_paraloc_tableized"),
                "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
        save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_17/synthetic_chrom17_para_con_split00",i,".paralogous_var_align.RDS"))
}
for (i in 10:99){
        info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_17/synthetic.vep.cov.table_chrom17_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
        if (info$size == 0) next
        info = file(paste0("./data/all_possible_mutations/synthetic_exome/chrom_17/synthetic.vep.cov.table_chrom17_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
        max_no_col = (max(count.fields(info, sep = "\t"))-5)
        if (max_no_col == 0) next
        tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_17/synthetic.vep.cov.table_chrom17_wIDs_proper_split0",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
        if (dim(tableized_file)[1] == 0) next
        p.normal_PA = Paralogous_var_align(
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_17/synthetic.vep.cov.table_chrom17_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"),
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_17/synthetic.vep.cov.table_chrom17_wIDs_proper_split0",i,".out_paraloc_tableized"),
                "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
        save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_17/synthetic_chrom17_para_con_split0",i,".paralogous_var_align.RDS"))
}
for (i in 100:209){
        info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_17/synthetic.vep.cov.table_chrom17_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"))
        if (info$size == 0) next
        info = file(paste0("./data/all_possible_mutations/synthetic_exome/chrom_17/synthetic.vep.cov.table_chrom17_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"))
        max_no_col = (max(count.fields(info, sep = "\t"))-5)
        if (max_no_col == 0) next
        tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_17/synthetic.vep.cov.table_chrom17_wIDs_proper_split",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
        if (dim(tableized_file)[1] == 0) next
        p.normal_PA = Paralogous_var_align(
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_17/synthetic.vep.cov.table_chrom17_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"),
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_17/synthetic.vep.cov.table_chrom17_wIDs_proper_split",i,".out_paraloc_tableized"),
                "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
        save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_17/synthetic_chrom17_para_con_split",i,".paralogous_var_align.RDS"))
}

#chrom18
for (i in 0:9){
        info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_18/synthetic.vep.cov.table_chrom18_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
        if (info$size == 0) next
        info = file(paste0("./data/all_possible_mutations/synthetic_exome/chrom_18/synthetic.vep.cov.table_chrom18_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
        max_no_col = (max(count.fields(info, sep = "\t"))-5)
        if (max_no_col == 0) next
        tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_18/synthetic.vep.cov.table_chrom18_wIDs_proper_split00",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
        if (dim(tableized_file)[1] == 0) next
        p.normal_PA = Paralogous_var_align(
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_18/synthetic.vep.cov.table_chrom18_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"),
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_18/synthetic.vep.cov.table_chrom18_wIDs_proper_split00",i,".out_paraloc_tableized"),
                "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
        save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_18/synthetic_chrom18_para_con_split00",i,".paralogous_var_align.RDS"))
}
for (i in 10:56){
        info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_18/synthetic.vep.cov.table_chrom18_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
        if (info$size == 0) next
        info = file(paste0("./data/all_possible_mutations/synthetic_exome/chrom_18/synthetic.vep.cov.table_chrom18_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
        max_no_col = (max(count.fields(info, sep = "\t"))-5)
        if (max_no_col == 0) next
        tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_18/synthetic.vep.cov.table_chrom18_wIDs_proper_split0",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
        if (dim(tableized_file)[1] == 0) next
        p.normal_PA = Paralogous_var_align(
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_18/synthetic.vep.cov.table_chrom18_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"),
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_18/synthetic.vep.cov.table_chrom18_wIDs_proper_split0",i,".out_paraloc_tableized"),
                "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
        save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_18/synthetic_chrom18_para_con_split0",i,".paralogous_var_align.RDS"))
}

#chrom19
for (i in 0:9){
        info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_19/synthetic.vep.cov.table_chrom19_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
        if (info$size == 0) next
        info = file(paste0("./data/all_possible_mutations/synthetic_exome/chrom_19/synthetic.vep.cov.table_chrom19_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
        max_no_col = (max(count.fields(info, sep = "\t"))-5)
        if (max_no_col == 0) next
        tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_19/synthetic.vep.cov.table_chrom19_wIDs_proper_split00",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
        if (dim(tableized_file)[1] == 0) next
        p.normal_PA = Paralogous_var_align(
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_19/synthetic.vep.cov.table_chrom19_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"),
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_19/synthetic.vep.cov.table_chrom19_wIDs_proper_split00",i,".out_paraloc_tableized"),
                "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
        save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_19/synthetic_chrom19_para_con_split00",i,".paralogous_var_align.RDS"))
}
for (i in 10:99){
        info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_19/synthetic.vep.cov.table_chrom19_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
        if (info$size == 0) next
        info = file(paste0("./data/all_possible_mutations/synthetic_exome/chrom_19/synthetic.vep.cov.table_chrom19_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
        max_no_col = (max(count.fields(info, sep = "\t"))-5)
        if (max_no_col == 0) next
        tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_19/synthetic.vep.cov.table_chrom19_wIDs_proper_split0",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
        if (dim(tableized_file)[1] == 0) next
        p.normal_PA = Paralogous_var_align(
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_19/synthetic.vep.cov.table_chrom19_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"),
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_19/synthetic.vep.cov.table_chrom19_wIDs_proper_split0",i,".out_paraloc_tableized"),
                "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
        save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_19/synthetic_chrom19_para_con_split0",i,".paralogous_var_align.RDS"))
}
for (i in 100:226){
        info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_19/synthetic.vep.cov.table_chrom19_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"))
        if (info$size == 0) next
        info = file(paste0("./data/all_possible_mutations/synthetic_exome/chrom_19/synthetic.vep.cov.table_chrom19_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"))
        max_no_col = (max(count.fields(info, sep = "\t"))-5)
        if (max_no_col == 0) next
        tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_19/synthetic.vep.cov.table_chrom19_wIDs_proper_split",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
        if (dim(tableized_file)[1] == 0) next
        p.normal_PA = Paralogous_var_align(
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_19/synthetic.vep.cov.table_chrom19_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"),
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_19/synthetic.vep.cov.table_chrom19_wIDs_proper_split",i,".out_paraloc_tableized"),
                "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
        save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_19/synthetic_chrom19_para_con_split",i,".paralogous_var_align.RDS"))
}

#chrom20
for (i in 0:9){
        info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_20/synthetic.vep.cov.table_chrom20_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
        if (info$size == 0) next
        info = file(paste0("./data/all_possible_mutations/synthetic_exome/chrom_20/synthetic.vep.cov.table_chrom20_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
        max_no_col = (max(count.fields(info, sep = "\t"))-5)
        if (max_no_col == 0) next
        tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_20/synthetic.vep.cov.table_chrom20_wIDs_proper_split00",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
        if (dim(tableized_file)[1] == 0) next
        p.normal_PA = Paralogous_var_align(
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_20/synthetic.vep.cov.table_chrom20_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"),
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_20/synthetic.vep.cov.table_chrom20_wIDs_proper_split00",i,".out_paraloc_tableized"),
                "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
        save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_20/synthetic_chrom20_para_con_split00",i,".paralogous_var_align.RDS"))
}
for (i in 10:85){
        info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_20/synthetic.vep.cov.table_chrom20_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
        if (info$size == 0) next
        info = file(paste0("./data/all_possible_mutations/synthetic_exome/chrom_20/synthetic.vep.cov.table_chrom20_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
        max_no_col = (max(count.fields(info, sep = "\t"))-5)
        if (max_no_col == 0) next
        tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_20/synthetic.vep.cov.table_chrom20_wIDs_proper_split0",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
        if (dim(tableized_file)[1] == 0) next
        p.normal_PA = Paralogous_var_align(
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_20/synthetic.vep.cov.table_chrom20_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"),
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_20/synthetic.vep.cov.table_chrom20_wIDs_proper_split0",i,".out_paraloc_tableized"),
                "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
        save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_20/synthetic_chrom20_para_con_split0",i,".paralogous_var_align.RDS"))
}

#chrom21
for (i in 0:9){
        info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_21/synthetic.vep.cov.table_chrom21_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
        if (info$size == 0) next
        info = file(paste0("./data/all_possible_mutations/synthetic_exome/chrom_21/synthetic.vep.cov.table_chrom21_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
        max_no_col = (max(count.fields(info, sep = "\t"))-5)
        if (max_no_col == 0) next
        tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_21/synthetic.vep.cov.table_chrom21_wIDs_proper_split00",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
        if (dim(tableized_file)[1] == 0) next
        p.normal_PA = Paralogous_var_align(
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_21/synthetic.vep.cov.table_chrom21_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"),
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_21/synthetic.vep.cov.table_chrom21_wIDs_proper_split00",i,".out_paraloc_tableized"),
                "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
        save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_21/synthetic_chrom21_para_con_split00",i,".paralogous_var_align.RDS"))
}
for (i in 10:36){
        info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_21/synthetic.vep.cov.table_chrom21_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
        if (info$size == 0) next
        info = file(paste0("./data/all_possible_mutations/synthetic_exome/chrom_21/synthetic.vep.cov.table_chrom21_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
        max_no_col = (max(count.fields(info, sep = "\t"))-5)
        if (max_no_col == 0) next
        tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_21/synthetic.vep.cov.table_chrom21_wIDs_proper_split0",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
        if (dim(tableized_file)[1] == 0) next
        p.normal_PA = Paralogous_var_align(
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_21/synthetic.vep.cov.table_chrom21_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"),
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_21/synthetic.vep.cov.table_chrom21_wIDs_proper_split0",i,".out_paraloc_tableized"),
                "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
        save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_21/synthetic_chrom21_para_con_split0",i,".paralogous_var_align.RDS"))
}

#chrom22
for (i in 0:9){
        info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_22/synthetic.vep.cov.table_chrom22_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
        if (info$size == 0) next
        info = file(paste0("./data/all_possible_mutations/synthetic_exome/chrom_22/synthetic.vep.cov.table_chrom22_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
        max_no_col = (max(count.fields(info, sep = "\t"))-5)
        if (max_no_col == 0) next
        tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_22/synthetic.vep.cov.table_chrom22_wIDs_proper_split00",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
        if (dim(tableized_file)[1] == 0) next
        p.normal_PA = Paralogous_var_align(
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_22/synthetic.vep.cov.table_chrom22_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"),
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_22/synthetic.vep.cov.table_chrom22_wIDs_proper_split00",i,".out_paraloc_tableized"),
                "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
        save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_22/synthetic_chrom22_para_con_split00",i,".paralogous_var_align.RDS"))
}
for (i in 10:77){
        info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_22/synthetic.vep.cov.table_chrom22_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
        if (info$size == 0) next
        info = file(paste0("./data/all_possible_mutations/synthetic_exome/chrom_22/synthetic.vep.cov.table_chrom22_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
        max_no_col = (max(count.fields(info, sep = "\t"))-5)
        if (max_no_col == 0) next
        tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_22/synthetic.vep.cov.table_chrom22_wIDs_proper_split0",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
        if (dim(tableized_file)[1] == 0) next
        p.normal_PA = Paralogous_var_align(
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_22/synthetic.vep.cov.table_chrom22_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"),
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_22/synthetic.vep.cov.table_chrom22_wIDs_proper_split0",i,".out_paraloc_tableized"),
                "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
        save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_22/synthetic_chrom22_para_con_split0",i,".paralogous_var_align.RDS"))
}

#chromX
for (i in 0:9){
        info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_X/synthetic.vep.cov.table_chromX_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
        if (info$size == 0) next
        info = file(paste0("./data/all_possible_mutations/synthetic_exome/chrom_X/synthetic.vep.cov.table_chromX_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
        max_no_col = (max(count.fields(info, sep = "\t"))-5)
        if (max_no_col == 0) next
        tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_X/synthetic.vep.cov.table_chromX_wIDs_proper_split00",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
        if (dim(tableized_file)[1] == 0) next
        p.normal_PA = Paralogous_var_align(
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_X/synthetic.vep.cov.table_chromX_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"),
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_X/synthetic.vep.cov.table_chromX_wIDs_proper_split00",i,".out_paraloc_tableized"),
                "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
        save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_X/synthetic_chromX_para_con_split00",i,".paralogous_var_align.RDS"))
}
for (i in 10:99){
        info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_X/synthetic.vep.cov.table_chromX_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
        if (info$size == 0) next
        info = file(paste0("./data/all_possible_mutations/synthetic_exome/chrom_X/synthetic.vep.cov.table_chromX_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
        max_no_col = (max(count.fields(info, sep = "\t"))-5)
        if (max_no_col == 0) next
        tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_X/synthetic.vep.cov.table_chromX_wIDs_proper_split0",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
        if (dim(tableized_file)[1] == 0) next
        p.normal_PA = Paralogous_var_align(
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_X/synthetic.vep.cov.table_chromX_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"),
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_X/synthetic.vep.cov.table_chromX_wIDs_proper_split0",i,".out_paraloc_tableized"),
                "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
        save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_X/synthetic_chromX_para_con_split0",i,".paralogous_var_align.RDS"))
}
for (i in 100:130){
        info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_X/synthetic.vep.cov.table_chromX_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"))
        if (info$size == 0) next
        info = file(paste0("./data/all_possible_mutations/synthetic_exome/chrom_X/synthetic.vep.cov.table_chromX_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"))
        max_no_col = (max(count.fields(info, sep = "\t"))-5)
        if (max_no_col == 0) next
        tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_X/synthetic.vep.cov.table_chromX_wIDs_proper_split",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
        if (dim(tableized_file)[1] == 0) next
        p.normal_PA = Paralogous_var_align(
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_X/synthetic.vep.cov.table_chromX_wIDs_proper_split",i,".out_paraloc_paralogs2.para_con"),
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_X/synthetic.vep.cov.table_chromX_wIDs_proper_split",i,".out_paraloc_tableized"),
                "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
        save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_X/synthetic_chromX_para_con_split",i,".paralogous_var_align.RDS"))
}

#chromY
for (i in 0:9){
        info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_Y/synthetic.vep.cov.table_chromY_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
        if (info$size == 0) next
        info = file(paste0("./data/all_possible_mutations/synthetic_exome/chrom_Y/synthetic.vep.cov.table_chromY_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"))
        max_no_col = (max(count.fields(info, sep = "\t"))-5)
        if (max_no_col == 0) next
        tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_Y/synthetic.vep.cov.table_chromY_wIDs_proper_split00",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
        if (dim(tableized_file)[1] == 0) next
        p.normal_PA = Paralogous_var_align(
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_Y/synthetic.vep.cov.table_chromY_wIDs_proper_split00",i,".out_paraloc_paralogs2.para_con"),
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_Y/synthetic.vep.cov.table_chromY_wIDs_proper_split00",i,".out_paraloc_tableized"),
                "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
        save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_Y/synthetic_chromY_para_con_split00",i,".paralogous_var_align.RDS"))
}
for (i in 10:11){
        info = file.info(paste0("./data/all_possible_mutations/synthetic_exome/chrom_Y/synthetic.vep.cov.table_chromY_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
        if (info$size == 0) next
        info = file(paste0("./data/all_possible_mutations/synthetic_exome/chrom_Y/synthetic.vep.cov.table_chromY_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"))
        max_no_col = (max(count.fields(info, sep = "\t"))-5)
        if (max_no_col == 0) next
        tableized_file = read.csv(file = paste0("./data/all_possible_mutations/synthetic_exome/chrom_Y/synthetic.vep.cov.table_chromY_wIDs_proper_split0",i,".out_paraloc_tableized"), header = T, stringsAsFactors = F, sep = "\t")
        if (dim(tableized_file)[1] == 0) next
        p.normal_PA = Paralogous_var_align(
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_Y/synthetic.vep.cov.table_chromY_wIDs_proper_split0",i,".out_paraloc_paralogs2.para_con"),
                paste0("./data/all_possible_mutations/synthetic_exome/chrom_Y/synthetic.vep.cov.table_chromY_wIDs_proper_split0",i,".out_paraloc_tableized"),
                "./data/clinvar/clinvar_20181028_GRCh37_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
        save(p.normal_PA, file = paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RDS_objects/chrom_Y/synthetic_chromY_para_con_split0",i,".paralogous_var_align.RDS"))
}
