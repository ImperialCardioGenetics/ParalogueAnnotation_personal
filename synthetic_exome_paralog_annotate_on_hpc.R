Packages = c("tidyverse", "plyr", "dplyr")
new.packages = Packages[!(Packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "https://cran.ma.imperial.ac.uk/")
lapply(Packages, library, character.only = TRUE)

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

for (i in 1:22){
        synthetic_para_con.paralogous_var_align = Paralogous_var_align(
                paste0("/work/nyl112/data/synthetic_exome/compressed_para_con_and_tableized/synthetic.vep.cov.table_chrom",i,"_wIDs_proper_total.out_paraloc_paralogs2.para_con.gz"),
                paste0("/work/nyl112/data/synthetic_exome/compressed_para_con_and_tableized/synthetic.vep.cov.table_chrom",i,"_wIDs_proper_total.out_paraloc_tableized.gz"),
                "/work/nyl112/ParalogueAnnotation_personal/data/clinvar/clinvar_20190114_GRCh38_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
        save(synthetic_para_con.paralogous_var_align, file = paste0("/work/nyl112/data/synthetic_exome/compressed_para_con_and_tableized/synthetic_para_con_chrom",i,".paralogous_var_align.RData"))
        rm(synthetic_para_con.paralogous_var_align)
} 

# chrom1.paralogous_var_align = Paralogous_var_align("/work/nyl112/data/synthetic_exome/compressed_para_con_and_tableized/synthetic.vep.cov.table_chrom1_wIDs_proper_total.out_paraloc_paralogs2.para_con.gz", "/work/nyl112/data/synthetic_exome/compressed_para_con_and_tableized/synthetic.vep.cov.table_chrom1_wIDs_proper_total.out_paraloc_tableized.gz","/work/nyl112/ParalogueAnnotation_personal/data/clinvar/clinvar_20190114_GRCh38_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
# save(chrom1.paralogous_var_align, file = "/work/nyl112/data/synthetic_exome/compressed_para_con_and_tableized/synthetic_para_con_chrom1.paralogous_var_align.RData")
# rm(chrom1.paralogous_var_align)

synthetic_para_con.paralogous_var_align = Paralogous_var_align(
        "/work/nyl112/data/synthetic_exome/compressed_para_con_and_tableized/synthetic.vep.cov.table_chromX_wIDs_proper_total.out_paraloc_paralogs2.para_con.gz", 
        "/work/nyl112/data/synthetic_exome/compressed_para_con_and_tableized/synthetic.vep.cov.table_chromX_wIDs_proper_total.out_paraloc_tableized.gz",
        "/work/nyl112/ParalogueAnnotation_personal/data/clinvar/clinvar_20190114_GRCh38_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
save(synthetic_para_con.paralogous_var_align, file = "/work/nyl112/data/synthetic_exome/compressed_para_con_and_tableized/synthetic_para_con_chromX.paralogous_var_align.RData")
rm(synthetic_para_con.paralogous_var_align)

synthetic_para_con.paralogous_var_align = Paralogous_var_align(
        "/work/nyl112/data/synthetic_exome/compressed_para_con_and_tableized/synthetic.vep.cov.table_chromY_wIDs_proper_total.out_paraloc_paralogs2.para_con.gz", 
        "/work/nyl112/data/synthetic_exome/compressed_para_con_and_tableized/synthetic.vep.cov.table_chromY_wIDs_proper_total.out_paraloc_tableized.gz",
        "/work/nyl112/ParalogueAnnotation_personal/data/clinvar/clinvar_20190114_GRCh38_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized")
save(synthetic_para_con.paralogous_var_align, file = "/work/nyl112/data/synthetic_exome/compressed_para_con_and_tableized/synthetic_para_con_chromY.paralogous_var_align.RData")
rm(synthetic_para_con.paralogous_var_align)


# save(chrom1.paralogous_var_align,
#      chrom2.paralogous_var_align,
#      chrom3.paralogous_var_align,
#      chrom4.paralogous_var_align,
#      chrom5.paralogous_var_align,
#      chrom6.paralogous_var_align,
#      chrom7.paralogous_var_align,
#      chrom8.paralogous_var_align,
#      chrom9.paralogous_var_align,
#      chrom10.paralogous_var_align,
#      chrom11.paralogous_var_align,
#      chrom12.paralogous_var_align,
#      chrom13.paralogous_var_align,
#      chrom14.paralogous_var_align,
#      chrom15.paralogous_var_align,
#      chrom16.paralogous_var_align,
#      chrom17.paralogous_var_align,
#      chrom18.paralogous_var_align,
#      chrom19.paralogous_var_align,
#      chrom20.paralogous_var_align,
#      chrom21.paralogous_var_align,
#      chrom22.paralogous_var_align,
#      chromX.paralogous_var_align,
#      chromY.paralogous_var_align,
#      file = "/work/nyl112/data/synthetic_exome/compressed_para_con_and_tableized/synthetic_all_chrom_para_con.paralogous_var_align.RData"
# )