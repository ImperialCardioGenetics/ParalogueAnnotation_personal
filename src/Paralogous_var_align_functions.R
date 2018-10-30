Packages = c("tidyverse", "dplyr", "ggplot2", "ggsignif", "biomaRt", "knitr", "png", "grid", "tinytex", "pander", "kableExtra", "clusterProfiler", "org.Hs.eg.db")
lapply(Packages, library, character.only = TRUE)

Paralogous_var_align = function(paralogs2_file, paralog_tableized_file, joining_tableized_file=paralog_tableized_file){ #Function for joining together variant paralogous locations
  # system(paste("python /media/nick/Data/Users/N/Documents/PhD/Paralogues/ParalogueAnnotation_personal/src/paralogs_file_remove_last_column.py ", paralogs2_file, sep = ""))
  paralog_data = file(paralogs2_file)
  max_no_col = (max(count.fields(paralog_data, sep = "\t"))-5) #-6 for "Variant_pos", "ID", "Gene", "Ref", "Alt", and "\n"; -5 for above python script
  paralog_data = read.csv(file=paralogs2_file, sep="\t", header=FALSE, col.names = c("Variant_pos", "ID", "Gene", "REF", "ALT", paste("paralog", 1:max_no_col, sep = "")))
  paralog_tableized_data = read.csv(file=paralog_tableized_file, sep = "\t", header=TRUE, stringsAsFactors=FALSE)
  paralog_tableized_data = paralog_tableized_data[!is.na(paralog_tableized_data$Amino_acids),]
  
  paralog_tableized_data$Variant_pos = paste(paralog_tableized_data$CHROM,paralog_tableized_data$POS, sep = " ")
  paralog_tableized_data$REF_Amino_acids = sapply(paralog_tableized_data[,"Amino_acids"],strsplit, "/")
  # paralog_tableized_data$REF_Amino_acids = sapply(paralog_tableized_data[,"REF_Amino_acids"],unlist)
  paralog_tableized_data$REF_Amino_acids = sapply(paralog_tableized_data[,"REF_Amino_acids"],function(x) x[1])
  paralog_tableized_data$ALT_Amino_acids = sapply(paralog_tableized_data[,"Amino_acids"],strsplit, "/")
  # paralog_tableized_data$ALT_Amino_acids = sapply(paralog_tableized_data[,"ALT_Amino_acids"],unlist)
  paralog_tableized_data$ALT_Amino_acids = sapply(paralog_tableized_data[,"ALT_Amino_acids"],function(x) x[2])
  paralog_data = left_join(paralog_data,paralog_tableized_data, by =  c("Variant_pos", "ID", "REF", "ALT", "Gene" = "SYMBOL"))
  gathered_paralog_data = filter(gather(paralog_data, paralog, paralog_pos, paste("paralog", 1:max_no_col, sep = ""), factor_key = TRUE), paralog_pos != "")
  
  joining_tableized_data = read.csv(file=joining_tableized_file, sep = "\t", header=TRUE, stringsAsFactors=FALSE)
  joining_tableized_data = joining_tableized_data[!is.na(joining_tableized_data$Amino_acids),]
  
  joining_tableized_data$Variant_pos = paste(joining_tableized_data$CHROM,joining_tableized_data$POS, sep = " ")
  joining_tableized_data$REF_Amino_acids = sapply(joining_tableized_data[,"Amino_acids"],strsplit, "/")
  # joining_tableized_data$REF_Amino_acids = sapply(joining_tableized_data[,"REF_Amino_acids"],unlist)
  joining_tableized_data$REF_Amino_acids = sapply(joining_tableized_data[,"REF_Amino_acids"],function(x) x[1])
  joining_tableized_data$ALT_Amino_acids = sapply(joining_tableized_data[,"Amino_acids"],strsplit, "/")
  # joining_tableized_data$ALT_Amino_acids = sapply(joining_tableized_data[,"ALT_Amino_acids"],unlist)
  joining_tableized_data$ALT_Amino_acids = sapply(joining_tableized_data[,"ALT_Amino_acids"],function(x) x[2])
  ref_data = joining_tableized_data
  
  Total_paralog_annotations = left_join(gathered_paralog_data, ref_data, by = c("paralog_pos" = "Variant_pos"))
  num_of_paralog_anno = sum(!is.na(Total_paralog_annotations$ID.y))
  return(list("paralog_data" = paralog_data, "gathered_paralog_data" = gathered_paralog_data, "Total_paralog_annotations" = Total_paralog_annotations, "num_of_paralog_anno" = num_of_paralog_anno, "ref_data" = ref_data, "max_no_col" = max_no_col))
}

Paralogous_var_align_compressed = function(paralogs2_file, paralog_tableized_file, joining_tableized_file=paralog_tableized_file){ #Function for joining together variant paralogous locations - compressed files
  paralog_data = gzfile(paralogs2_file)
  max_no_col = (max(count.fields(paralog_data, sep = "\t"))-5) #-6 for "Variant_pos", "ID", "Gene", "Ref", "Alt", and "\n"; -5 for above python script
  paralog_data = read.csv(file=gzfile(paralogs2_file), sep="\t", header=FALSE, col.names = c("Variant_pos", "ID", "Gene", "REF", "ALT", paste("paralog", 1:max_no_col, sep = "")))
  paralog_tableized_data = read.csv(file=gzfile(paralog_tableized_file), sep = "\t", header=TRUE, stringsAsFactors=FALSE)
  paralog_tableized_data = paralog_tableized_data[!is.na(paralog_tableized_data$Amino_acids),]
  
  paralog_tableized_data$Variant_pos = paste(paralog_tableized_data$CHROM,paralog_tableized_data$POS, sep = " ")
  paralog_tableized_data$REF_Amino_acids = sapply(paralog_tableized_data[,"Amino_acids"],strsplit, "/")
  # paralog_tableized_data$REF_Amino_acids = sapply(paralog_tableized_data[,"REF_Amino_acids"],unlist)
  paralog_tableized_data$REF_Amino_acids = sapply(paralog_tableized_data[,"REF_Amino_acids"],function(x) x[1])
  paralog_tableized_data$ALT_Amino_acids = sapply(paralog_tableized_data[,"Amino_acids"],strsplit, "/")
  # paralog_tableized_data$ALT_Amino_acids = sapply(paralog_tableized_data[,"ALT_Amino_acids"],unlist)
  paralog_tableized_data$ALT_Amino_acids = sapply(paralog_tableized_data[,"ALT_Amino_acids"],function(x) x[2])
  paralog_data = left_join(paralog_data,paralog_tableized_data, by =  c("Variant_pos", "ID", "REF", "ALT", "Gene" = "SYMBOL"))
  gathered_paralog_data = filter(gather(paralog_data, paralog, paralog_pos, paste("paralog", 1:max_no_col, sep = ""), factor_key = TRUE), paralog_pos != "")
  
  joining_tableized_data = read.csv(file=joining_tableized_file, sep = "\t", header=TRUE, stringsAsFactors=FALSE)
  joining_tableized_data = joining_tableized_data[!is.na(joining_tableized_data$Amino_acids),]
  
  joining_tableized_data$Variant_pos = paste(joining_tableized_data$CHROM,joining_tableized_data$POS, sep = " ")
  joining_tableized_data$REF_Amino_acids = sapply(joining_tableized_data[,"Amino_acids"],strsplit, "/")
  # joining_tableized_data$REF_Amino_acids = sapply(joining_tableized_data[,"REF_Amino_acids"],unlist)
  joining_tableized_data$REF_Amino_acids = sapply(joining_tableized_data[,"REF_Amino_acids"],function(x) x[1])
  joining_tableized_data$ALT_Amino_acids = sapply(joining_tableized_data[,"Amino_acids"],strsplit, "/")
  # joining_tableized_data$ALT_Amino_acids = sapply(joining_tableized_data[,"ALT_Amino_acids"],unlist)
  joining_tableized_data$ALT_Amino_acids = sapply(joining_tableized_data[,"ALT_Amino_acids"],function(x) x[2])
  ref_data = joining_tableized_data
  
  Total_paralog_annotations = left_join(gathered_paralog_data, ref_data, by = c("paralog_pos" = "Variant_pos"))
  num_of_paralog_anno = sum(!is.na(Total_paralog_annotations$ID.y))
  return(list("paralog_data" = paralog_data, "gathered_paralog_data" = gathered_paralog_data, "Total_paralog_annotations" = Total_paralog_annotations, "num_of_paralog_anno" = num_of_paralog_anno, "ref_data" = ref_data, "max_no_col" = max_no_col))
}

Subset_var_align = function(gene_subset, paralogs2_file, paralog_tableized_file, joining_tableized_file=paralog_tableized_file){ #Function for joining together variant paralogous locations
  # system(paste("python /media/nick/Data/Users/N/Documents/PhD/Paralogues/ParalogueAnnotation_personal/src/paralogs_file_remove_last_column.py ", paralogs2_file, sep = ""))
  paralog_data = file(paralogs2_file)
  max_no_col = (max(count.fields(paralog_data, sep = "\t"))-5) #-6 for "Variant_pos", "ID", "Gene", "Ref", "Alt", and "\n"; -5 for above python script
  paralog_data = read.csv(file=paralogs2_file, sep="\t", header=FALSE, col.names = c("Variant_pos", "ID", "Gene", "REF", "ALT", paste("paralog", 1:max_no_col, sep = "")))
  
  paralog_data = dplyr::filter(paralog_data, Gene %in% gene_subset)
  
  paralog_tableized_data = read.csv(file=paralog_tableized_file, sep = "\t", header=TRUE, stringsAsFactors=FALSE)
  paralog_tableized_data = paralog_tableized_data[!is.na(paralog_tableized_data$Amino_acids),]
  
  paralog_tableized_data$Variant_pos = paste(paralog_tableized_data$CHROM,paralog_tableized_data$POS, sep = " ")
  paralog_tableized_data$REF_Amino_acids = sapply(paralog_tableized_data[,"Amino_acids"],strsplit, "/")
  paralog_tableized_data$REF_Amino_acids = sapply(paralog_tableized_data[,"REF_Amino_acids"],unlist)
  paralog_tableized_data$REF_Amino_acids = sapply(paralog_tableized_data[,"REF_Amino_acids"],function(x) x[1])
  paralog_tableized_data$ALT_Amino_acids = sapply(paralog_tableized_data[,"Amino_acids"],strsplit, "/")
  paralog_tableized_data$ALT_Amino_acids = sapply(paralog_tableized_data[,"ALT_Amino_acids"],unlist)
  paralog_tableized_data$ALT_Amino_acids = sapply(paralog_tableized_data[,"ALT_Amino_acids"],function(x) x[2])
  paralog_data = left_join(paralog_data,paralog_tableized_data, by =  c("Variant_pos", "ID", "REF", "ALT", "Gene" = "SYMBOL"))
  gathered_paralog_data = filter(gather(paralog_data, paralog, paralog_pos, paste("paralog", 1:max_no_col, sep = ""), factor_key = TRUE), paralog_pos != "")
  
  joining_tableized_data = read.csv(file=joining_tableized_file, sep = "\t", header=TRUE, stringsAsFactors=FALSE)
  joining_tableized_data = joining_tableized_data[!is.na(joining_tableized_data$Amino_acids),]
  
  joining_tableized_data$Variant_pos = paste(joining_tableized_data$CHROM,joining_tableized_data$POS, sep = " ")
  joining_tableized_data$REF_Amino_acids = sapply(joining_tableized_data[,"Amino_acids"],strsplit, "/")
  joining_tableized_data$REF_Amino_acids = sapply(joining_tableized_data[,"REF_Amino_acids"],unlist)
  joining_tableized_data$REF_Amino_acids = sapply(joining_tableized_data[,"REF_Amino_acids"],function(x) x[1])
  joining_tableized_data$ALT_Amino_acids = sapply(joining_tableized_data[,"Amino_acids"],strsplit, "/")
  joining_tableized_data$ALT_Amino_acids = sapply(joining_tableized_data[,"ALT_Amino_acids"],unlist)
  joining_tableized_data$ALT_Amino_acids = sapply(joining_tableized_data[,"ALT_Amino_acids"],function(x) x[2])
  ref_data = joining_tableized_data
  
  Total_paralog_annotations = left_join(gathered_paralog_data, ref_data, by = c("paralog_pos" = "Variant_pos"))
  num_of_paralog_anno = sum(!is.na(Total_paralog_annotations$ID.y))
  return(list("paralog_data" = paralog_data, "gathered_paralog_data" = gathered_paralog_data, "Total_paralog_annotations" = Total_paralog_annotations, "num_of_paralog_anno" = num_of_paralog_anno, "ref_data" = ref_data, "max_no_col" = max_no_col))
}

ParaZ_var_align = function(paraz_cutoff,paralogs2_file, paralog_tableized_file, joining_tableized_file=paralog_tableized_file){ #Function for joining together variant paralogous locations with Para Z scores
  # system(paste("python /media/nick/Data/Users/N/Documents/PhD/Paralogues/ParalogueAnnotation_personal/src/paralogs_file_remove_last_column.py ", paralogs2_file, sep = ""))
  paralog_data = file(paralogs2_file)
  max_no_col = (max(count.fields(paralog_data, sep = "\t"))-5) #-6 for "Variant_pos", "ID", "Gene", "Ref", "Alt", and "\n"; -5 for above python script
  paralog_data = read.csv(file=paralogs2_file, sep="\t", header=FALSE, col.names = c("Variant_pos", "ID", "Gene", "REF", "ALT", paste("paralog", 1:max_no_col, sep = "")))
  paralog_tableized_data = read.csv(file=paralog_tableized_file, sep = "\t", header=TRUE, stringsAsFactors=FALSE)
  paralog_tableized_data = paralog_tableized_data[!is.na(paralog_tableized_data$Amino_acids),]
  
  paralog_tableized_data$Variant_pos = paste(paralog_tableized_data$CHROM,paralog_tableized_data$POS, sep = " ")
  paralog_tableized_data$REF_Amino_acids = sapply(paralog_tableized_data[,"Amino_acids"],strsplit, "/")
  paralog_tableized_data$REF_Amino_acids = sapply(paralog_tableized_data[,"REF_Amino_acids"],unlist)
  paralog_tableized_data$REF_Amino_acids = sapply(paralog_tableized_data[,"REF_Amino_acids"],function(x) x[1])
  paralog_tableized_data$ALT_Amino_acids = sapply(paralog_tableized_data[,"Amino_acids"],strsplit, "/")
  paralog_tableized_data$ALT_Amino_acids = sapply(paralog_tableized_data[,"ALT_Amino_acids"],unlist)
  paralog_tableized_data$ALT_Amino_acids = sapply(paralog_tableized_data[,"ALT_Amino_acids"],function(x) x[2])
  paralog_data = left_join(paralog_data,paralog_tableized_data, by =  c("Variant_pos", "ID", "REF", "ALT", "Gene" = "SYMBOL"))
  
  paraz_paralog_data = paralog_data[!is.na(paralog_data$Para_Z_score) & paralog_data$Para_Z_score >= paraz_cutoff,]
  
  gathered_paralog_data = filter(gather(paraz_paralog_data, paralog, paralog_pos, paste("paralog", 1:max_no_col, sep = ""), factor_key = TRUE), paralog_pos != "")
  
  joining_tableized_data = read.csv(file=joining_tableized_file, sep = "\t", header=TRUE, stringsAsFactors=FALSE)
  joining_tableized_data = joining_tableized_data[!is.na(joining_tableized_data$Amino_acids),]
  
  joining_tableized_data$Variant_pos = paste(joining_tableized_data$CHROM,joining_tableized_data$POS, sep = " ")
  joining_tableized_data$REF_Amino_acids = sapply(joining_tableized_data[,"Amino_acids"],strsplit, "/")
  joining_tableized_data$REF_Amino_acids = sapply(joining_tableized_data[,"REF_Amino_acids"],unlist)
  joining_tableized_data$REF_Amino_acids = sapply(joining_tableized_data[,"REF_Amino_acids"],function(x) x[1])
  joining_tableized_data$ALT_Amino_acids = sapply(joining_tableized_data[,"Amino_acids"],strsplit, "/")
  joining_tableized_data$ALT_Amino_acids = sapply(joining_tableized_data[,"ALT_Amino_acids"],unlist)
  joining_tableized_data$ALT_Amino_acids = sapply(joining_tableized_data[,"ALT_Amino_acids"],function(x) x[2])
  ref_data = joining_tableized_data
  
  Total_paralog_annotations = left_join(gathered_paralog_data, ref_data, by = c("paralog_pos" = "Variant_pos"))
  num_of_paralog_anno = sum(!is.na(Total_paralog_annotations$ID.y))
  return(list("paralog_data" = paralog_data, "gathered_paralog_data" = gathered_paralog_data, "Total_paralog_annotations" = Total_paralog_annotations, "num_of_paralog_anno" = num_of_paralog_anno, "ref_data" = ref_data, "max_no_col" = max_no_col))
}

conf_matrix = function(ptop.num_of_paralog_anno, p.paralog_data, btop.num_of_paralog_anno, b.paralog_data){ #function for calculating confusion matrix and stats
  con_table_TP = ptop.num_of_paralog_anno
  con_table_FP = btop.num_of_paralog_anno
  con_table_TN = nrow(b.paralog_data)-btop.num_of_paralog_anno
  con_table_FN = nrow(p.paralog_data)-ptop.num_of_paralog_anno
  con_table_PPV = con_table_TP/(con_table_TP+con_table_FP)
  con_table_Sensitivty = con_table_TP/(con_table_TP+con_table_FN)
  con_table_FPR = con_table_FP/(con_table_FP+con_table_TN)
  con_table = matrix(
    c(nrow(p.paralog_data),
      ptop.num_of_paralog_anno,
      nrow(b.paralog_data),
      btop.num_of_paralog_anno
    ), ncol = 2
  )
  colnames(con_table) = c("Pathogenic", "Benign")
  rownames(con_table) = c("Number of variants in total", "Number of variants predicted as pathogenic")
  con_table_p_value = fisher.test(con_table)
  return(list("con_table" = con_table, "PPV" = con_table_PPV, "Sensitivity" = con_table_Sensitivty, "FPR" = con_table_FPR,"Pvalue" = con_table_p_value, "TP" = con_table_TP, "FP" = con_table_FP, "TN" = con_table_TN, "FN" = con_table_FN))
}

conf_matrix_benign = function(ptob.num_of_paralog_anno, p.paralog_data, btob.num_of_paralog_anno, b.paralog_data){ #function for calculating confusion matrix and stats
  con_table_TP = ptob.num_of_paralog_anno
  con_table_FP = btob.num_of_paralog_anno
  con_table_TN = nrow(b.paralog_data)-btob.num_of_paralog_anno
  con_table_FN = nrow(p.paralog_data)-ptob.num_of_paralog_anno
  con_table_PPV = con_table_TP/(con_table_TP+con_table_FP)
  con_table_Sensitivty = con_table_TP/(con_table_TP+con_table_FN)
  con_table_FPR = con_table_FP/(con_table_FP+con_table_TN)
  con_table = matrix(
    c(nrow(p.paralog_data),
      ptob.num_of_paralog_anno,
      nrow(b.paralog_data),
      btob.num_of_paralog_anno
    ), ncol = 2
  )
  colnames(con_table) = c("Pathogenic", "Benign")
  rownames(con_table) = c("Number of variants in total", "Number of variants predicted as benign")
  con_table_p_value = fisher.test(con_table)
  return(list("con_table" = con_table, "PPV" = con_table_PPV, "Sensitivity" = con_table_Sensitivty, "FPR" = con_table_FPR,"Pvalue" = con_table_p_value, "TP" = con_table_TP, "FP" = con_table_FP, "TN" = con_table_TN, "FN" = con_table_FN))
}

var_rem_matrix = function(con_table1, con_table2, p.paralog_data, b.paralog_data){ #function for calculating confusion matrix and stats for variants removed after filtering at each QC stage
  var_rem_con_TP = con_table1[2]-con_table2[2]
  var_rem_con_FP = con_table1[4]-con_table2[4]
  var_rem_con_FN = length(p.paralog_data$ID)-var_rem_con_TP
  var_rem_con_PPV = var_rem_con_TP/(var_rem_con_TP+var_rem_con_FP)
  var_rem_con_Sensitivity = var_rem_con_TP/(var_rem_con_TP+var_rem_con_FN)
  var_rem_con = matrix(
    c(length(p.paralog_data$ID),
      var_rem_con_TP,
      length(b.paralog_data$ID),
      var_rem_con_FP
    ), ncol = 2
  )
  colnames(var_rem_con) = c("Pathogenic", "Benign")
  rownames(var_rem_con) = c("Number of variants in total", "Number of variants predicted as pathogenic")
  var_rem_con_p_value = fisher.test(var_rem_con)
  return(list("con_table" = var_rem_con, "PPV" = var_rem_con_PPV, "Sensitivity" = var_rem_con_Sensitivity, "Pvalue" = var_rem_con_p_value, "TP" = var_rem_con_TP, "FP" = var_rem_con_FP, "FN" = var_rem_con_FN))
}