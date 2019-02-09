Packages = c("tidyverse", "dplyr", "ggplot2", "ggsignif", "biomaRt", "knitr", "png", "grid", "tinytex", "pander", "kableExtra", "clusterProfiler", "org.Hs.eg.db", "DiagrammeR")
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
  input_variants = distinct(paralog_data[c("Variant_pos", "ID", "REF", "ALT")])
  
  if (Overlap == 0){
    n_occur = data.frame(table(paralog_data$ID))
    n_occur = n_occur[n_occur$Freq==1,]
    paralog_data = filter(paralog_data, ID %in% n_occur$Var1)
  } else if (Overlap == 1){
    n_occur = data.frame(table(paralog_data$ID))
    n_occur = n_occur[n_occur$Freq>1,]
    paralog_data = filter(paralog_data, ID %in% n_occur$Var1)
  } else {
    n_occur = NA
  }
  
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
  ref_data = joining_tableized_data
  
  Total_paralog_annotations = left_join(gathered_paralog_data, ref_data, by = c("paralog_pos" = "Variant_pos"))
  Total_paralog_annotations = distinct(Total_paralog_annotations)
  Unique_variant_annotations = Total_paralog_annotations[!is.na(Total_paralog_annotations$ID.y),]
  Unique_variant_annotations = distinct(Unique_variant_annotations[c("Variant_pos", "ID.x", "REF.x", "ALT.x")])
  
  true_num_of_paralog_anno = length(Unique_variant_annotations$Variant_pos)
  num_of_paralog_anno = sum(!is.na(Total_paralog_annotations$ID.y))
  return(list("paralog_data" = paralog_data, 
              "input_variants" = input_variants, 
              "Unique_variant_annotations" = Unique_variant_annotations,
              "gathered_paralog_data" = gathered_paralog_data, 
              "Total_paralog_annotations" = Total_paralog_annotations, 
              "num_of_paralog_anno" = num_of_paralog_anno, 
              "true_num_of_paralog_anno" = true_num_of_paralog_anno, 
              "ref_data" = ref_data, 
              "max_no_col" = max_no_col,
              "Subset" = Subset,
              "n_occur" = n_occur))
}

# Paralogous_var_align_compressed is depreciated after implementing features into Paralogous_var_align
Paralogous_var_align_compressed = function(paralogs2_file, paralog_tableized_file, joining_tableized_file=paralog_tableized_file){ #Function for joining together variant paralogous locations - compressed files
  paralog_data = gzfile(paralogs2_file)
  max_no_col = (max(count.fields(paralog_data, sep = "\t"))-5) #-6 for "Variant_pos", "ID", "Gene", "Ref", "Alt", and "\n"; -5 for above python script
  paralog_data = read.csv(file=gzfile(paralogs2_file), sep="\t", header=FALSE, col.names = c("Variant_pos", "ID", "Gene", "REF", "ALT", paste("paralog", 1:max_no_col, sep = "")))
  paralog_tableized_data = read.csv(file=gzfile(paralog_tableized_file), sep = "\t", header=TRUE, stringsAsFactors=FALSE)
  paralog_tableized_data = paralog_tableized_data[!is.na(paralog_tableized_data$Amino_acids),]
  
  paralog_tableized_data$Variant_pos = paste(paralog_tableized_data$CHROM,paralog_tableized_data$POS, sep = " ")
  paralog_tableized_data$REF_Amino_acids = sapply(paralog_tableized_data[,"Amino_acids"],strsplit, "/")
  paralog_tableized_data$REF_Amino_acids = sapply(paralog_tableized_data[,"REF_Amino_acids"],function(x) x[1])
  paralog_tableized_data$ALT_Amino_acids = sapply(paralog_tableized_data[,"Amino_acids"],strsplit, "/")
  paralog_tableized_data$ALT_Amino_acids = sapply(paralog_tableized_data[,"ALT_Amino_acids"],function(x) x[2])
  paralog_data = left_join(paralog_data,paralog_tableized_data, by =  c("Variant_pos", "ID", "REF", "ALT", "Gene" = "SYMBOL"))
  
  paralog_data = distinct(paralog_data)
  
  gathered_paralog_data = filter(gather(paralog_data, paralog, paralog_pos, paste("paralog", 1:max_no_col, sep = ""), factor_key = TRUE), paralog_pos != "")
  
  joining_tableized_data = read.csv(file=joining_tableized_file, sep = "\t", header=TRUE, stringsAsFactors=FALSE)
  joining_tableized_data = joining_tableized_data[!is.na(joining_tableized_data$Amino_acids),]
  
  joining_tableized_data$Variant_pos = paste(joining_tableized_data$CHROM,joining_tableized_data$POS, sep = " ")
  joining_tableized_data$REF_Amino_acids = sapply(joining_tableized_data[,"Amino_acids"],strsplit, "/")
  joining_tableized_data$REF_Amino_acids = sapply(joining_tableized_data[,"REF_Amino_acids"],function(x) x[1])
  joining_tableized_data$ALT_Amino_acids = sapply(joining_tableized_data[,"Amino_acids"],strsplit, "/")
  joining_tableized_data$ALT_Amino_acids = sapply(joining_tableized_data[,"ALT_Amino_acids"],function(x) x[2])
  ref_data = joining_tableized_data
  
  Total_paralog_annotations = left_join(gathered_paralog_data, ref_data, by = c("paralog_pos" = "Variant_pos"))
  Total_paralog_annotations = distinct(Total_paralog_annotations)
  
  num_of_paralog_anno = sum(!is.na(Total_paralog_annotations$ID.y))
  return(list("paralog_data" = paralog_data, "gathered_paralog_data" = gathered_paralog_data, "Total_paralog_annotations" = Total_paralog_annotations, "num_of_paralog_anno" = num_of_paralog_anno, "ref_data" = ref_data, "max_no_col" = max_no_col))
}
#

# Paralogous_var_align_no_overlap and Paralogous_var_align_overlap are depreciated after implementing features into Paralogous_var_align
Paralogous_var_align_no_overlap = function(paralogs2_file, paralog_tableized_file, joining_tableized_file=paralog_tableized_file){ #Function for joining together variant paralogous locations, BUT ignores variants belonging to overlapping genes
  paralog_data = file(paralogs2_file)
  max_no_col = (max(count.fields(paralog_data, sep = "\t"))-5) #-6 for "Variant_pos", "ID", "Gene", "Ref", "Alt", and "\n"; -5 for above python script
  paralog_data = read.csv(file=paralogs2_file, sep="\t", header=FALSE, col.names = c("Variant_pos", "ID", "Gene", "REF", "ALT", paste("paralog", 1:max_no_col, sep = "")))
  paralog_tableized_data = read.csv(file=paralog_tableized_file, sep = "\t", header=TRUE, stringsAsFactors=FALSE)
  paralog_tableized_data = paralog_tableized_data[!is.na(paralog_tableized_data$Amino_acids),]
  
  paralog_tableized_data$Variant_pos = paste(paralog_tableized_data$CHROM,paralog_tableized_data$POS, sep = " ")
  paralog_tableized_data$REF_Amino_acids = sapply(paralog_tableized_data[,"Amino_acids"],strsplit, "/")
  paralog_tableized_data$REF_Amino_acids = sapply(paralog_tableized_data[,"REF_Amino_acids"],function(x) x[1])
  paralog_tableized_data$ALT_Amino_acids = sapply(paralog_tableized_data[,"Amino_acids"],strsplit, "/")
  paralog_tableized_data$ALT_Amino_acids = sapply(paralog_tableized_data[,"ALT_Amino_acids"],function(x) x[2])
  paralog_data = left_join(paralog_data,paralog_tableized_data, by =  c("Variant_pos", "ID", "REF", "ALT", "Gene" = "SYMBOL"))
  
  paralog_data = paralog_data[!duplicated(paralog_data),]
  
  n_occur = data.frame(table(paralog_data$ID))
  n_occur = n_occur[n_occur$Freq==1,]
  paralog_data = filter(paralog_data, ID %in% n_occur$Var1)
  
  gathered_paralog_data = filter(gather(paralog_data, paralog, paralog_pos, paste("paralog", 1:max_no_col, sep = ""), factor_key = TRUE), paralog_pos != "")
  
  joining_tableized_data = read.csv(file=joining_tableized_file, sep = "\t", header=TRUE, stringsAsFactors=FALSE)
  joining_tableized_data = joining_tableized_data[!is.na(joining_tableized_data$Amino_acids),]
  
  joining_tableized_data$Variant_pos = paste(joining_tableized_data$CHROM,joining_tableized_data$POS, sep = " ")
  joining_tableized_data$REF_Amino_acids = sapply(joining_tableized_data[,"Amino_acids"],strsplit, "/")
  joining_tableized_data$REF_Amino_acids = sapply(joining_tableized_data[,"REF_Amino_acids"],function(x) x[1])
  joining_tableized_data$ALT_Amino_acids = sapply(joining_tableized_data[,"Amino_acids"],strsplit, "/")
  joining_tableized_data$ALT_Amino_acids = sapply(joining_tableized_data[,"ALT_Amino_acids"],function(x) x[2])
  ref_data = joining_tableized_data
  
  Total_paralog_annotations = left_join(gathered_paralog_data, ref_data, by = c("paralog_pos" = "Variant_pos"))
  num_of_paralog_anno = sum(!is.na(Total_paralog_annotations$ID.y))
  return(list("paralog_data" = paralog_data, "gathered_paralog_data" = gathered_paralog_data, "Total_paralog_annotations" = Total_paralog_annotations, "num_of_paralog_anno" = num_of_paralog_anno, "ref_data" = ref_data, "max_no_col" = max_no_col, "n_occur" = n_occur))
}
Paralogous_var_align_overlap = function(paralogs2_file, paralog_tableized_file, joining_tableized_file=paralog_tableized_file){ #Function for joining together variant paralogous locations, BUT ignores variants belonging to overlapping genes
  paralog_data = file(paralogs2_file)
  max_no_col = (max(count.fields(paralog_data, sep = "\t"))-5) #-6 for "Variant_pos", "ID", "Gene", "Ref", "Alt", and "\n"; -5 for above python script
  paralog_data = read.csv(file=paralogs2_file, sep="\t", header=FALSE, col.names = c("Variant_pos", "ID", "Gene", "REF", "ALT", paste("paralog", 1:max_no_col, sep = "")))
  paralog_tableized_data = read.csv(file=paralog_tableized_file, sep = "\t", header=TRUE, stringsAsFactors=FALSE)
  paralog_tableized_data = paralog_tableized_data[!is.na(paralog_tableized_data$Amino_acids),]
  
  paralog_tableized_data$Variant_pos = paste(paralog_tableized_data$CHROM,paralog_tableized_data$POS, sep = " ")
  paralog_tableized_data$REF_Amino_acids = sapply(paralog_tableized_data[,"Amino_acids"],strsplit, "/")
  paralog_tableized_data$REF_Amino_acids = sapply(paralog_tableized_data[,"REF_Amino_acids"],function(x) x[1])
  paralog_tableized_data$ALT_Amino_acids = sapply(paralog_tableized_data[,"Amino_acids"],strsplit, "/")
  paralog_tableized_data$ALT_Amino_acids = sapply(paralog_tableized_data[,"ALT_Amino_acids"],function(x) x[2])
  paralog_data = left_join(paralog_data,paralog_tableized_data, by =  c("Variant_pos", "ID", "REF", "ALT", "Gene" = "SYMBOL"))
  
  paralog_data = paralog_data[!duplicated(paralog_data),]
  
  n_occur = data.frame(table(paralog_data$ID))
  n_occur = n_occur[n_occur$Freq>1,]
  paralog_data = filter(paralog_data, ID %in% n_occur$Var1)
  
  gathered_paralog_data = filter(gather(paralog_data, paralog, paralog_pos, paste("paralog", 1:max_no_col, sep = ""), factor_key = TRUE), paralog_pos != "")
  
  joining_tableized_data = read.csv(file=joining_tableized_file, sep = "\t", header=TRUE, stringsAsFactors=FALSE)
  joining_tableized_data = joining_tableized_data[!is.na(joining_tableized_data$Amino_acids),]
  
  joining_tableized_data$Variant_pos = paste(joining_tableized_data$CHROM,joining_tableized_data$POS, sep = " ")
  joining_tableized_data$REF_Amino_acids = sapply(joining_tableized_data[,"Amino_acids"],strsplit, "/")
  joining_tableized_data$REF_Amino_acids = sapply(joining_tableized_data[,"REF_Amino_acids"],function(x) x[1])
  joining_tableized_data$ALT_Amino_acids = sapply(joining_tableized_data[,"Amino_acids"],strsplit, "/")
  joining_tableized_data$ALT_Amino_acids = sapply(joining_tableized_data[,"ALT_Amino_acids"],function(x) x[2])
  ref_data = joining_tableized_data
  
  Total_paralog_annotations = left_join(gathered_paralog_data, ref_data, by = c("paralog_pos" = "Variant_pos"))
  num_of_paralog_anno = sum(!is.na(Total_paralog_annotations$ID.y))
  return(list("paralog_data" = paralog_data, "gathered_paralog_data" = gathered_paralog_data, "Total_paralog_annotations" = Total_paralog_annotations, "num_of_paralog_anno" = num_of_paralog_anno, "ref_data" = ref_data, "max_no_col" = max_no_col, "n_occur" = n_occur))
}
#

# Subset_var_align is depreciated after implementing features into Paralogous_var_align
Subset_var_align_gene = function(gene_subset, paralogs2_file, paralog_tableized_file, joining_tableized_file=paralog_tableized_file){ #Function for joining together variant paralogous locations
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
Subset_var_align_id = function(id_subset, paralogs2_file, paralog_tableized_file, joining_tableized_file=paralog_tableized_file){ #Function for joining together variant paralogous locations
  # system(paste("python /media/nick/Data/Users/N/Documents/PhD/Paralogues/ParalogueAnnotation_personal/src/paralogs_file_remove_last_column.py ", paralogs2_file, sep = ""))
  paralog_data = file(paralogs2_file)
  max_no_col = (max(count.fields(paralog_data, sep = "\t"))-5) #-6 for "Variant_pos", "ID", "Gene", "Ref", "Alt", and "\n"; -5 for above python script
  paralog_data = read.csv(file=paralogs2_file, sep="\t", header=FALSE, col.names = c("Variant_pos", "ID", "Gene", "REF", "ALT", paste("paralog", 1:max_no_col, sep = "")))
  
  paralog_data = dplyr::filter(paralog_data, ID %in% id_subset)
  
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
#

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
  con_table_Specificity = con_table_TN/(con_table_TN+con_table_FP)
  con_table_FPR = con_table_FP/(con_table_FP+con_table_TN)
  con_table_ACC = (con_table_TP + con_table_TN)/(con_table_TP + con_table_TN + con_table_FP + con_table_FN)
  con_table = matrix(
    c(nrow(p.paralog_data)-ptop.num_of_paralog_anno,
      ptop.num_of_paralog_anno,
      nrow(b.paralog_data)-btop.num_of_paralog_anno,
      btop.num_of_paralog_anno
    ), ncol = 2
  )
  colnames(con_table) = c("Pathogenic", "Benign")
  rownames(con_table) = c("Number of variants not predicted as pathogenic", "Number of variants predicted as pathogenic")
  con_table_p_value = fisher.test(con_table)
  return(list("con_table" = con_table, "Accuracy" = con_table_ACC, "PPV" = con_table_PPV, "Sensitivity" = con_table_Sensitivty, "Specificity" = con_table_Specificity, "FPR" = con_table_FPR, "Pvalue" = con_table_p_value, "TP" = con_table_TP, "FP" = con_table_FP, "TN" = con_table_TN, "FN" = con_table_FN))
}

conf_matrix_benign = function(ptob.num_of_paralog_anno, p.paralog_data, btob.num_of_paralog_anno, b.paralog_data){ #function for calculating confusion matrix and stats
  con_table_TP = btob.num_of_paralog_anno
  con_table_FP = ptob.num_of_paralog_anno
  con_table_FN = nrow(b.paralog_data)-btob.num_of_paralog_anno
  con_table_TN = nrow(p.paralog_data)-ptob.num_of_paralog_anno
  con_table_PPV = con_table_TP/(con_table_TP+con_table_FP)
  con_table_Sensitivty = con_table_TP/(con_table_TP+con_table_FN)
  con_table_Specificity = con_table_TN/(con_table_TN+con_table_FP)
  con_table_FPR = con_table_FP/(con_table_FP+con_table_TN)
  con_table_ACC = (con_table_TP + con_table_TN)/(con_table_TP + con_table_TN + con_table_FP + con_table_FN)
  con_table = matrix(
    c(nrow(p.paralog_data)-ptob.num_of_paralog_anno,
      ptob.num_of_paralog_anno,
      nrow(b.paralog_data)-btob.num_of_paralog_anno,
      btob.num_of_paralog_anno
    ), ncol = 2
  )
  colnames(con_table) = c("Pathogenic", "Benign")
  rownames(con_table) = c("Number of variants not predicted as benign", "Number of variants predicted as benign")
  con_table_p_value = fisher.test(con_table)
  return(list("con_table" = con_table, "Accuracy" = con_table_ACC, "PPV" = con_table_PPV, "Sensitivity" = con_table_Sensitivty, "Specificity" = con_table_Specificity, "FPR" = con_table_FPR,"Pvalue" = con_table_p_value, "TP" = con_table_TP, "FP" = con_table_FP, "TN" = con_table_TN, "FN" = con_table_FN))
}

var_rem_matrix = function(con_table1, con_table2, p.paralog_data, b.paralog_data){ #function for calculating confusion matrix and stats for variants removed after filtering at each QC stage
  var_rem_con_TP = con_table1[2]-con_table2[2]
  var_rem_con_FP = con_table1[4]-con_table2[4]
  var_rem_con_TN = length(b.paralog_data$ID)-var_rem_con_FP
  var_rem_con_FN = length(p.paralog_data$ID)-var_rem_con_TP
  var_rem_con_PPV = var_rem_con_TP/(var_rem_con_TP+var_rem_con_FP)
  var_rem_con_Sensitivity = var_rem_con_TP/(var_rem_con_TP+var_rem_con_FN)
  var_rem_con_Specificity = var_rem_con_TN/(var_rem_con_TN+var_rem_con_FP)
  var_rem_con_FPR = var_rem_con_FP/(var_rem_con_FP+var_rem_con_TN)
  var_rem_con_ACC = (var_rem_con_TP + var_rem_con_TN)/(var_rem_con_TP + var_rem_con_TN + var_rem_con_FP + var_rem_con_FN)
  var_rem_con = matrix(
    c(nrow(p.paralog_data)-var_rem_con_TP,
      var_rem_con_TP,
      nrow(b.paralog_data)-var_rem_con_FP,
      var_rem_con_FP
    ), ncol = 2
  )
  colnames(var_rem_con) = c("Pathogenic", "Benign")
  rownames(var_rem_con) = c("Number of variants in total", "Number of variants predicted as pathogenic")
  var_rem_con_p_value = fisher.test(var_rem_con)
  return(list("con_table" = var_rem_con, "Accuracy" = var_rem_con_ACC, "PPV" = var_rem_con_PPV, "Sensitivity" = var_rem_con_Sensitivity, "Specificity" = var_rem_con_Specificity, "FPR" = var_rem_con_FPR, "Pvalue" = var_rem_con_p_value, "TP" = var_rem_con_TP, "FP" = var_rem_con_FP, "FN" = var_rem_con_FN))
}

var_rem_matrix_benign = function(con_table1, con_table2, p.paralog_data, b.paralog_data){ #function for calculating confusion matrix and stats for variants removed after filtering at each QC stage
  var_rem_con_TP = con_table1[4]-con_table2[4]
  var_rem_con_FP = con_table1[2]-con_table2[2]
  var_rem_con_TN = length(p.paralog_data$ID)-var_rem_con_FP
  var_rem_con_FN = length(b.paralog_data$ID)-var_rem_con_TP
  var_rem_con_PPV = var_rem_con_TP/(var_rem_con_TP+var_rem_con_FP)
  var_rem_con_Sensitivity = var_rem_con_TP/(var_rem_con_TP+var_rem_con_FN)
  var_rem_con_Specificity = var_rem_con_TN/(var_rem_con_TN+var_rem_con_FP)
  var_rem_con_FPR = var_rem_con_FP/(var_rem_con_FP+var_rem_con_TN)
  var_rem_con_ACC = (var_rem_con_TP + var_rem_con_TN)/(var_rem_con_TP + var_rem_con_TN + var_rem_con_FP + var_rem_con_FN)
  var_rem_con = matrix(
    c(nrow(p.paralog_data)-var_rem_con_TP,
      var_rem_con_TP,
      nrow(b.paralog_data)-var_rem_con_FP,
      var_rem_con_FP
    ), ncol = 2
  )
  colnames(var_rem_con) = c("Pathogenic", "Benign")
  rownames(var_rem_con) = c("Number of variants in total", "Number of variants predicted as pathogenic")
  var_rem_con_p_value = fisher.test(var_rem_con)
  return(list("con_table" = var_rem_con, "Accuracy" = var_rem_con_ACC, "PPV" = var_rem_con_PPV, "Sensitivity" = var_rem_con_Sensitivity, "Specificity" = var_rem_con_Specificity, "FPR" = var_rem_con_FPR, "Pvalue" = var_rem_con_p_value, "TP" = var_rem_con_TP, "FP" = var_rem_con_FP, "FN" = var_rem_con_FN))
}


raw_conf_matrix = function(num_patho_pred_patho, p.tableized_data, num_benign_pred_patho, b.tableized_data){ #function for calculating confusion matrix and stats from raw TP, TN, FP and FN numbers from anything
  con_table_TP = num_patho_pred_patho
  con_table_FP = num_benign_pred_patho
  con_table_TN = nrow(b.tableized_data)-num_benign_pred_patho
  con_table_FN = nrow(p.tableized_data)-num_patho_pred_patho
  con_table_PPV = con_table_TP/(con_table_TP+con_table_FP)
  con_table_Sensitivty = con_table_TP/(con_table_TP+con_table_FN)
  con_table_Specificity = con_table_TN/(con_table_TN+con_table_FP)
  con_table_FPR = con_table_FP/(con_table_FP+con_table_TN)
  con_table_ACC = (con_table_TP + con_table_TN)/(con_table_TP + con_table_TN + con_table_FP + con_table_FN)
  con_table = matrix(
    c(nrow(p.tableized_data),
      num_patho_pred_patho,
      nrow(b.tableized_data),
      num_benign_pred_patho
    ), ncol = 2
  )
  colnames(con_table) = c("Pathogenic", "Benign")
  rownames(con_table) = c("Number of variants in total", "Number of variants predicted as pathogenic")
  con_table_p_value = fisher.test(con_table)
  return(list("con_table" = con_table, "Accuracy" = con_table_ACC, "PPV" = con_table_PPV, "Sensitivity" = con_table_Sensitivty, "Specificity" = con_table_Specificity, "FPR" = con_table_FPR,"Pvalue" = con_table_p_value, "TP" = con_table_TP, "FP" = con_table_FP, "TN" = con_table_TN, "FN" = con_table_FN))
}

calc_EF = function(a, b, c, d){ #function for calculating Odds Ratios and Etiological Fractions
  if (a == 0 | b == 0 | c == 0 | d == 0){
    # a = a + 0.5
    # b = b + 0.5
    # c = c + 0.5
    # d = d + 0.5
    OR = "NA"
    OR_CI = "NA"
    EF = "NA"
    EF_CI = "NA"
  } else {
  
  OR = (a/b)/(c/d)
  
  SE_ln_OR = sqrt((1/a)+(1/b)+(1/c)+(1/d))
  
  OR_CI = c(
    exp(log(OR) - 1.96 * SE_ln_OR),
    exp(log(OR) + 1.96 * SE_ln_OR)
    )
  
  EF = (OR-1)/OR
  
  pi_hat_1 = a/(b + a)
  pi_hat_0 = c/(d + c)
  N_0 = c + d
  N_1 = a + b
  EF_hat = 1 - (pi_hat_0/pi_hat_1)
  Phi_hat_squared = (pi_hat_0/pi_hat_1)^2
  VAR_hat_Phi_hat = Phi_hat_squared * (((1-pi_hat_0)/(N_0 * pi_hat_0)) + ((1-pi_hat_1)/(N_1 * pi_hat_1)))
  
  EF_CI = c(
    EF_hat - 1.96 * sqrt(VAR_hat_Phi_hat),
    min(EF_hat + 1.96 * sqrt(VAR_hat_Phi_hat), 1)
  )
  }
  return(list(
              # "VAR_hat_Phi_hat" = VAR_hat_Phi_hat,
              # "Phi_hat_squared" = Phi_hat_squared,
              # "EF_hat" = EF_hat,
              # "pi_hat_1" = pi_hat_1,
              # "pi_hat_0" = pi_hat_0,
              # "N_0" = N_0,
              # "N_1" = N_1,
              "OR" = OR, 
              "OR_CI" = OR_CI, 
              "EF" = EF, 
              "EF_CI" = EF_CI
              ))
}

case_control_gene_split = function(data_df, patho_var_ids, gene){ #function for splitting genes in the case control study to allow gene by gene analysis of EFs
  if (!missing(gene)){
    gene_data = data_df[data_df$gene_symbol==gene,]
  } else {
    gene_data = data_df
  }
  patho_gene_data = gene_data[gene_data$mut_id %in% patho_var_ids,]
  if ("sum.vs_case_count." %in% colnames(data_df)){
    patho_gene_cases = sum(patho_gene_data$sum.vs_case_count.)
    Total_gene_cases = sum(gene_data$sum.vs_case_count.)
  } else {
    patho_gene_cases = sum(patho_gene_data$mut_exac_count)
    Total_gene_cases = sum(gene_data$mut_exac_count)
  }
  return(list("pa_predicted_gene_cases" = patho_gene_cases, "Total_gene_cases" = Total_gene_cases, "gene_data" = gene_data, "pa_predicted_gene_data" = patho_gene_data))
}