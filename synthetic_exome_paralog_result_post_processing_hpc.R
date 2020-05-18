##RUN ON HPC - for trimming data down to only what's needed for paralog app website and also extracting predicted pathogenic variants for GeL##
## will need to export PBS_ARRAY_INDEX and bash for loop to run chromosomes in parallel##

args = commandArgs(trailingOnly=TRUE)

Packages = c("tidyverse", "plyr", "dplyr", "stringr","httr","jsonlite","xml2")
new.packages = Packages[!(Packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "https://cran.ma.imperial.ac.uk/")
# lapply(Packages, library, character.only = TRUE)
library("plyr")
library("dplyr")
# library("tidyverse")

##Split job by chrom##
for (j in args[1]){
        for (k in c("noQC","para_con","all_con")){
                Total_annotations = NULL
                Paraloc = NULL
                files = list.files(path=paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RData_objects/chrom_",j,"/",k), pattern="*.RData", full.names=TRUE, recursive=FALSE)
                for (i in files){
                        print(i)
                        load(i)
                        tmp_Total_annotations = p.normal_PA$Total_paralog_annotations
                        tmp_Total_annotations = tmp_Total_annotations[,!(names(tmp_Total_annotations) %in% c("Variant_pos","FILTER.x","BIOTYPE.x","Paralogue_Vars.x","REF_Amino_acids.x","ALT_Amino_acids.x","paralog","paralog_pos","FILTER.y","BIOTYPE.y","Paralogue_Vars.y","REF_Amino_acids.y","ALT_Amino_acids.y"))]
                        tmp_Total_annotations = tmp_Total_annotations[!is.na(tmp_Total_annotations$POS.x),]
                        tmp_Total_annotations$Para_Z_score.x = as.numeric(tmp_Total_annotations$Para_Z_score.x)
                        tmp_Total_annotations$Para_Z_score.y = as.numeric(tmp_Total_annotations$Para_Z_score.y)
                        
                        split_number = unlist(strsplit(i, "_split", fixed = TRUE))[length(unlist(strsplit(i, "_split", fixed = TRUE)))] #split to get split number of file
                        split_number = unlist(strsplit(split_number, ".", fixed = TRUE))[1]
                        new_tableized_file = read.csv(file=paste0("/rds/general/project/lms-ware-analysis/live/nick/RBH-work/Paralog_Anno/data_files/all_possible_mutation/synthetic_exome/synthetic_exome_GRCh37_renamed/chrom_",j,"/synthetic.vep.cov.table_chrom",j,"_wIDs_proper_split",split_number,".out_paraloc_tableized_for_shinyapp"), sep = "\t", header=TRUE, stringsAsFactors=FALSE)
                        
                        
                        if (is.null(Total_annotations)){
                                Total_annotations = tmp_Total_annotations
                        } else {
                                Total_annotations = base::rbind(Total_annotations, dplyr::setdiff(tmp_Total_annotations, Total_annotations))
                        }
                        if (k == "noQC"){
                                tmp_Paraloc = p.normal_PA$paralog_data
                                # tmp_Paraloc = tmp_Paraloc[,!(names(tmp_Paraloc) %in% c("Variant_pos","FILTER","BIOTYPE","Paralogue_Vars","Para_Z_score","REF_Amino_acids","ALT_Amino_acids","Protein_position","Amino_acids","Codons"))]
                                # tmp_Paraloc = dplyr::select(tmp_Paraloc, CHROM, POS, dplyr::everything())
                                tmp_Paraloc = tmp_Paraloc[,names(tmp_Paraloc) %in% c("CHROM","POS","REF","Gene","Paralogue_Vars")]
                                #tmp_Paraloc = subset(tmp_Paraloc, select=c(CHROM,POS,REF,Gene,Paralogue_Vars)) #IF NOT COMBINING INTO VAR THEN NEED TO CHANGE HOW WE LOOK UP DATA
                                tmp_Paraloc = tmp_Paraloc[!is.na(tmp_Paraloc$POS),]
                                tmp_Paraloc = dplyr::distinct(tmp_Paraloc)
                                tmp_Paraloc$Paralogue_Vars = as.character(tmp_Paraloc$Paralogue_Vars)
                                tmp_Paraloc$Paralogue_Vars = sapply(tmp_Paraloc$Paralogue_Vars, stringr::str_replace, "&", "") #PROBABLY A GOOD IDEA TO DO THIS IN POST-PROCESSING BEFORE LOADING DATA IN 
                                tmp_Paraloc$Paralogue_Vars = sapply(tmp_Paraloc$Paralogue_Vars, stringr::str_replace_all, "&", " ")
                                
                                
                                tmp_Paraloc$CHROM = as.character(tmp_Paraloc$CHROM)
                                # if (all(is.na(tmp_Paraloc$POS))){
                                #         tmp_df = data.frame(QC = character())
                                #         tmp_Paraloc = base::cbind(tmp_Paraloc, tmp_df)
                                # } else {
                                #         tmp_Paraloc$QC = toString(k)
                                # }
                                if (is.null(Paraloc)){
                                        Paraloc = tmp_Paraloc
                                } else {
                                        # Paraloc = plyr::rbind.fill(Paraloc, dplyr::anti_join(x=tmp_Paraloc, y=Paraloc, by=c("CHROM","POS","ID","Gene","REF","ALT")))
                                        Paraloc = base::rbind(Paraloc, dplyr::setdiff(tmp_Paraloc, Paraloc))
                                }
                        }
                }
                
                
                save(Total_annotations, file = paste0("/work/nyl112/data/synthetic_exome/Synthetic_exome_paralog_result_post_processed/chrom_",j,"/Total_annotations_chrom_",j,"_",k,".RData"))
                Total_annotations = Total_annotations[,c("CHROM.x","POS.x","ID.x","REF.x","ALT.x")]
                names(Total_annotations) = c("CHROM","POS","ID","REF","ALT")
                write.table(Total_annotations,file=paste0("/work/nyl112/data/synthetic_exome/Synthetic_exome_paralog_result_post_processed/chrom_",j,"/Total_annotations_chrom_",j,"_",k,"_predicted_pathogenic.vcf"), na="", row.names = FALSE, col.names = TRUE, sep = ",", quote=FALSE)
                if (!is.null(Paraloc)){
                        save(Paraloc, file = paste0("/work/nyl112/data/synthetic_exome/Synthetic_exome_paralog_result_post_processed/chrom_",j,"/Para_locations_chrom_",j,"_",k,".RData"))
                }
        }
}

##Using Ensembl's REST API to create flat file hash tables for ENSG to ENSP ids.
mart_export = read.delim("/work/nyl112/Paralogue_Annotation_App/paralog_app/data/mart_export.txt", quote="", stringsAsFactors=F)
map=setNames(mart_export$HGNC.symbol, mart_export$Gene.stable.ID)
for (j in c(1:22,"X","Y")){
  files = list.files(path=paste0("/rds/general/project/lms-ware-analysis/live/nick/RBH-work/Paralog_Anno/data_files/all_possible_mutation/synthetic_exome/synthetic_exome_GRCh37_renamed/chrom_",j), pattern="*.out_paraloc_tableized_for_shinyapp", full.names=TRUE, recursive=FALSE)
  for (i in files){
    tableized_for_shinyapp = read.csv(file=i, sep = "\t", header=TRUE, stringsAsFactors=FALSE)
  }
}