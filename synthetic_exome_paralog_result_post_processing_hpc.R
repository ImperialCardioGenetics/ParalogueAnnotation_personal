##RUN ON HPC - for trimming data down to only what's needed for paralog app website and also extracting predicted pathogenic variants for GeL##
## will need to export PBS_ARRAY_INDEX and bash for loop to run chromosomes in parallel##

args = commandArgs(trailingOnly=TRUE)

Packages = c("tidyverse", "plyr", "dplyr")
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
                        # if (all(is.na(tmp_Total_annotations$POS.x))){
                        #         tmp_df = data.frame(QC = character())
                        #         tmp_Total_annotations = base::cbind(tmp_Total_annotations, tmp_df)
                        # } else {
                        #         tmp_Total_annotations$QC = toString(k)
                        # }
                        if (is.null(Total_annotations)){
                                Total_annotations = tmp_Total_annotations
                        } else {
                                Total_annotations = base::rbind(Total_annotations, dplyr::setdiff(tmp_Total_annotations, Total_annotations))
                        }
                        if (k == "noQC"){
                                tmp_Paraloc = p.normal_PA$paralog_data
                                tmp_Paraloc = tmp_Paraloc[,!(names(tmp_Paraloc) %in% c("Variant_pos","FILTER","BIOTYPE","Paralogue_Vars","Para_Z_score","REF_Amino_acids","ALT_Amino_acids","Protein_position","Amino_acids","Codons"))]
                                # tmp_Paraloc = tmp_Paraloc[,names(tmp_Paraloc) %in% c("CHROM","POS","ID","Gene","REF","ALT","Paralogue_Vars")]
                                tmp_Paraloc = tmp_Paraloc[!is.na(tmp_Paraloc$POS),]
                                tmp_Paraloc = dplyr::select(tmp_Paraloc, CHROM, POS, dplyr::everything())
                                # if (all(is.na(tmp_Paraloc$POS))){
                                #         tmp_df = data.frame(QC = character())
                                #         tmp_Paraloc = base::cbind(tmp_Paraloc, tmp_df)
                                # } else {
                                #         tmp_Paraloc$QC = toString(k)
                                # }
                                if (is.null(Paraloc)){
                                        Paraloc = tmp_Paraloc
                                } else {
                                        Paraloc = plyr::rbind.fill(Paraloc, dplyr::anti_join(x=tmp_Paraloc, y=Paraloc, by=c("CHROM","POS","ID","Gene","REF","ALT")))
                                        # Paraloc = base::rbind(Paraloc, dplyr::setdiff(tmp_Paraloc, Paraloc))
                                }
                        }
                }
                save(Total_annotations, file = paste0("/work/nyl112/data/synthetic_exome/Synthetic_exome_paralog_result_post_processed/Total_annotations_chrom_",j,"_",k,".RData"))
                Total_annotations = Total_annotations[,c("CHROM.x","POS.x","ID.x","REF.x","ALT.x")]
                names(Total_annotations) = c("CHROM","POS","ID","REF","ALT")
                write.table(Total_annotations,file=paste0("/work/nyl112/data/synthetic_exome/Synthetic_exome_paralog_result_post_processed/Total_annotations_chrom_",j,"_",k,"_predicted_pathogenic.vcf"), na="", row.names = FALSE, col.names = TRUE, sep = ",", quote=FALSE)
                if (!is.null(Paraloc)){
                        save(Paraloc, file = paste0("/work/nyl112/data/synthetic_exome/Synthetic_exome_paralog_result_post_processed/Para_locations_chrom_",j,"_",k,".RData"))
                }
        }
}