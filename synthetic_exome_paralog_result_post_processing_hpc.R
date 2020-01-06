##RUN ON HPC - for trimming data down to only what's needed for paralog app website and also extracting predicted pathogenic variants for GeL##

Packages = c("tidyverse", "plyr", "dplyr")
new.packages = Packages[!(Packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "https://cran.ma.imperial.ac.uk/")
# lapply(Packages, library, character.only = TRUE)
# library("plyr")
library("dplyr")
# library("tidyverse")

# max_no_paralogous_pos = 0

##Chr 1 to 22##
for (j in 1:22){
        for (k in c("noQC","para_con","all_con")){
                Total_annotations = NULL
                files = list.files(path=paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RData_objects/chrom_",j,"/",k), pattern="*.RData", full.names=TRUE, recursive=FALSE)
                for (i in files){
                        print(i)
                        load(i)
                        tmp_Total_annotations = p.normal_PA$Total_paralog_annotations
                        tmp_Total_annotations = tmp_Total_annotations[,!(names(tmp_Total_annotations) %in% c("Variant_pos","FILTER.x","BIOTYPE.x","Paralogue_Vars.x","REF_Amino_acids.x","ALT_Amino_acids.x","paralog","paralog_pos","FILTER.y","BIOTYPE.y","Paralogue_Vars.y","REF_Amino_acids.y","ALT_Amino_acids.y"))]
                        tmp_Total_annotations = tmp_Total_annotations[!is.na(tmp_Total_annotations$POS.x),]
                        if (is.null(Total_annotations)){
                                Total_annotations = tmp_Total_annotations
                        } else {
                                Total_annotations = base::rbind(Total_annotations, dplyr::setdiff(tmp_Total_annotations, Total_annotations))
                        }
                }
                save(Total_annotations, file = paste0("/work/nyl112/data/synthetic_exome/Total_annotations_chrom_",j,"_",k,".RData"))
                Total_annotations = Total_annotations[,c("CHROM.x","POS.x","ID.x","REF.x","ALT.x")]
                names(Total_annotations) = c("CHROM","POS","ID","REF","ALT")
                write.table(Total_annotations,file=paste0("/work/nyl112/data/synthetic_exome/Total_annotations_chrom_X_",k,"_predicted_pathogenic.vcf"), na="", row.names = FALSE, col.names = FALSE, sep = ",", quote=FALSE)
        }
}

##ChrX##
for (k in c("noQC","para_con","all_con")){
        Total_annotations = NULL
        files = list.files(path=paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RData_objects/chrom_X/",k), pattern="*.RData", full.names=TRUE, recursive=FALSE)
        # print(files)
        for (i in files){
                print(i)
                load(i)
                # All_possible_query_genes = c(All_possible_query_genes, unique(p.normal_PA$Left_joined_gathered_paralog_data$Gene))
                if (is.null(Total_annotations)){
                        Total_annotations = p.normal_PA$Unique_variant_gene_annotations
                } else {
                        Total_annotations = base::rbind(Total_annotations, dplyr::setdiff(p.normal_PA$Unique_variant_gene_annotations, Total_annotations))
                }                
                # if (p.normal_PA$max_no_col > max_no_paralogous_pos){
                #         max_no_paralogous_pos = p.normal_PA$max_no_col
                # }
        }
        Genes = unique(Total_annotations$Gene)
        chr_pos = as.numeric(sapply(strsplit(as.character(Total_annotations$Variant_pos), split = " "), "[", 2))
        # print(Total_annotations)
        var_positions_data = data.frame(Chrom = "X", Position = chr_pos)
        write.table(Genes,file=paste0("/work/nyl112/data/synthetic_exome/synthetic_exome_chrom_X_",k,"_genes.txt"), na="", row.names = FALSE, col.names = FALSE, sep = ",", quote=FALSE)
        save(var_positions_data, file = paste0("/work/nyl112/data/synthetic_exome/var_positions_data_chrom_X_",k,".RData"))
}

##ChrY##
for (k in c("noQC","para_con","all_con")){
        Total_annotations = NULL
        files = list.files(path=paste0("/work/nyl112/data/synthetic_exome/paralogous_var_align.RData_objects/chrom_Y/",k), pattern="*.RData", full.names=TRUE, recursive=FALSE)
        for (i in files){
                print(i)
                load(i)
                # All_possible_query_genes = c(All_possible_query_genes, unique(p.normal_PA$Left_joined_gathered_paralog_data$Gene))
                if (is.null(Total_annotations)){
                        Total_annotations = p.normal_PA$Unique_variant_gene_annotations
                } else {
                        Total_annotations = base::rbind(Total_annotations, dplyr::setdiff(p.normal_PA$Unique_variant_gene_annotations, Total_annotations))
                }               
                # if (p.normal_PA$max_no_col > max_no_paralogous_pos){
                #         max_no_paralogous_pos = p.normal_PA$max_no_col
                # }
        }
        Genes = unique(Total_annotations$Gene)
        chr_pos = as.numeric(sapply(strsplit(as.character(Total_annotations$Variant_pos), split = " "), "[", 2))
        var_positions_data = data.frame(Chrom = "Y", Position = chr_pos)
        write.table(Genes,file=paste0("/work/nyl112/data/synthetic_exome/synthetic_exome_chrom_Y_",k,"_genes.txt"), na="", row.names = FALSE, col.names = FALSE, sep = ",", quote=FALSE)
        save(var_positions_data, file = paste0("/work/nyl112/data/synthetic_exome/var_positions_data_chrom_Y_",k,".RData"))
}