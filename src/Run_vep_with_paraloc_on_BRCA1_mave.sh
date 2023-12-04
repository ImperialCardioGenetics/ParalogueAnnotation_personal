#PBS -lwalltime=72:0:0
#PBS -lselect=1:ncpus=1:mem=16gb

module load anaconda3/personal
source /home/nyl112/.bashrc

/home/nyl112/perl5/perlbrew/bin/perlbrew switch perl-5.22.0

perl -I /work/nyl112/paralogueAnnotator/ /work/nyl112/ensembl-vep/vep --force_overwrite --vcf --allele_number --canonical --offline --cache --dir_cache /work/nyl112/ -assembly GRCh38 -i /rds/general/project/lms-ware-analysis/live/nick/ParalogueAnnotation_personal/data/otherMaveDB_score/BRCA1_mave_benign_b38.vcf -o /rds/general/project/lms-ware-analysis/live/nick/ParalogueAnnotation_personal/data/otherMaveDB_score/BRCA1_mave_benign_b38.out_paraloc --plugin ParalogueAnnotation,paraloc

perl -I /work/nyl112/paralogueAnnotator/ /work/nyl112/ensembl-vep/vep --force_overwrite --vcf --allele_number --canonical --offline --cache --dir_cache /work/nyl112/ -assembly GRCh38 -i /rds/general/project/lms-ware-analysis/live/nick/ParalogueAnnotation_personal/data/otherMaveDB_score/BRCA1_mave_patho_b38.vcf -o /rds/general/project/lms-ware-analysis/live/nick/ParalogueAnnotation_personal/data/otherMaveDB_score/BRCA1_mave_patho_b38.out_paraloc --plugin ParalogueAnnotation,paraloc
