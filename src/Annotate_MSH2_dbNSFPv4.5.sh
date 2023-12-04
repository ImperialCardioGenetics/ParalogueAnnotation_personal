#PBS -lwalltime=72:0:0
#PBS -lselect=1:ncpus=5:mem=10gb

module load anaconda3/personal
source /home/nyl112/.bashrc

cd /rds/general/project/lms-ware-analysis/live/nick/dbNSFP4.5 
java -Xmx5g search_dbNSFP45a -i /rds/general/project/lms-ware-analysis/live/nick/ParalogueAnnotation_personal/data/otherMaveDB_score/MSH2_mave_deleterious_b38.vcf -o /rds/general/project/lms-ware-analysis/live/nick/dbNSFP4.5/dbNSFP45a_MSH2_mave_deleterious_b38.vcf.out -p -w 1-6,38-40,59-61,77-84,96-101,108-110,122-140,147-149

java -Xmx5g search_dbNSFP45a -i /rds/general/project/lms-ware-analysis/live/nick/ParalogueAnnotation_personal/data/otherMaveDB_score/MSH2_mave_benign_b38.vcf -o /rds/general/project/lms-ware-analysis/live/nick/dbNSFP4.5/dbNSFP45a_MSH2_mave_benign_b38.vcf.out -p -w 1-6,38-40,59-61,77-84,96-101,108-110,122-140,147-149