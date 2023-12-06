#PBS -lwalltime=72:0:0
#PBS -lselect=1:ncpus=5:mem=10gb

module load anaconda3/personal
source /home/nyl112/.bashrc

cd /rds/general/project/lms-ware-analysis/live/nick/dbNSFP4.5 
java -Xmx5g search_dbNSFP45a -v hg19 -i /rds/general/project/lms-ware-analysis/live/nick/ParalogueAnnotation_personal/data/J_Li_data/PPARG_patho.vcf -o /rds/general/project/lms-ware-analysis/live/nick/dbNSFP4.5/dbNSFP45a_PPARG_patho.vcf.out -p -w 1-6,28,38-40,59-61,77-84,96-101,108-110,122-140,147-149

java -Xmx5g search_dbNSFP45a -v hg19 -i /rds/general/project/lms-ware-analysis/live/nick/ParalogueAnnotation_personal/data/J_Li_data/PPARG_benign.vcf -o /rds/general/project/lms-ware-analysis/live/nick/dbNSFP4.5/dbNSFP45a_PPARG_benign.vcf.out -p -w 1-6,28,38-40,59-61,77-84,96-101,108-110,122-140,147-149