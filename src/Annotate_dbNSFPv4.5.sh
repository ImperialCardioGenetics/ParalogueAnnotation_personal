#PBS -lwalltime=72:0:0
#PBS -lselect=1:ncpus=1:mem=4gb

module load anaconda3/personal
source /home/nyl112/.bashrc

python /work/nyl112/Pfams/Pfam_meta_domains/src/Create_pfam_domain_mapping_objects.py /rds/general/project/lms-ware-analysis/live/nick/metadomains