#! /bin/bash -login
#SBATCH -J CTCF_snakemake
#SBATCH -t 3-00:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -p high
#SBATCH --mem=2gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scpeng@ucdavis.edu


# -J job name
# -t time (minutes)
# -N number of nodes requested
# -n number of tasks to be run
# -c number of spus per task
# --mem minimum memory requested

# activate conda in general
conda activate base
cd /home/pengsc/projects/CTCF/

# Do something
mkdir -p logs_slurm
snakemake -j 12 --cluster-config /home/pengsc/projects/CTCF/cluster.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -N {cluster.nodes} -n {cluster.tasks}  -c {cluster.cpus} -J {cluster.name} -o {cluster.output} -e {cluster.output} --mail-type={cluster.email_type} --mail-user={cluster.email}" -s /home/pengsc/projects/CTCF/Snakefile -p --use-conda --jobs 12
