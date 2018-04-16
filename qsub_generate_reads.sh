
#!/bin/bash
#
#$ -l athena=true
#$ -j y                         # join std out & err
#$ -cwd                         # run job and store logs in cwd
#$ -N art_read_sim             # job name
#$ -l h_rt=10:00:00      # wall time request of 5 hour
#$ -pe smp 1             # parallel environment, 18 cores
#$ -l vf=5G              # memory needed 5G per core
#$ -l h_vmem=10G          # max mem request of 10 GB per core
#$ -V                    # import environment

art_illumina -sam -i $1 -p -l 150 -ss HS25 -f 20 -m 200 -s 10 -o ${1}.read_sim
