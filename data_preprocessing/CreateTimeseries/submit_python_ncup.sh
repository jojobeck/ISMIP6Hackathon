#PBS -S /bin/bash                                                                                                                                                                                                                                                           
#PBS -P au88
#PBS -q normal
#PBS -l ncpus=48
#PBS -l walltime=12:00:00
#PBS -l mem=190GB
#PBS -l jobfs=200GB
#PBS -M johanna.beckmann@monash.edu
#PBS -l wd
#PBS -l software=matlab_monash
#PBS -o submit_python.outlog
#PBS -e submit_python.errlog
#PBS -l storage=gdata/au88

# load data
module purge
module load openmpi/4.0.2
module load python3/3.9.2
source /home/565/jb1863/datascience/bin/activate  # commented out by conda initialize
# source /home/565/jb1863/datascience/bin/activate
# module load openmpi/4.0.2
# module load intel-mkl/2020.3.304
# module load python3/3.9.2
# set number of open OMP threads
export OMP_NUM_THREADS=$PBS_NCPUS

# run python application
# python plot_2d_fiels_thick.py  > submit_python.log
python3 Calculate_and_save_data_dynSLR_anom.py > submit_python.log
