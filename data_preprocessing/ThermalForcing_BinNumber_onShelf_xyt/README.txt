'sampling_Thermal_foricng_files.py' creates a table and runs of ice sheet models for a certain exp/climate forcing.
Change i_e it you want to change the exp/climate forcing.
So far I have done all these only for 8km resolution only.
OUTPUT: - thermal forcing (no depth interpolation) on shelf (x,y,t)
OUTPUT: - Bin number on shelf (x,y,t)
all stored in
`/g/data/au88/ismip6_hackathon/dataset_2HD/modelname/exp_resolution/`
`submit_python.sh` is my bash script to submit this as a job on the cluster. 
I use a virtual enviroments in python as well, so the installed packeges of that are listet in `packages.txt`

have not coded anything for parallel python  :(
 
