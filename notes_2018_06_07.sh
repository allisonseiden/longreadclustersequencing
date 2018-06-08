module purge
module load python/3.5.0 py_packages/3.5 tabix
virtualenv i_venv
source i_venv/bin/activate
pip install --upgrade pip
pip install whatshap
ssh interactive5
module load python/3.5.0 py_packages/3.5 tabix
source i_venv/bin/activate

# check if whatshap is loaded
whatshap 
