
screen -R -D venv_screen
module avail 2>&1 | grep python
module load python/3.5.0
module load py_packages/3.5
module load tabix

# starting the virtual environment
virtualenv myvenv
source myvenv/bin/activate


# installing whatshap
module purge
module load python/3.5.0 py_packages/3.5 tabix
virtualenv venv
source venv/bin/activate
pip install --upgrade pip
pip install whatshap
