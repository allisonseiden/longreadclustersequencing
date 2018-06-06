
screen -R -D venv_screen
module avail 2>&1 | grep python
module load python/3.5.0
module load py_packages/3.5
module load tabix

 virtualenv myvenv
