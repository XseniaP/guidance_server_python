#!/usr/bin/python3  

import sys
import site
import os

path_to_project_folder = "/Users/kpolonsky/PycharmProjects/guidance_server_python/"
path_to_conda_env = "/Users/kpolonsky/miniconda3/envs/Guidance_py_env/"
path_to_conda_env_packages = "/Users/kpolonsky/miniconda3/envs/Guidance_py_env/lib/python3.9/site-packages"

#site.addsitedir('/var/www/flask/guidance/venv/lib/python3.8/site-packages')
#site.addsitedir('/var/www/flask/guidance/venv/lib/python3.8/site-packages/flask')

# sys.path.insert(0,"/var/www/flask/guidance")
# sys.path.insert(0, '/var/www/flask/apptest/venv/lib/python3.10/site-packages/')

sys.path.insert(0,path_to_project_folder)
sys.path.insert(0,path_to_conda_env_packages)

# # Set the Python executable path
os.environ['PYTHONHOME'] = path_to_conda_env
os.environ['PYTHONPATH'] = path_to_conda_env_packages

print()
from app import app as application
