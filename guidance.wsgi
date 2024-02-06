#!/usr/bin/python3  

import sys
import site

#site.addsitedir('/var/www/flask/guidance/venv/lib/python3.8/site-packages') 
#site.addsitedir('/var/www/flask/guidance/venv/lib/python3.8/site-packages/flask')

sys.path.insert(0,"/var/www/flask/guidance")
sys.path.insert(0, '/var/www/flask/apptest/venv/lib/python3.10/site-packages/')

print()
from app import app as application
