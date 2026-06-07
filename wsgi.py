import os
import sys

APP_ROOT = "/var/www/vhosts/dev.guidance.tau.ac.il/httpdocs"
sys.path.insert(0, APP_ROOT)
sys.path.insert(1, "/guidance")  # keep this once the mount is back

# Force utils.py to treat the project root as the base
sys.argv[0] = os.path.join(APP_ROOT, "wsgi.py")
os.chdir(APP_ROOT)

from dotenv import load_dotenv
load_dotenv(os.path.join(APP_ROOT, ".env"))

from app import app as application

