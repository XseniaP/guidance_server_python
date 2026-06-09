import os
import sys
import stat
from enum import Enum
import platform


class State(Enum):
    """Job lifecycle states — mirrors the definition in the root utils.py."""
    Running = 1
    Finished = 2
    Crashed = 3
    Waiting = 4
    Init = 5
    Queue = 6
    NotExists = 7
    Error = 8

    def str(self):
        return {
            State.Running: 'RUNNING',
            State.Crashed: 'FAILED',
            State.Error: 'ERROR',
            State.Finished: 'FINISHED',
            State.Waiting: 'WAITING',
            State.Init: 'INIT',
            State.Queue: 'IN QUEUE',
            State.NotExists: 'NOT EXISTS',
        }[self]

# constants to use when sending e-mails using the server admin's email address.
ADMIN_EMAIL = 'evolseq@gmail.com' #'TAU Evolseq <evolseq@tauex.tau.ac.il>' #TODO: Josef move the credentials to dot-env
DEV_EMAIL = 'josefspr@gmail.com' #TODO: Josef move the credentials to dot-env
SMTP_SERVER = 'smtp.gmail.com:587' #'mxout.tau.ac.il' #TODO: Josef move the credentials to dot-env
ADMIN_USER_NAME = 'evolseq' #TODO: Josef move the credentials to dot-env
ADMIN_PASSWORD = '' #TODO: Josef move the password to dot-env
SEND_EMAIL_DIR_IBIS = '/home/josefspr/bioseq/bioSequence_scripts_and_constants/sendEmail'

OWNER_EMAIL = 'josefspr@gmail.com' #TODO: Josef move the credentials to dot-env

# general variables Ksenia
# __file__ is guidance3/constants.py; BIN_DIR is the project root two levels up
_PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))  # guidance3/
BIN_DIR = os.path.dirname(_PACKAGE_DIR)                    # project root
Bin = os.path.join(BIN_DIR, "script")                      # script/ dir (programs live here)
SERVERS_RESULTS_DIR = os.path.join(BIN_DIR, 'results')
SERVERS_LOGS_DIR = os.path.join(BIN_DIR, 'logs')

RELOAD_INTERVAL = 30
RELOAD_TAGS = f'<META HTTP-EQUIV="REFRESH" CONTENT={RELOAD_INTERVAL}>'
NO_CACHE_TAGS = f'<META HTTP-EQUIV="PRAGMA" CONTENT="NO-CACHE">'
RUNNING_STATUS_TAG = '''<H1 align=center>Pepitope Job Status Page</h1>'''
FINISHED_STATUS_TAG = '''<H1 align=center>Pepitope Job Status Page - <font color='red'>FINISHED</font></h1>\n
<a href=#finish><H2 align=center>Go to the results</font></H2></a>\n'''
FAILED_STATUS_TAG = '''<H1 align=center>Pepitope Job Status Page - <font color='red'>FAILED</font></h1>\n
<a href=#finish><H2 align=center>Go to the results</font></H2></a>\n'''

# relevant modules

WEBSERVER_NAME_CAPITAL = 'Guidance'
WEBSERVER_NAME = 'guidance'
WEBSERVER_URL = '/guidance' #f'http://{WEBSERVER_NAME}.tau.ac.il'
WEBSERVER_URL_EXT = 'https://taux.evolseq.net/guidance'
WEBSERVER_OLD_URL = f'http://{WEBSERVER_NAME}-old.tau.ac.il'
WEBSERVER_TITLE = '<b>Server for a multiple sequence alignment confidence score calculation</b>'

WEBSERVER_RESULTS_DIR = os.path.join(SERVERS_RESULTS_DIR, 'Guidance')
WEBSERVER_LOGS_DIR = os.path.join(SERVERS_LOGS_DIR, 'Guidance')
WEBSERVER_HTML_DIR = f'/var/www/html/{WEBSERVER_NAME}/ver2'

WEBSERVER_RESULTS_URL = os.path.join(WEBSERVER_URL, 'results')
WEBSERVER_LOG_URL = os.path.join(WEBSERVER_URL, 'logs')
WEBSERVER_RESULTS_URL_EXT = os.path.join(WEBSERVER_URL_EXT, 'results')
WEBSERVER_PROCESS_STATE_URL = os.path.join(WEBSERVER_URL, 'process_state')
WEBSERVER_PROCESS_STATE_URL_EXT = os.path.join(WEBSERVER_URL_EXT, 'process_state')
WEBSERVER_RESULTS_OLD_URL = os.path.join(WEBSERVER_OLD_URL, 'results')
SOURCES = '/source.php'


#############################################################################################################
# Ksenia: START; SCRIPTs and Programs
#############################################################################################################

# Package-relative data and programs directories (work after pip install)
_DATA_DIR = os.path.join(_PACKAGE_DIR, 'data')
_PROGS_DIR = os.path.join(_PACKAGE_DIR, 'programs',
                           'mac' if platform.system() == 'Darwin' else 'linux')

# Script-like invocations — use the installed guidance3 package via -m so no script/ path needed
SCRIPTS_DIR = os.path.join(BIN_DIR, 'script')  # kept for any legacy callers
MASK_SCRIPT = f"{sys.executable} -m guidance3.sequences.filters mask"
REMOVE_POS_SCRIPT = f"{sys.executable} -m guidance3.sequences.filters remove_pos"
REMOVE_SEQ_SCRIPT = f"{sys.executable} -m guidance3.sequences.filters remove_seq"
CONCAT_SCRIPT = f"{sys.executable} -m guidance3.sequences.concat"
# MAIN_SCRIPT must remain a file path — deployment files (not ours) call it as `python3 {MAIN_SCRIPT}`
MAIN_SCRIPT = os.path.join(_PACKAGE_DIR, 'pipeline', 'main.py')
HOT_PROGRAM = f"{sys.executable} -m guidance3.hot_cos.main"
HOT_GUIDANCE3_PROGRAM = HOT_PROGRAM

# Data files — inside the guidance3 package
MAFFT_OP_DIST = os.path.join(_DATA_DIR, 'balibase.mafft_7123_mafft.op.Dist20bins.txt')
MAFFT_OP_DIST_0_25 = os.path.join(_DATA_DIR, 'balibase.mafft_7123_mafft.op2.Dist25bins.txt')
MAFFT_EP_DIST_0_25 = os.path.join(_DATA_DIR, 'balibase.mafft_7123_mafft.ep2.Dist20bins.txt')
MIDPOINT_ROOTING_R = os.path.join(_DATA_DIR, 'MidPoint_Rooting.R')
MidPoint_Rooting_R = MIDPOINT_ROOTING_R
FEATURES_EXTRACTION_MATRIX_DIR = os.path.join(_DATA_DIR, 'features_extraction', 'input_config_files')

MSA_Score_CSS = "/static/css/MSA_Colored.NEW.EM.css"

# Compiled binaries — inside the guidance3 package, platform-selected
MSA_SET_SCORE = os.path.join(_PROGS_DIR, 'msa_set_score')
REMOVE_TAXA = os.path.join(_PROGS_DIR, 'removeTaxa')
isEqualTopologyProg = os.path.join(_PROGS_DIR, 'isEqualTree')
SEMPHY = os.path.join(_PROGS_DIR, 'semphy')
SEMPHY_BBL = SEMPHY
CLUSTAL_OMEGA = os.path.join(_PROGS_DIR, 'clustalo')
IQTREE = os.path.join(_PROGS_DIR, 'iqtree', 'bin', 'iqtree2')
FEATURES_EXTRACTION_PROG = os.path.join(_PROGS_DIR, 'features_for_msas')

# DL model — inside the guidance3 package
DL_MODEL_PREDICT_SCRIPT = os.path.join(_PACKAGE_DIR, 'dl_model', 'scripts', 'predict_pretrained_main.py')
# Amino-acid model weights
DL_MODEL_PATH = os.path.join(_PACKAGE_DIR, 'dl_model', 'input', 'orthomam_model2', 'regressor_model_0_mode1_dseq_from_true.keras')
DL_MODEL_SCALER_PATH = os.path.join(_PACKAGE_DIR, 'dl_model', 'input', 'orthomam_model2', 'scaler_0_mode1_dseq_from_true.pkl')
# Nucleotide model weights
DL_MODEL_NUC_PATH = os.path.join(_PACKAGE_DIR, 'dl_model', 'input', 'nucleotides_model2', 'regressor_model_0_mode1_dseq_from_true.keras')
DL_MODEL_NUC_SCALER_PATH = os.path.join(_PACKAGE_DIR, 'dl_model', 'input', 'nucleotides_model2', 'scaler_0_mode1_dseq_from_true.pkl')

if platform.system() == 'Darwin':
    # Force arm64 slice — the webserver may run under Rosetta (x86_64), but numpy/TF are arm64-only.
    # Use sys.executable so venv/conda environments are respected.
    # DL_MODEL_PYTHON = ['arch', '-arm64', sys.executable]
    DL_MODEL_PYTHON = ['arch', '-arm64', '/Library/Frameworks/Python.framework/Versions/3.10/bin/python3']
else:
    DL_MODEL_PYTHON = os.environ.get('GUIDANCE_DL_PYTHON', sys.executable).split()

MUSCLE = "muscle"
MAFFT_GUIDANCE = "mafft"
PRANK_LECS = "prank"
PRANK = "prank"
PAGAN_LECS = "/share/apps/pagan-msa/bin/pagan"
PAGAN = 'pagan'
RUBY = 'ruby'

#############################################################################################################
# Ksenia: END;
#############################################################################################################

REQUIRED_MODULES = ['miniconda/miniconda3-4.7.12','python/python-3.6.7']

# MOVE AND UPDATE
EXTERAL_SCRIPTS_PATH = '/home/josefspr/bioseq/bioSequence_scripts_and_constants' #TODO: Josef to update if needed
Q_SUBMITTER_SCRIPT = f'{EXTERAL_SCRIPTS_PATH}/q_submitter_power_flask.py' #TODO: Josef to update if needed
GUIDANCE_RUNNING_JOBS = f'{EXTERAL_SCRIPTS_PATH}/guidance_running_jobs.list' #TODO: Josef to update if needed
JOB_QUEUE_NAME =  'pupkowebr@power9' #'itaym' #TODO: Josef to update the queue name


SUBMISSIONS_LOG = os.path.join(SERVERS_LOGS_DIR,WEBSERVER_NAME_CAPITAL, 'guidance.logv2')
RESULT_WEBPAGE_NAME = 'output.html'
EMAIL_FILE_NAME = 'email.txt'

MODE_0755 = stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH

CONTAINER_WIDTH = 'width: 850px'
CONTAINER_NO_MARGIN = 'margin: 0 auto'
CONTAINER_FONT = 'font-size: 20px'

CONTAINER_STYLE = f'{CONTAINER_WIDTH}; {CONTAINER_NO_MARGIN}; {CONTAINER_FONT}'

PROCESSING_MSG = f'<i>{WEBSERVER_NAME.upper()}</i> is now processing your request. This page will be automatically ' \
    f'updated every few seconds (until the job is done). You can also reload it manually. Once the job has finished, ' \
    f'several links to the output files will appear below. '

SYS_ERROR_MSG = "SYSTEM ERROR - GUIDANCE session has been terminated!"

PROGRESS_BAR_ANCHOR = '''<!--progress_bar_anchor-->'''
PROGRESS_BAR_TAG = '''<div class="progress"><div class="progress-bar progress-bar-striped active" role="progressbar" style="width:100%"></div></div>'''

MODE_0755 = stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH

DAILY_TEST_SEQUENCE = ">NC_001802.HXB2\nMQPIPIVAIVALVVAIIIAIVVWSIVIIEYRKILRQRKIDRLIDRLIERAEDSGNESEGEISALVEMGVEMGHHAPWDVDDL\n>EF637049.B\nMQSLQIVAIVALVVTAIIAIVVWSIVLIEYRKLLRQRKIDRLIDRIRERAEDSGNESEGDQEELAGLVERGHLAPWDVDDL\n>EF514700.B\nMQPLEILAIVALVVAIILAIVVWTIVFIEYKKILRQRKIDRLIDRIAERAEDSGNESEGDQEELSALVDMGHDAPWVVVDQ\n>DQ056417.C\nMLESIDYRLGVAALLLALIIAIIVWIIAYLEYRKLLRQRRIDKLIKRIRERAEDSGNESEGDIEELSTMVDVEHLRLLDVNNL\n>AY463217.C\nMVDLLAGVDYRVGVGALIIALIIAIIVWIWVYIEYRKLLRQRKIDWLIKRLREREEDSGNESEGDTEELATMVDMGHLRLLDDNNV\n>DQ011165.C\nMLNFLAGVDYRIGVGALIVGLIIAIVVWIIVYLEYRKLVKQRKIDWLIERIRERAEDSGNESEGDTEELATMVDMGHLRLLDAYDL\n>AB254142.C\nMINFAARVDYRVGVAAFTIALIIAIVVWIIVYLELVRQRKIDQLIIRIREREEDSGNESEGDIEELSTMVDMGQLRLLDGNGL\n>AY901969.C\nMVNLLEKVNLFEKVDYRLGVGALLIALVIAIIVWTIAYIEYRKLVRQRKIDWLVKRIRERAEDSGNESDGDTEELSTMVDLGHLRLLDVAEL\n>EU110088.A1\nMNQLQILAIXGLVVALILAIVVWTIVGIEYRKLLRQRRIDRLIKRISERAEDSGNESDGDTEELSQLVEMGNYNLGFDDNL\n>AB253428.A1\nMQLLEICAVVGLVVALIIAIVVWTIVGIEYKKLLKQRKIDRLVDRIRERAEDSGNESDGDREELSLLVDMGDYDLGDDNNL\n>AF457052.A1\nMLSALEICAIAGLVIALIIAIVVWTIVGIEYRRLLKQRKIDRLIERIRERAEDSGNESDGDTEELAALIEMGNYDLGDANDL\n>AF077336.F1\nMSYLLAIGIAALIVALIIAIVVWTIVYIEYKKLVRQRKINKLYKRIRERAEDSGNESEGDAEELAALGEMGPFIPGDINNL\n>DQ168575.G\nMKSLEISAIVGLIVAFIAAIVVWTIVLIEYRKIRKQKRIDKILDRIRERAEDSGNESEGDTEELATLVDMVDFEPWVGDNL\n>AY795907.D\nMQTLEILSIVALVIAAIIAIIVWTIVYIEYRKIRRQRKIDQLIDRIRERAEDSGNESEGDEEELSTLMEMGHAAPWNVADDL\n"
WRITE_DAILY_TEST_SCRIPT = f"{Bin}/write_daily_test_flask.py"
DAILY_TEST_DIR = "/home/josefspr/bioseq/bioSequence_scripts_and_constants/daily_tests/"

# from SharedConstants
class EMAIL_CONSTS:
    def create_title(state, job_name):
        if state == State.Finished:
            if job_name != "":
                return f'{WEBSERVER_NAME_CAPITAL} {job_name} - Job Finished'
            return f'{WEBSERVER_NAME_CAPITAL} - Job finished'
        elif state == State.Crashed:
            if job_name != "":
                return f'{WEBSERVER_NAME_CAPITAL} {job_name} - Job Crashed'
            return f'{WEBSERVER_NAME_CAPITAL} - Job Crashed'
        else:
            return f'unknown state in create_title at EMAIL_CONSTS'

    FINISHED_TITLE = f'{WEBSERVER_NAME_CAPITAL} - Job Finished'
    FINISHED_CONTENT = '''Thanks, for using Guidance\nYour results are at:\n{results_url}/{process_id}\nPlease, remember to cite us'''
    CRASHED_TITLE = f'{WEBSERVER_NAME_CAPITAL} - Job Failed'
    CRASHED_CONTENT =  '''Thanks, Your job has failed\nView your run at:\n{results_url}/{process_id}\n'''
    INIT_TITLE = f'{WEBSERVER_NAME_CAPITAL} - Your job has been submitted'
    INIT_CONTENT = '''Once the analysis will be ready, we will let you know! \nMeanwhile, you can track the progress of your job at:\n{results_url}/{process_id}'''

GUIDANCE_JOB_PREFIX = 'guidance'
MAIN_JOB_PREFIX = GUIDANCE_JOB_PREFIX
POSTPROCESS_JOB_PREFIX = 'PP'


# Job listener and management function naming
INTERVAL_BETWEEN_LISTENER_SAMPLES = 5  # in seconds
INTERVAL_BETWEEN_CLEANING_THE_PROCESSES_DICT = 24  # in hours
TIME_TO_SAVE_PROCESSES_IN_THE_PROCESSES_DICT = 7  # in days
LONG_RUNNING_JOBS_NAME = 'LongRunning'
QUEUE_JOBS_NAME = 'Queue'
NEW_RUNNING_JOBS_NAME = 'NewRunning'
FINISHED_JOBS_NAME = 'Finished'
ERROR_JOBS_NAME = 'Error'
WEIRD_BEHAVIOR_JOB_TO_CHECK = ''
PATH2SAVE_PROCESS_DICT = r'SavedObjects/processes.dict'
PATH2SAVE_WAITING_LIST = r'SavedObjects/waiting.lst'
PATH2SAVE_PREVIOUS_DF = r'SavedObjects/previous_processes.csv'

# PBS Listener consts
JOB_RUNNING_TIME_LIMIT_IN_HOURS = 10
JOB_NUMBER_COL = 'job_number'
JOB_NAME_COL = 'job_name'
JOB_STATUS_COL = 'job_status'
JOB_ELAPSED_TIME = 'elapsed_time'
JOB_CHANGE_COLS = [JOB_NUMBER_COL, JOB_NAME_COL, JOB_STATUS_COL]
QstatDataColumns = [JOB_NUMBER_COL, 'username', 'queue', JOB_NAME_COL, 'session_id', 'nodes', 'cpus', 'req_mem',
                    'req_time', JOB_STATUS_COL, JOB_ELAPSED_TIME]
SRVER_USERNAME = 'bioseq'

# Monitor consts
SEPERATOR_FOR_MONITOR_DF = '###'
#PATH2SAVE_MONITOR_DATA = r'SavedObjects/monitored_data'
# PATH2SAVE_MONITOR_DATA = r'/home/josefspr/results/Guidance'

# KSENIA
PATH2SAVE_MONITOR_DATA = f'{BIN_DIR}/results/Guidance'

class UI_CONSTS:

    states_text_dict = {
        State.Running: "Your process is running",
        State.Finished: "Your process finished... Redirecting to results page", #TODO is needed??
        State.Crashed: "Your process crashed\n we suggest you rerun the process.", #TODO finish
        State.Waiting: "We currently run other processes :( \n Your process will start soon",
        State.Init: "We are verifing your input, your process will start shortly",
        State.Queue: "Job is queued",
    }

    global allowed_files_str
    ALLOWED_EXTENSIONS = {'fasta', 'fastqc', 'gz'}
    allowed_files_str = ', '.join(ALLOWED_EXTENSIONS) #better to path string than list

    class UI_Errors(Enum):
        UNKNOWN_PROCESS_ID = 'The provided process id does not exist'
        INVALID_EXPORT_PARAMS ='invalid paramters for export'
        POSTPROCESS_CRASH = 'can\'t postprocess'
        INVALID_MAIL = 'invalid mail'
        CANT_ADD_PROCESS = 'can\'t add search process'
        INVALID_FILE = f'invalid file or file extenstion, please use a valid: {allowed_files_str} file'
        EXPORT_FILE_UNAVAILABLE = f'failed to export file, try to rerun the file'
        PAGE_NOT_FOUND = 'The requested page does not exist'

