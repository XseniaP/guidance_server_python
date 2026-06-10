#!/powerapps/share/centos7/python-anaconda3.6.5/bin/python

#############################################################################################################
# this file should be saved as part of the pipeline and the cgi should import it rather than copy it twice! #
#############################################################################################################

import os, sys, stat
from utils import State
from enum import Enum
import platform

# constants to use when sending e-mails using the server admin's email address.
ADMIN_EMAIL = 'evolseq@tauex.tau.ac.il'
DEV_EMAIL = 'josefspr@gmail.com'
SMTP_SERVER = 'mxout.tau.ac.il'
ADMIN_USER_NAME = 'evolseq'
ADMIN_PASSWORD = 'elana'
SEND_EMAIL_DIR_IBIS = '/guidance/perl/'

OWNER_EMAIL = 'josefspr@gmail.com'

# general variables Ksenia
script_path = os.path.abspath(__file__)
Bin = f"{os.path.dirname(script_path)}/script" #path to the script folder
BIN_DIR = '/guidance/guidance_server_python' #os.path.dirname(Bin)  #main project folder
SERVERS_RESULTS_DIR = '/guidance/results' #os.path.join(BIN_DIR, 'results')
SERVERS_LOGS_DIR = '/guidance/results/logs' #os.path.join(BIN_DIR, 'logs')

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
WEBSERVER_URL = 'http://dev.guidance.tau.ac.il' #f'http://{WEBSERVER_NAME}.tau.ac.il'
WEBSERVER_URL_EXT = 'http://dev.guidance.tau.ac.il'
WEBSERVER_OLD_URL = f'http://{WEBSERVER_NAME}-old.tau.ac.il'
WEBSERVER_TITLE = '<b>Server for a multiple sequence alignment confidence score calculation</b>'

WEBSERVER_RESULTS_DIR = SERVERS_RESULTS_DIR #os.path.join(SERVERS_RESULTS_DIR, 'guidance')
WEBSERVER_LOGS_DIR = SERVERS_LOGS_DIR #os.path.join(SERVERS_RESULTS_DIR, 'guidance', 'logs')
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

SCRIPTS_DIR = os.path.join(BIN_DIR, 'script')
MAIN_SCRIPT = os.path.join(SCRIPTS_DIR, 'guidance_main.py')
EMAIL_SCRIPT = os.path.join(SCRIPTS_DIR, 'email_sender.py')
MASK_SCRIPT = os.path.join(SCRIPTS_DIR, 'mask_low_score_residues_webserver.py')
REMOVE_POS_SCRIPT = os.path.join(SCRIPTS_DIR, 'remove_pos_below_cutoff.py')
REMOVE_SEQ_SCRIPT = os.path.join(SCRIPTS_DIR, 'remove_seq_below_cutoff.py')
CONCAT_SCRIPT = os.path.join(SCRIPTS_DIR,'concat_aln_filelist_web.py')
HOT_PROGRAM = os.path.join(SCRIPTS_DIR, 'hot_cos_main.py')
MAFFT_OP_DIST = os.path.join(SCRIPTS_DIR, 'balibase.mafft_7123_mafft.op.Dist20bins.txt')
MAFFT_OP_DIST_0_25 = os.path.join(SCRIPTS_DIR, 'balibase.mafft_7123_mafft.op2.Dist25bins.txt')
MAFFT_EP_DIST_0_25 = os.path.join(SCRIPTS_DIR, 'balibase.mafft_7123_mafft.ep2.Dist20bins.txt')
HOT_GUIDANCE2_PROGRAM = os.path.join(SCRIPTS_DIR, 'hot_cos_main.py')
MIDPOINT_ROOTING_R = os.path.join(SCRIPTS_DIR, 'programs', 'MidPoint_Rooting.R')
MSA_Score_CSS = "http://dev.guidance.tau.ac.il/static/css/MSA_Colored.NEW.EM.css"
MidPoint_Rooting_R = os.path.join(SCRIPTS_DIR, 'programs', 'MidPoint_Rooting.R')

if platform.system()=='Darwin': # macOS platform
    MSA_SET_SCORE = os.path.join(SCRIPTS_DIR, 'programs', 'msa_set_score', 'msa_set_score')
    REMOVE_TAXA = os.path.join(SCRIPTS_DIR, 'programs','removeTaxa','removeTaxa')
    isEqualTopologyProg = os.path.join(SCRIPTS_DIR, 'programs', 'isEqualTree', 'isEqualTree')
    SEMPHY = os.path.join(BIN_DIR, 'script/programs/semphy/semphy')
    SEMPHY_BBL = os.path.join(BIN_DIR, 'script/programs/semphy/semphy')
    CLUSTAL_OMEGA = os.path.join(BIN_DIR, 'script/programs/clustalo')
    IQTREE = os.path.join(BIN_DIR, 'script/programs/iqtree/bin/iqtree2')

elif platform.system()=='Linux': # Linux platform
    MSA_SET_SCORE = os.path.join(SCRIPTS_DIR, 'programs', 'Linux','msa_set_score', 'msa_set_score')
    REMOVE_TAXA = os.path.join(SCRIPTS_DIR, 'programs', 'Linux', 'removeTaxa', 'removeTaxa')
    isEqualTopologyProg = os.path.join(SCRIPTS_DIR, 'programs', 'Linux', 'isEqualTree', 'isEqualTree')
    SEMPHY = os.path.join(BIN_DIR, 'script/programs/Linux/semphy/semphy')
    SEMPHY_BBL = os.path.join(BIN_DIR, 'script/programs/Linux/semphy/semphy')
    CLUSTAL_OMEGA = os.path.join(BIN_DIR, 'script/programs/Linux/clustalo')
    IQTREE = os.path.join(BIN_DIR, 'script/programs/Linux/iqtree/bin/iqtree2')

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
EXTERAL_SCRIPTS_PATH = '/home/josefspr/bioseq/bioSequence_scripts_and_constants'
Q_SUBMITTER_SCRIPT = f'{EXTERAL_SCRIPTS_PATH}/q_submitter_power_flask.py'
GUIDANCE_RUNNING_JOBS = f'{EXTERAL_SCRIPTS_PATH}/guidance_running_jobs.list'
JOB_QUEUE_NAME =  'pupkowebr@power9' #'itaym'


SUBMISSIONS_LOG =  os.path.join(WEBSERVER_LOGS_DIR, 'guidance.logv2')
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
WRITE_DAILY_TEST_SCRIPT = f"/guidance/script/write_daily_test_flask.py"
DAILY_TEST_DIR = "/guidance/daily_tests"

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
PATH2SAVE_PROCESS_DICT = r'/var/www/vhosts/dev.guidance.tau.ac.il/httpdocs/SavedObjects/processes.dict'
PATH2SAVE_WAITING_LIST = r'/var/www/vhosts/dev.guidance.tau.ac.il/httpdocs/SavedObjects/waiting.lst'
PATH2SAVE_PREVIOUS_DF = r'/var/www/vhosts/dev.guidance.tau.ac.il/httpdocs/SavedObjects/previous_processes.csv'

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
PATH2SAVE_MONITOR_DATA = '/guidance/results/Guidance' #f'{BIN_DIR}/results/Guidance'

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

