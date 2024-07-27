# import os, sys, stat
# from utils import State
# from enum import Enum
#
# # constants to use when sending e-mails using the server admin's email address.
# ADMIN_EMAIL = 'evolseq@gmail.com'  # 'TAU Evolseq <evolseq@tauex.tau.ac.il>'
# DEV_EMAIL = 'josefspr@gmail.com'
# SMTP_SERVER = 'smtp.gmail.com:587'  # 'mxout.tau.ac.il'
# ADMIN_USER_NAME = 'evolseq'
# ADMIN_PASSWORD = 'yipnqomnsofhytqp'
# SEND_EMAIL_DIR_IBIS = '/home/josefspr/bioseq/bioSequence_scripts_and_constants/sendEmail'
#
# OWNER_EMAIL = 'josefspr@gmail.com'
#
# # general variables
# script_path = os.path.abspath(__file__)
# Bin = os.path.dirname(script_path) #path to the script folder
# # Bin = os.path.dirname(sys.argv[0])
# BIN_DIR = os.path.dirname(Bin)    #path to the main project folder
# SERVERS_RESULTS_DIR = os.path.join(BIN_DIR, 'results')
# SERVERS_LOGS_DIR = os.path.join(BIN_DIR, 'logs')
#
# RELOAD_INTERVAL = 30
# RELOAD_TAGS = f'<META HTTP-EQUIV="REFRESH" CONTENT={RELOAD_INTERVAL}>'
# NO_CACHE_TAGS = f'<META HTTP-EQUIV="PRAGMA" CONTENT="NO-CACHE">'
# RUNNING_STATUS_TAG = '''<H1 align=center>Pepitope Job Status Page</h1>'''
# FINISHED_STATUS_TAG = '''<H1 align=center>Pepitope Job Status Page - <font color='red'>FINISHED</font></h1>\n
# <a href=#finish><H2 align=center>Go to the results</font></H2></a>\n'''
# FAILED_STATUS_TAG = '''<H1 align=center>Pepitope Job Status Page - <font color='red'>FAILED</font></h1>\n
# <a href=#finish><H2 align=center>Go to the results</font></H2></a>\n'''
#
# # relevant modules
#
# WEBSERVER_NAME_CAPITAL = 'Guidance'
# WEBSERVER_NAME = 'guidance'
# WEBSERVER_URL = '/guidance'  # f'http://{WEBSERVER_NAME}.tau.ac.il'
# WEBSERVER_URL_EXT = 'https://taux.evolseq.net/guidance'
# WEBSERVER_OLD_URL = f'http://{WEBSERVER_NAME}-old.tau.ac.il'
# WEBSERVER_TITLE = '<b>Server for epitope mapping using affinity-selected peptides</b>'
#
# WEBSERVER_RESULTS_DIR = os.path.join(SERVERS_RESULTS_DIR, 'Guidance')
# WEBSERVER_LOGS_DIR = os.path.join(SERVERS_LOGS_DIR, 'Guidance')
# WEBSERVER_HTML_DIR = f'/var/www/html/{WEBSERVER_NAME}/ver2'
#
# WEBSERVER_RESULTS_URL = os.path.join(WEBSERVER_URL, 'results')
# WEBSERVER_LOG_URL = os.path.join(WEBSERVER_URL, 'logs')
# WEBSERVER_RESULTS_URL_EXT = os.path.join(WEBSERVER_URL_EXT, 'results')
# WEBSERVER_PROCESS_STATE_URL = os.path.join(WEBSERVER_URL, 'process_state')
# WEBSERVER_PROCESS_STATE_URL_EXT = os.path.join(WEBSERVER_URL_EXT, 'process_state')
# WEBSERVER_RESULTS_OLD_URL = os.path.join(WEBSERVER_OLD_URL, 'results')
# SOURCES = '/source.php'
#
# SCRIPTS_DIR = os.path.join(BIN_DIR, 'script')
#
# MAIN_SCRIPT = os.path.join(SCRIPTS_DIR, 'guidance_main.py')
#
# MASK_SCRIPT = os.path.join(SCRIPTS_DIR, 'mask_low_score_residues_webserver.py')
#
# REMOVE_POS_SCRIPT = os.path.join(SCRIPTS_DIR, 'remove_pos_below_cutoff.py')
#
# REMOVE_SEQ_SCRIPT = os.path.join(SCRIPTS_DIR, 'remove_seq_below_cutoff.py')
#
# CONCAT_SCRIPT = os.path.join(SCRIPTS_DIR, 'concat_aln_filelist_web.py')
#
# MSA_SET_SCORE = os.path.join(SCRIPTS_DIR, 'programs', 'msa_set_score', 'msa_set_score')
# HOT_PROGRAM = os.path.join(SCRIPTS_DIR, 'hot_cos_main.py')
# MAFFT_OP_DIST = os.path.join(SCRIPTS_DIR, 'balibase.mafft_7123_mafft.op.Dist20bins.txt')
# MAFFT_OP_DIST_0_25 = os.path.join(SCRIPTS_DIR, 'balibase.mafft_7123_mafft.op2.Dist25bins.txt')
# MAFFT_EP_DIST_0_25 = os.path.join(SCRIPTS_DIR, 'balibase.mafft_7123_mafft.ep2.Dist20bins.txt')
# HOT_GUIDANCE2_PROGRAM = os.path.join(SCRIPTS_DIR, 'hot_cos_main.py')
# MIDPOINT_ROOTING_R = os.path.join(SCRIPTS_DIR, 'programs', 'MidPoint_Rooting.R')
# MSA_Score_CSS = "https://taux.evolseq.net/guidance/static/css/MSA_Colored.NEW.EM.css"
# MidPoint_Rooting_R = os.path.join(SCRIPTS_DIR, 'programs', 'MidPoint_Rooting.R')
# isEqualTopologyProg = os.path.join(SCRIPTS_DIR, 'programs', 'isEqualTree', 'isEqualTree')
#
#
#
# REQUIRED_MODULES = ['miniconda/miniconda3-4.7.12', 'python/python-3.6.7']
#
# # move to project folder
# EXTERAL_SCRIPTS_PATH = '/Users/kpolonsky/Documents/GUIDANCE-guidance.v2.02/www/bioSequence_scripts_and_constants'
# # EXTERAL_SCRIPTS_PATH = '/home/josefspr/bioseq/bioSequence_scripts_and_constants'
# Q_SUBMITTER_SCRIPT = f'{EXTERAL_SCRIPTS_PATH}/q_submitter_power_flask.py'
# GUIDANCE_RUNNING_JOBS = f'{EXTERAL_SCRIPTS_PATH}/guidance_running_jobs.list'
# JOB_QUEUE_NAME = 'pupkowebr@power9'  # 'itaym'
#
# SUBMISSIONS_LOG = os.path.join(SERVERS_LOGS_DIR, WEBSERVER_NAME_CAPITAL, 'guidance.logv2')
# RESULT_WEBPAGE_NAME = 'output.html'
# EMAIL_FILE_NAME = 'email.txt'
#
# MODE_0755 = stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH
#
# CONTAINER_WIDTH = 'width: 850px'
# CONTAINER_NO_MARGIN = 'margin: 0 auto'
# CONTAINER_FONT = 'font-size: 20px'
#
# CONTAINER_STYLE = f'{CONTAINER_WIDTH}; {CONTAINER_NO_MARGIN}; {CONTAINER_FONT}'
#
# PROCESSING_MSG = f'<i>{WEBSERVER_NAME.upper()}</i> is now processing your request. This page will be automatically ' \
#                  f'updated every few seconds (until the job is done). You can also reload it manually. Once the job has finished, ' \
#                  f'several links to the output files will appear below. '
#
# SYS_ERROR_MSG = "SYSTEM ERROR - GUIDANCE session has been terminated!"
#
# PROGRESS_BAR_ANCHOR = '''<!--progress_bar_anchor-->'''
# PROGRESS_BAR_TAG = '''<div class="progress"><div class="progress-bar progress-bar-striped active" role="progressbar" style="width:100%"></div></div>'''
#
# MODE_0755 = stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH
#
# DAILY_TEST_SEQUENCE = ">NC_001802.HXB2\nMQPIPIVAIVALVVAIIIAIVVWSIVIIEYRKILRQRKIDRLIDRLIERAEDSGNESEGEISALVEMGVEMGHHAPWDVDDL\n>EF637049.B\nMQSLQIVAIVALVVTAIIAIVVWSIVLIEYRKLLRQRKIDRLIDRIRERAEDSGNESEGDQEELAGLVERGHLAPWDVDDL\n>EF514700.B\nMQPLEILAIVALVVAIILAIVVWTIVFIEYKKILRQRKIDRLIDRIAERAEDSGNESEGDQEELSALVDMGHDAPWVVVDQ\n>DQ056417.C\nMLESIDYRLGVAALLLALIIAIIVWIIAYLEYRKLLRQRRIDKLIKRIRERAEDSGNESEGDIEELSTMVDVEHLRLLDVNNL\n>AY463217.C\nMVDLLAGVDYRVGVGALIIALIIAIIVWIWVYIEYRKLLRQRKIDWLIKRLREREEDSGNESEGDTEELATMVDMGHLRLLDDNNV\n>DQ011165.C\nMLNFLAGVDYRIGVGALIVGLIIAIVVWIIVYLEYRKLVKQRKIDWLIERIRERAEDSGNESEGDTEELATMVDMGHLRLLDAYDL\n>AB254142.C\nMINFAARVDYRVGVAAFTIALIIAIVVWIIVYLELVRQRKIDQLIIRIREREEDSGNESEGDIEELSTMVDMGQLRLLDGNGL\n>AY901969.C\nMVNLLEKVNLFEKVDYRLGVGALLIALVIAIIVWTIAYIEYRKLVRQRKIDWLVKRIRERAEDSGNESDGDTEELSTMVDLGHLRLLDVAEL\n>EU110088.A1\nMNQLQILAIXGLVVALILAIVVWTIVGIEYRKLLRQRRIDRLIKRISERAEDSGNESDGDTEELSQLVEMGNYNLGFDDNL\n>AB253428.A1\nMQLLEICAVVGLVVALIIAIVVWTIVGIEYKKLLKQRKIDRLVDRIRERAEDSGNESDGDREELSLLVDMGDYDLGDDNNL\n>AF457052.A1\nMLSALEICAIAGLVIALIIAIVVWTIVGIEYRRLLKQRKIDRLIERIRERAEDSGNESDGDTEELAALIEMGNYDLGDANDL\n>AF077336.F1\nMSYLLAIGIAALIVALIIAIVVWTIVYIEYKKLVRQRKINKLYKRIRERAEDSGNESEGDAEELAALGEMGPFIPGDINNL\n>DQ168575.G\nMKSLEISAIVGLIVAFIAAIVVWTIVLIEYRKIRKQKRIDKILDRIRERAEDSGNESEGDTEELATLVDMVDFEPWVGDNL\n>AY795907.D\nMQTLEILSIVALVIAAIIAIIVWTIVYIEYRKIRRQRKIDQLIDRIRERAEDSGNESEGDEEELSTLMEMGHAAPWNVADDL\n"
# WRITE_DAILY_TEST_SCRIPT = "/home/josefspr/bioseq/guidance/guidance.v2.02/www/Guidance/write_daily_test_flask.py"
# DAILY_TEST_DIR = "/home/josefspr/bioseq/bioSequence_scripts_and_constants/daily_tests/"
#
#
# # from SharedConstants
# class EMAIL_CONSTS:
#     def create_title(state, job_name):
#         if state == State.Finished:
#             if job_name != "":
#                 return f'{WEBSERVER_NAME_CAPITAL} {job_name} - Job Finished'
#             return f'{WEBSERVER_NAME_CAPITAL} - Job finished'
#         elif state == State.Crashed:
#             if job_name != "":
#                 return f'{WEBSERVER_NAME_CAPITAL} {job_name} - Job Crashed'
#             return f'{WEBSERVER_NAME_CAPITAL} - Job Crashed'
#         else:
#             return f'unknown state in create_title at EMAIL_CONSTS'
#
#     FINISHED_TITLE = f'{WEBSERVER_NAME_CAPITAL} - Job Finished'
#     FINISHED_CONTENT = '''Thanks, for using Guidance\nYour results are at:\n{results_url}/{process_id}\nPlease, remember to cite us'''
#     CRASHED_TITLE = f'{WEBSERVER_NAME_CAPITAL} - Job Failed'
#     CRASHED_CONTENT = '''Thanks, Your job has failed\nView your run at:\n{results_url}/{process_id}\n'''
#     INIT_TITLE = f'{WEBSERVER_NAME_CAPITAL} - Your job has been submitted'
#     INIT_CONTENT = '''Once the analysis will be ready, we will let you know! \nMeanwhile, you can track the progress of your job at:\n{results_url}/{process_id}'''
#
#
# GUIDANCE_JOB_PREFIX = 'guidance'
# MAIN_JOB_PREFIX = GUIDANCE_JOB_PREFIX
# POSTPROCESS_JOB_PREFIX = 'PP'
#
# # Job listener and management function naming
# INTERVAL_BETWEEN_LISTENER_SAMPLES = 5  # in seconds
# INTERVAL_BETWEEN_CLEANING_THE_PROCESSES_DICT = 24  # in hours
# TIME_TO_SAVE_PROCESSES_IN_THE_PROCESSES_DICT = 7  # in days
# LONG_RUNNING_JOBS_NAME = 'LongRunning'
# QUEUE_JOBS_NAME = 'Queue'
# NEW_RUNNING_JOBS_NAME = 'NewRunning'
# FINISHED_JOBS_NAME = 'Finished'
# ERROR_JOBS_NAME = 'Error'
# WEIRD_BEHAVIOR_JOB_TO_CHECK = ''
# PATH2SAVE_PROCESS_DICT = r'SavedObjects/processes.dict'
# PATH2SAVE_WAITING_LIST = r'SavedObjects/waiting.lst'
# PATH2SAVE_PREVIOUS_DF = r'SavedObjects/previous_processes.csv'
#
# # PBS Listener consts
# JOB_RUNNING_TIME_LIMIT_IN_HOURS = 10
# JOB_NUMBER_COL = 'job_number'
# JOB_NAME_COL = 'job_name'
# JOB_STATUS_COL = 'job_status'
# JOB_ELAPSED_TIME = 'elapsed_time'
# JOB_CHANGE_COLS = [JOB_NUMBER_COL, JOB_NAME_COL, JOB_STATUS_COL]
# QstatDataColumns = [JOB_NUMBER_COL, 'username', 'queue', JOB_NAME_COL, 'session_id', 'nodes', 'cpus', 'req_mem',
#                     'req_time', JOB_STATUS_COL, JOB_ELAPSED_TIME]
# SRVER_USERNAME = 'bioseq'
#
# # Monitor consts
# SEPERATOR_FOR_MONITOR_DF = '###'
# # PATH2SAVE_MONITOR_DATA = r'SavedObjects/monitored_data'
# # PATH2SAVE_MONITOR_DATA = r'/home/josefspr/results/Guidance'
# PATH2SAVE_MONITOR_DATA = r'/Users/kpolonsky/PycharmProjects/guidance_server/results/Guidance'
#
#
# class UI_CONSTS:
#     states_text_dict = {
#         State.Running: "Your process is running",
#         State.Finished: "Your process finished... Redirecting to results page",  # TODO is needed??
#         State.Crashed: "Your process crashed\n we suggest you rerun the process.",  # TODO finish
#         State.Waiting: "We currently run other processes :( \n Your process will start soon",
#         State.Init: "We are verifing your input, your process will start shortly",
#         State.Queue: "Job is queued",
#     }
#
#     global allowed_files_str
#     ALLOWED_EXTENSIONS = {'fasta', 'fastqc', 'gz'}
#     allowed_files_str = ', '.join(ALLOWED_EXTENSIONS)  # better to path string than list
#
#     class UI_Errors(Enum):
#         UNKNOWN_PROCESS_ID = 'The provided process id does not exist'
#         INVALID_EXPORT_PARAMS = 'invalid paramters for export'
#         POSTPROCESS_CRASH = 'can\'t postprocess'
#         INVALID_MAIL = 'invalid mail'
#         CANT_ADD_PROCESS = 'can\'t add search process'
#         INVALID_FILE = f'invalid file or file extenstion, please use a valid: {allowed_files_str} file'
#         EXPORT_FILE_UNAVAILABLE = f'failed to export file, try to rerun the file'
#         PAGE_NOT_FOUND = 'The requested page does not exist'
#
# #
# # # constants to use when sending e-mails using the server admin's email address.
# # ADMIN_EMAIL = "TAU BioSequence <bioSequence@tauex.tau.ac.il>"
# # ADMIN_USER_NAME = "bioSequence"
# # ADMIN_PASSWORD = ""
# # SMTP_SERVER = ""
# #
# # # the name of the list of all running processes
# # QUEUING_JOBS = "/bioseq/bioSequence_scripts_and_constants/queuing_jobs.list"
# # RUNNING_JOBS = "/bioseq/bioSequence_scripts_and_constants/running_jobs.list"
# # SUBMITTED_JOBS = "/bioseq/bioSequence_scripts_and_constants/submitted_jobs.list"
# # JOBS_ON_BIOSEQ_NODE = "/bioseq/bioSequence_scripts_and_constants/jobs_on_bioc.01_node.list"
# # JOBS_WAITING_BIOSEQ_NODE = "/bioseq/bioSequence_scripts_and_constants/jobs_waiting_bioc.01_node.list"
# # CONSURF_RUNNING_JOBS = "/bioseq/bioSequence_scripts_and_constants/consurf_running_jobs.list"
# # SELECTON_RUNNING_JOBS = "/bioseq/bioSequence_scripts_and_constants/selecton_running_jobs.list"
# # CONSEQ_RUNNING_JOBS = "/bioseq/bioSequence_scripts_and_constants/conseq_running_jobs.list"
# # PEPITOPE_RUNNING_JOBS = "/bioseq/bioSequence_scripts_and_constants/pepitope_running_jobs.list"
# #
# # # Databases urls
# # PROTEOPEDIA = "http://proteopedia.org/wiki/index.php/"
# # PDB_DB = "http://www.rcsb.org/pdb/explore/explore.do?structureId="
# # RCSB_WGET = "wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/"
# # RCSB = "http://www.rcsb.org/"
# # PISA_WGET = "wget http://www.ebi.ac.uk/msd-srv/pisa/cgi-bin/multimer.pdb?"
# #
# # # CGIs paths
# # CONSURF_CGI_DIR = "/var/www/cgi-bin/ConSurf"
# #
# # #general paths
# # SERVERS_RESULTS_DIR = "/bioseq/data/results/"
# # SERVERS_LOGS_DIR = "/bioseq/data/logs/"
# # SEND_EMAIL_DIR = "/bioseq/bioSequence_scripts_and_constants/sendEmail"
# # SEND_EMAIL_DIR_IBIS = "/bioseq/bioSequence_scripts_and_constants/sendEmail"
# # DAEMON_LOG_FILE = "/bioseq/bioSequence_scripts_and_constants/daemon.log"
# # UPDATE_RUN_TIME_LOG_FILE = "/bioseq/bioSequence_scripts_and_constants/update_runTime.log"
# # CONSURF_CGI = "/var/www/cgi-bin/ConSurf"
# # BIOSEQ_TEMP = "/bioseq/temp/"
# #
# # # servers urls:
# # SELECTON_URL = "http://selecton.tau.ac.il"
# # CONSEQ_URL = "http://conseq.tau.ac.il/"
# # CONSURF_URL = "http://consurf.tau.ac.il/"
# # NEW_CONSURF_URL = "http://consurf.tau.ac.il/"
# # EPITOPIA_URL = "http://epitopia.tau.ac.il/"
# # PEPITOPE_URL = "http://pepitope.tau.ac.il/"
# # QMF_URL = "http://quasimotifinder.tau.ac.il/"
# # PATCHFINDER_URL = "http://patchfinder.tau.ac.il/"
# # FASTML_URL = "http://fastml.tau.ac.il/"
# # RECONST_URL = "http://fastml.tau.ac.il/reconst/"
# # GAIN_LOSS_URL = "http://gloome.tau.ac.il/"
# # CONSURF_DB_URL = "http://consurfdb.tau.ac.il/"
# # GILAD_SERVER_URL = "http://mud.tau.ac.il/"
# # MCPep_URL = "http://bental.tau.ac.il/MCPep/"
# # GUIDANCE_URL = "http://guidance.tau.ac.il/"
# # GUIDANCE_INDELS_URL = "http://guidance.tau.ac.il/indels/"
# # SPECBOOST_URL = "http://bental.tau.ac.il/specBoost/"
# # PROMAYA_URL = "http://bental.tau.ac.il/ProMaya/"
# # HOMOLOGY_SEARCH_URL = "http://fastml.tau.ac.il/HomologySearch/"
# # COPAP_URL = "http://copap.tau.ac.il/"
# #
# # #servers logs:
# # CONSURF_LOG = "/bioseq/ConSurf_old/consurf.log"
# # CONSURF_NEW_LOG = "/bioseq/ConSurf/consurf.log"
# # SELECTON_LOG = "/bioseq/Selecton/selecton.log"
# # EPITOPIA_LOG = "/bioseq/epitopia/epitopia.log"
# # CONSEQ_LOG = "/bioseq/ConSeq/conseq.log"
# # PEPITOPE_LOG = "/bioseq/pepitope/pepitope.log"
# # RECONST_LOG = "/bioseq/ReConst_Server/reconst.log"
# # MCPep_LOG = "/bioseq/MCPep/mcpep.log"
# # Guidance_LOG = "/bioseq/Guidance/guidance.log"
# # Guidance_Indels_LOG = "/bioseq/GuidanceIndels/guidance_Indels.log"
# # MuD_LOG = "/bioseq/Gilad_Server/MuD.log"
# # FASTML_LOG = "/bioseq/FastML/fastml.log"
# # SPECBOOST_LOG = "/bioseq/specBoost/specBoost.log"
# # GAIN_LOSS_LOG = "/bioseq/GainLoss/GainLoss.log"
# # PROMAYA_LOG = "/bioseq/ProMaya/ProMaya.log"
# # COPAP_LOG = "/bioseq/CoPAP/CoPAP.log"
# #
# # #servers results urls:
# # # servers urls:
# # SELECTON_RESULTS_URL = SELECTON_URL + "/results/"
# #
# # #external databases
# # PQS = "/biodb/PQS/"
# # PDB_DIVIDED = "/biodb/PDB/data/structures/divided/pdb/"
# # SWISSPROT_DB = "/biodb/BLAST/Proteins/swissprot"
# # UNIPROT_DB = "/biodb/BLAST/Proteins/uniprot"
# # CLEAN_UNIPROT_DB = "/biodb/BLAST/Proteins/clean_uniprot"
# # UNIREF90_DB = "/biodb/BLAST/Proteins/uniref90"
# # PDBAA_NCBI = "/biodb/BLAST/Proteins/pdbaa"
# # CULLED_PDB = "/groups/bioseq.home/HAIM/PDBAA/pdbaaent"
# # PDB_DUNBRACK = "/groups/bioseq.home/HAIM/PDBAA/pdbaa"
# # NR_PROT_DB = "/biodb/BLAST/Proteins/nr"
# # NR_NUC_DB = "/biodb/BLAST/Nucleotides/nt"
# # UNIPROT_DAT_INDEX = "/bioseq/data/results/GB_CDS/uniprot.dat.bp_index"
# # PDB_TO_UNIPROT = "/bioseq/data/results/PDB_to_UNIPROT/idmapping_PDB_UNIPROTKB.dat"
# # PDB_TO_UNIPROT_test = "/biodb/idmapping_PDB_UNIPROTKB.dat"
# #
# #
# # #internal databases
# # # EPITOPIA_DATA = "/bioseq/epitopia/data"
# #
# # #external programs
# #
# # # Paths to different BLAST programs
# # BLASTALL = "/opt/bio/ncbi/bin/blastall"
# # BLASTPGP = "blastpgp"
# # CS_BLAST = "/share/apps/csblast-2.1.0-linux64/csblast_static"
# #
# # # Paths to MUSCLE program
# # MUSCLE_LECS = "/share/apps/bin/muscle"
# # # MUSCLE = "/usr/local/bin/muscle"
# # # MUSCLE_3_6 = "/bioseq/Programs/muscle_3.6_from_BIOCLUSTER/muscle3.6/muscle"
# # # MUSCLE_LECS = "muscle"
# # MUSCLE = "muscle"
# #
# # # Paths to different versions of ClustalW program
# # CLUSTALW_LECS = "/share/apps/bin/clustalw"
# # CLUSTALW = "/usr/local/bin/clustalw"
# # CLUSTALW_1_82 = "/bioseq/Programs/ClustalW_1.82/clustalw1.82/clustalw"
# # CLUSTALW_1_81 = "/bioseq/Programs/ClustalW_1.81/clustalw1.81/clustalw"
# # CLUSTALW_2_0_10 = "/bioseq/Programs/ClustalW_2.0.10/clustalw-2.0.10-linux-i386-libcppstatic/clustalw2"
# # CLUSTAL_OMEGA = os.path.join(BIN_DIR, 'script/programs/clustalo')
# #
# # # Paths to different versions of MAFFT program
# # MAFFT_LINSI = "/usr/local/bin/mafft-linsi"
# # MAFFT = "/usr/local/bin/mafft"
# # # MAFFT_GUIDANCE = "/bioseq/Programs/MAFFT_6.833/bin/mafft"
# # MAFFT_GUIDANCE = "mafft"
# # MAFFT_LINSI_GUIDANCE = "/bioseq/Programs/MAFFT_6.833/bin/mafft --localpair --maxiterate 1000"
# #
# # # Paths to PRANK program for phylogenetic analysis
# # # PRANK_LECS = "/share/apps/bin/prank"
# # # PRANK = "/usr/local/bin/prank"
# # PRANK_LECS = "prank"
# # PRANK = "prank"
# #
# # # Path to T-Coffee program
# # T_COFFEE = "/share/apps/T-COFFEE-8.47/bin/binaries/linux/t_coffee"
# #
# # # Path to PAGAN program for phylogenetic analysis
# # PAGAN_LECS = "/share/apps/pagan-msa/bin/pagan"
# #
# # # Path to the tree viewer directory
# # TREE_VIEWER_DIR = "/bioseq/ConSurf_old/treeViewer/"
# #
# # # Path to PACC scripts directory
# # PACC_path = "/bioseq/ConSeq/external_scripts/PACC/"
# #
# # # Paths to different versions of Rate4Site program
# # RATE4SITE_BIOC_VER = "/bioseq/rate4site/BioCluster_Nov_06_dev/rate4site.exe"
# # RATE4SITE_SLOW_BIOC_VER = "/bioseq/rate4site/BioCluster_Nov_06_dev/rate4siteSlow.exe"
# # RATE4SITE = "/db1/Local/src/Rate4SiteSource/r4s_Nov_06_dev/rate4site.exe"
# # RATE4SITE_SLOW = "/db1/Local/src/Rate4SiteSource/r4s_Nov_06_dev/rate4siteSlow.exe"
# # RATE4SITE_SLOW_LECS = "/share/apps/bin/rate4site_slow"
# # RATE4SITE_LOCAL = "/bioseq/rate4site/rate4site"
# # RATE4SITE_SLOW_LOCAL = "/bioseq/rate4site/rate4site.doubleRep"
# # RATE4SITE_WITH_LG = "/bioseq/rate4site/With_LG/rate4site"
# # RATE4SITE_WITH_LG_SLOW = "/bioseq/rate4site/With_LG/rate4site.doubleRep"
# #
# # # Path to Ruby programming language executable
# # RUBY = "/share/apps/bin/ruby"  # or "/usr/bin/ruby"
# #
# # # Path to CD-HIT program
# # CD_HIT_DIR = "/bioseq/cd_hit/"
# #
# # # Paths to PACC prediction and MSA conversion scripts
# # PREDICT_PACC = "/bioseq/ConSeq/external_scripts/PACC/run.sh"
# # MSA_to_HSSP = "/bioseq/ConSeq/external_scripts/PACC/MSA2hssp.pl"
# #
# # # Path to Semphy program
# # SEMPHY = "/bioseq/Programs/Semphy/semphy.doubleRep"
# #
# # # Path to IQ_TREE program
# # # IQTREE = "/bioseq/Programs/iqtree/iqtree2"
# #
# # # Path to Semphy BBL program
# # # SEMPHY_BBL = ""
# #
# # # Path to Epitopia executable programs
# # # EPITOPIA_EXECUTABLES = "/bioseq/epitopia/executables"
# #
# # # Constant values
# # BLAST_MAX_HOMOLOGUES_TO_DISPLAY = 500
# # BLAST_PDB_MAX_HOMOLOGUES_TO_DISPLAY = 25
# # CONSURF_PIPE_FORM = "/bioseq/ConSurf_old/consurf_pipe.form"
# # SELECTON_MAX_NUCLEOTIDE = 15000
# # MAX_WALLTIME = "96:00:00"
# #
# # # Queue Details
# # BIOSEQ_NODE = "bioc01.tau.ac.il"
# # MAX_QUEUE_RUNS = 999
# #
# # # External links
# # # RCSB_WEB = "http://www.rcsb.org/"
# # # PYMOL_WEB = "http://pymol.sourceforge.net/"
# # # CHIMERA_WEB = 'http://www.rbvi.ucsf.edu/chimera/'
# # # CHIMERA_SAVING_FIGURE = 'http://www.cgl.ucsf.edu/chimera/current/docs/UsersGuide/print.html'
# # # CHIMERA_DOWNLOAD = CHIMERA_WEB + "download.html"
# # # MSA_CONVERT = 'http://www.ebi.ac.uk/cgi-bin/readseq.cgi'
# # # MSA_FORMATS = 'http://www.ebi.ac.uk/help/formats.html'
# # #
# # # # Redirect pages
# # # CONSURF_REDIRECT_PAGE = "http://consurf.tau.ac.il/too_many_runs.html"
# # # SELECTON_REDIRECT_PAGE = "http://selecton.tau.ac.il/too_many_runs.html"
# # # CONSEQ_REDIRECT_PAGE = "http://conseq.tau.ac.il/too_many_runs.html"
# # # PEPITOPE_REDIRECT_PAGE = "http://pepitope.tau.ac.il/too_many_runs.html"
# # #
# # # # FAQ pages
# # # CONSURF_TREE_FAQ = "http://consurf.tau.ac.il/quick_help.html#note5"
# # #
# # # # Files Name Conventions
# # # TEMPLATES_LIST_FILE = "List_of_Templates"
# # # PISA_ERRORS_FILE = "PISA_Errors"
