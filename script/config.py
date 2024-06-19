import os
import sys

Bin = os.path.dirname(sys.argv[0])
BIN_DIR = os.path.dirname(Bin)

# constants to use when sending e-mails using the server admin's email address.
ADMIN_EMAIL = "TAU BioSequence <bioSequence@tauex.tau.ac.il>"
ADMIN_USER_NAME = "bioSequence"
ADMIN_PASSWORD = ""
SMTP_SERVER = ""

# the name of the list of all running processes
QUEUING_JOBS = "/bioseq/bioSequence_scripts_and_constants/queuing_jobs.list"
RUNNING_JOBS = "/bioseq/bioSequence_scripts_and_constants/running_jobs.list"
SUBMITTED_JOBS = "/bioseq/bioSequence_scripts_and_constants/submitted_jobs.list"
JOBS_ON_BIOSEQ_NODE = "/bioseq/bioSequence_scripts_and_constants/jobs_on_bioc.01_node.list"
JOBS_WAITING_BIOSEQ_NODE = "/bioseq/bioSequence_scripts_and_constants/jobs_waiting_bioc.01_node.list"
CONSURF_RUNNING_JOBS = "/bioseq/bioSequence_scripts_and_constants/consurf_running_jobs.list"
SELECTON_RUNNING_JOBS = "/bioseq/bioSequence_scripts_and_constants/selecton_running_jobs.list"
CONSEQ_RUNNING_JOBS = "/bioseq/bioSequence_scripts_and_constants/conseq_running_jobs.list"
PEPITOPE_RUNNING_JOBS = "/bioseq/bioSequence_scripts_and_constants/pepitope_running_jobs.list"

# Databases urls
PROTEOPEDIA = "http://proteopedia.org/wiki/index.php/"
PDB_DB = "http://www.rcsb.org/pdb/explore/explore.do?structureId="
RCSB_WGET = "wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/"
RCSB = "http://www.rcsb.org/"
PISA_WGET = "wget http://www.ebi.ac.uk/msd-srv/pisa/cgi-bin/multimer.pdb?"

# CGIs paths
CONSURF_CGI_DIR = "/var/www/cgi-bin/ConSurf"

#general paths
SERVERS_RESULTS_DIR = "/bioseq/data/results/"
SERVERS_LOGS_DIR = "/bioseq/data/logs/"
SEND_EMAIL_DIR = "/bioseq/bioSequence_scripts_and_constants/sendEmail"
SEND_EMAIL_DIR_IBIS = "/bioseq/bioSequence_scripts_and_constants/sendEmail"
DAEMON_LOG_FILE = "/bioseq/bioSequence_scripts_and_constants/daemon.log"
UPDATE_RUN_TIME_LOG_FILE = "/bioseq/bioSequence_scripts_and_constants/update_runTime.log"
CONSURF_CGI = "/var/www/cgi-bin/ConSurf"
BIOSEQ_TEMP = "/bioseq/temp/"

# servers urls:
SELECTON_URL = "http://selecton.tau.ac.il"
CONSEQ_URL = "http://conseq.tau.ac.il/"
CONSURF_URL = "http://consurf.tau.ac.il/"
NEW_CONSURF_URL = "http://consurf.tau.ac.il/"
EPITOPIA_URL = "http://epitopia.tau.ac.il/"
PEPITOPE_URL = "http://pepitope.tau.ac.il/"
QMF_URL = "http://quasimotifinder.tau.ac.il/"
PATCHFINDER_URL = "http://patchfinder.tau.ac.il/"
FASTML_URL = "http://fastml.tau.ac.il/"
RECONST_URL = "http://fastml.tau.ac.il/reconst/"
GAIN_LOSS_URL = "http://gloome.tau.ac.il/"
CONSURF_DB_URL = "http://consurfdb.tau.ac.il/"
GILAD_SERVER_URL = "http://mud.tau.ac.il/"
MCPep_URL = "http://bental.tau.ac.il/MCPep/"
GUIDANCE_URL = "http://guidance.tau.ac.il/"
GUIDANCE_INDELS_URL = "http://guidance.tau.ac.il/indels/"
SPECBOOST_URL = "http://bental.tau.ac.il/specBoost/"
PROMAYA_URL = "http://bental.tau.ac.il/ProMaya/"
HOMOLOGY_SEARCH_URL = "http://fastml.tau.ac.il/HomologySearch/"
COPAP_URL = "http://copap.tau.ac.il/"

#servers logs:
CONSURF_LOG = "/bioseq/ConSurf_old/consurf.log"
CONSURF_NEW_LOG = "/bioseq/ConSurf/consurf.log"
SELECTON_LOG = "/bioseq/Selecton/selecton.log"
EPITOPIA_LOG = "/bioseq/epitopia/epitopia.log"
CONSEQ_LOG = "/bioseq/ConSeq/conseq.log"
PEPITOPE_LOG = "/bioseq/pepitope/pepitope.log"
RECONST_LOG = "/bioseq/ReConst_Server/reconst.log"
MCPep_LOG = "/bioseq/MCPep/mcpep.log"
Guidance_LOG = "/bioseq/Guidance/guidance.log"
Guidance_Indels_LOG = "/bioseq/GuidanceIndels/guidance_Indels.log"
MuD_LOG = "/bioseq/Gilad_Server/MuD.log"
FASTML_LOG = "/bioseq/FastML/fastml.log"
SPECBOOST_LOG = "/bioseq/specBoost/specBoost.log"
GAIN_LOSS_LOG = "/bioseq/GainLoss/GainLoss.log"
PROMAYA_LOG = "/bioseq/ProMaya/ProMaya.log"
COPAP_LOG = "/bioseq/CoPAP/CoPAP.log"

#servers results urls:
# servers urls:
SELECTON_RESULTS_URL = SELECTON_URL + "/results/"

#external databases
PQS = "/biodb/PQS/"
PDB_DIVIDED = "/biodb/PDB/data/structures/divided/pdb/"
SWISSPROT_DB = "/biodb/BLAST/Proteins/swissprot"
UNIPROT_DB = "/biodb/BLAST/Proteins/uniprot"
CLEAN_UNIPROT_DB = "/biodb/BLAST/Proteins/clean_uniprot"
UNIREF90_DB = "/biodb/BLAST/Proteins/uniref90"
PDBAA_NCBI = "/biodb/BLAST/Proteins/pdbaa"
CULLED_PDB = "/groups/bioseq.home/HAIM/PDBAA/pdbaaent"
PDB_DUNBRACK = "/groups/bioseq.home/HAIM/PDBAA/pdbaa"
NR_PROT_DB = "/biodb/BLAST/Proteins/nr"
NR_NUC_DB = "/biodb/BLAST/Nucleotides/nt"
UNIPROT_DAT_INDEX = "/bioseq/data/results/GB_CDS/uniprot.dat.bp_index"
PDB_TO_UNIPROT = "/bioseq/data/results/PDB_to_UNIPROT/idmapping_PDB_UNIPROTKB.dat"
PDB_TO_UNIPROT_test = "/biodb/idmapping_PDB_UNIPROTKB.dat"


#internal databases
EPITOPIA_DATA = "/bioseq/epitopia/data"

#external programs

# Paths to different BLAST programs
BLASTALL = "/opt/bio/ncbi/bin/blastall"
BLASTPGP = "blastpgp"
CS_BLAST = "/share/apps/csblast-2.1.0-linux64/csblast_static"

# Paths to MUSCLE program
MUSCLE_LECS = "/share/apps/bin/muscle"
# MUSCLE = "/usr/local/bin/muscle"
# MUSCLE_3_6 = "/bioseq/Programs/muscle_3.6_from_BIOCLUSTER/muscle3.6/muscle"
# MUSCLE_LECS = "muscle"
MUSCLE = "muscle"

# Paths to different versions of ClustalW program
CLUSTALW_LECS = "/share/apps/bin/clustalw"
CLUSTALW = "/usr/local/bin/clustalw"
CLUSTALW_1_82 = "/bioseq/Programs/ClustalW_1.82/clustalw1.82/clustalw"
CLUSTALW_1_81 = "/bioseq/Programs/ClustalW_1.81/clustalw1.81/clustalw"
CLUSTALW_2_0_10 = "/bioseq/Programs/ClustalW_2.0.10/clustalw-2.0.10-linux-i386-libcppstatic/clustalw2"
CLUSTAL_OMEGA = os.path.join(BIN_DIR, 'script/programs/clustalo')

# Paths to different versions of MAFFT program
MAFFT_LINSI = "/usr/local/bin/mafft-linsi"
MAFFT = "/usr/local/bin/mafft"
# MAFFT_GUIDANCE = "/bioseq/Programs/MAFFT_6.833/bin/mafft"
MAFFT_GUIDANCE = "mafft"
MAFFT_LINSI_GUIDANCE = "/bioseq/Programs/MAFFT_6.833/bin/mafft --localpair --maxiterate 1000"

# Paths to PRANK program for phylogenetic analysis
# PRANK_LECS = "/share/apps/bin/prank"
# PRANK = "/usr/local/bin/prank"
PRANK_LECS = "prank"
PRANK = "prank"

# Path to T-Coffee program
T_COFFEE = "/share/apps/T-COFFEE-8.47/bin/binaries/linux/t_coffee"

# Path to PAGAN program for phylogenetic analysis
PAGAN_LECS = "/share/apps/pagan-msa/bin/pagan"

# Path to the tree viewer directory
TREE_VIEWER_DIR = "/bioseq/ConSurf_old/treeViewer/"

# Path to PACC scripts directory
PACC_path = "/bioseq/ConSeq/external_scripts/PACC/"

# Paths to different versions of Rate4Site program
RATE4SITE_BIOC_VER = "/bioseq/rate4site/BioCluster_Nov_06_dev/rate4site.exe"
RATE4SITE_SLOW_BIOC_VER = "/bioseq/rate4site/BioCluster_Nov_06_dev/rate4siteSlow.exe"
RATE4SITE = "/db1/Local/src/Rate4SiteSource/r4s_Nov_06_dev/rate4site.exe"
RATE4SITE_SLOW = "/db1/Local/src/Rate4SiteSource/r4s_Nov_06_dev/rate4siteSlow.exe"
RATE4SITE_SLOW_LECS = "/share/apps/bin/rate4site_slow"
RATE4SITE_LOCAL = "/bioseq/rate4site/rate4site"
RATE4SITE_SLOW_LOCAL = "/bioseq/rate4site/rate4site.doubleRep"
RATE4SITE_WITH_LG = "/bioseq/rate4site/With_LG/rate4site"
RATE4SITE_WITH_LG_SLOW = "/bioseq/rate4site/With_LG/rate4site.doubleRep"

# Path to Ruby programming language executable
RUBY = "/share/apps/bin/ruby"  # or "/usr/bin/ruby"

# Path to CD-HIT program
CD_HIT_DIR = "/bioseq/cd_hit/"

# Paths to PACC prediction and MSA conversion scripts
PREDICT_PACC = "/bioseq/ConSeq/external_scripts/PACC/run.sh"
MSA_to_HSSP = "/bioseq/ConSeq/external_scripts/PACC/MSA2hssp.pl"

# Path to Semphy program
SEMPHY = "/bioseq/Programs/Semphy/semphy.doubleRep"

# Path to IQ_TREE program
# IQTREE = "/bioseq/Programs/iqtree/iqtree2"

# Path to Semphy BBL program
SEMPHY_BBL = ""

# Path to Epitopia executable programs
EPITOPIA_EXECUTABLES = "/bioseq/epitopia/executables"

# Constant values
BLAST_MAX_HOMOLOGUES_TO_DISPLAY = 500
BLAST_PDB_MAX_HOMOLOGUES_TO_DISPLAY = 25
CONSURF_PIPE_FORM = "/bioseq/ConSurf_old/consurf_pipe.form"
SELECTON_MAX_NUCLEOTIDE = 15000
MAX_WALLTIME = "96:00:00"

# Queue Details
BIOSEQ_NODE = "bioc01.tau.ac.il"
MAX_QUEUE_RUNS = 999

# External links
RCSB_WEB = "http://www.rcsb.org/"
PYMOL_WEB = "http://pymol.sourceforge.net/"
CHIMERA_WEB = 'http://www.rbvi.ucsf.edu/chimera/'
CHIMERA_SAVING_FIGURE = 'http://www.cgl.ucsf.edu/chimera/current/docs/UsersGuide/print.html'
CHIMERA_DOWNLOAD = CHIMERA_WEB + "download.html"
MSA_CONVERT = 'http://www.ebi.ac.uk/cgi-bin/readseq.cgi'
MSA_FORMATS = 'http://www.ebi.ac.uk/help/formats.html'

# Redirect pages
CONSURF_REDIRECT_PAGE = "http://consurf.tau.ac.il/too_many_runs.html"
SELECTON_REDIRECT_PAGE = "http://selecton.tau.ac.il/too_many_runs.html"
CONSEQ_REDIRECT_PAGE = "http://conseq.tau.ac.il/too_many_runs.html"
PEPITOPE_REDIRECT_PAGE = "http://pepitope.tau.ac.il/too_many_runs.html"

# FAQ pages
CONSURF_TREE_FAQ = "http://consurf.tau.ac.il/quick_help.html#note5"

# Files Name Conventions
TEMPLATES_LIST_FILE = "List_of_Templates"
PISA_ERRORS_FILE = "PISA_Errors"
