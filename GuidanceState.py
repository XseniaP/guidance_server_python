#module load python/python-3.6.7

import os
import shutil
import sys
import stat
import cgi
import cgitb
import logging
import subprocess
import socket
from storable import retrieve
from time import time, ctime, sleep
from random import randint
import re
import json
from InputValidator import InputValidator
from utils import *

Bin = os.path.dirname(sys.argv[0])
BIN_DIR = os.path.dirname(Bin)

# if os.path.exists('/home/josefspr/bioseq'):  # remote run
#     sys.path.insert(0, '/home/josefspr/bioseq/guidance/guidance.v2.02/www/Guidance')
#     sys.path.insert(1, '/home/josefspr/bioseq/bioSequence_scripts_and_constants')
# else:
#     sys.path.insert(0, '/Users/kpolonsky/Documents/GUIDANCE-guidance.v2.02/www/Guidance')
#     # sys.path.insert(0, sys.path.join(BIN_DIR, '')
#     sys.path.insert(1, '/Users/kpolonsky/Documents/GUIDANCE-guidance.v2.02/www/bioSequence_scripts_and_constants')
#     # sys.path.insert(1, '/Users/kpolonsky/Documents/GUIDANCE-guidance.v2.02/www/bioSequence_scripts_and_constants')

import SharedConsts as CONSTS  

#from /home/josefspr/bioseq/bioSequence_scripts_and_constants/
#from email_sender import send_email  
#def send_email(smtp_server, sender, receiver, subject, content):
#    return
                   
class GuidanceState:

    def __init__(self, jobId: str, form: dict = None, files = [], isRequest = True):
        
        if form != None: 
            try: 

                # working_dir, output url
                results_url = os.path.join(CONSTS.WEBSERVER_RESULTS_URL, jobId)
                output_url = os.path.join(results_url, CONSTS.RESULT_WEBPAGE_NAME)
                wd = os.path.join(CONSTS.WEBSERVER_RESULTS_DIR, jobId)
                os.makedirs(wd)
                
                # setup logger
                job_logger = get_job_logger(jobId)
                job_logger.info(f'{"#" * 100}\n{currentTime()}: entering initial GuidanceState.__init__\n')
                
                # initialize var dictionary
                var = {}
                var['WorkingDir'] = wd
                var['SeqsFile'] = 'Seqs.Orig.fas' #generic fixed name for the sequence file. the user file will be copied to this file.
                var['SeqsFile_Codons'] = 'Seqs.Orig_DNA.fas' #generic fixed name for the DNA CODONS sequence file. the user file will be copied to this file.
                var['dataset'] = 'MSA'
                #if dict_file_defined_not_empty('userMSA_File', files): 
                if dict_key_defined_not_empty('BACK_FROM_MAFFT', form):
                    var['Alignment_File'] = 'UserMSA'
                else: 
                    var['Alignment_File'] = f"{var['dataset']}.{form['MSA_Program']}.aln"
                var['align_param'] = ''
                var['run_url'] = results_url + '/'
                var['output_page'] = CONSTS.RESULT_WEBPAGE_NAME
                var['run_number'] = jobId
                var['code_fileName'] = 'Seqs.Codes'
                
                # write paramters to log file
                peek_form(form, files, job_logger) 
                
                if isRequest: 
                    self.form = form.to_dict()
                else: 
                    self.form = form
                
                if dict_key_defined_not_empty('BACK_FROM_MAFFT', form):
                    self.form['Redirect_From_MAFFT'] = form['BACK_FROM_MAFFT']
                    del self.form['BACK_FROM_MAFFT']
                else:
                    self.form['Redirect_From_MAFFT'] = '0'

                
                self.files = files
                self.var = var
                
            except:
            
                raise Exception ('GuidanceState.__init__ failed', 'system')
            
        else:
    
            try: 
                
                wd = os.path.join(CONSTS.WEBSERVER_RESULTS_DIR, jobId)
                
                job_logger = get_job_logger(jobId)
                job_logger.info(f'{"#" * 100}\nGuidanceState.__init__: jobId = {jobId}, wd = {wd}')
                
                if not os.path.exists (wd):
                    raise Exception ('GuidanceState.__init__: jobId does not exist', 'system')
                
                form_path = os.path.join( wd, 'FORM.json')
                with open (form_path, 'r') as f:
                    self.form = json.load (f)
                f.close()
                
                var_path = os.path.join( wd, 'VARS.json')
                with open (var_path, 'r') as f:
                    self.var = json.load (f)
                f.close()
                 
                self.files = {}
                
            except:
                
                raise Exception (f'GuidanceState.__init__: failed to load job {jobId}', 'system')
    
    def write_hello_world(self):
        
        output_path = os.path.join(self.var['WorkingDir'], CONSTS.RESULT_WEBPAGE_NAME)
        with open(output_path, "w") as f:
            f.write('Hello World Python<br>\n')
            f.write(f"Run number {self.var['run_number']}<br>\n")
            for key in sorted(self.form.keys()):
                f.write(f'{key} = {self.form[key]}<br>\n')
        f.close()
       
    def save_state(self, warning_messages = '', crash_flag = False): 
    
        try: 
        
            #FORM
            cgi_form = self.form
            files = self.files
            var = self.var
            
            form = {}
            form['MSA_Program'] = cgi_form['MSA_Program']
            form['JOB_TITLE'] = cgi_form['JOB_TITLE']
            form['Bootstraps'] = int(cgi_form['Bootstraps'])
            form['user_mail'] = cgi_form['email_add']
            form['PROGRAM'] = cgi_form['PROGRAM']
            form['Redirect_From_MAFFT'] = cgi_form['Redirect_From_MAFFT']

            form['Seq_Type'] = cgi_form['Seq_Type']
            
            if form['Redirect_From_MAFFT'] == '0':
                form['FASTA_Seqs_Text'] = cgi_form['FASTA_txt']
            form['CALLING_SERVER']= 'GUIDANCE2'
            form['Run_Number'] = var['run_number']
            
            if dict_key_defined_not_empty('MAFFT_OUT_RUN_NEMBER', cgi_form):
                form['MAFFT_RUN_OUTPUT_NAME'] = cgi_form['MAFFT_OUT_RUN_NEMBER']
                
            if dict_key_defined_not_empty('MAFFT_ALIGN', cgi_form):
                form['MAFFT_ALIGN'] = cgi_form['MAFFT_ALIGN']
            
            if dict_key_defined_not_empty('userMSA_File', cgi_form):
                form['userMSA_File'] = cgi_form['userMSA_File']
            
            if dict_key_defined_not_empty('outorder', cgi_form):
                form['Align_Order'] = cgi_form['outorder']
                
            if dict_key_defined_not_empty('GENCODE', cgi_form):
                form['CodonTable'] = cgi_form['GENCODE'] 
            
            if dict_key_defined_not_empty('F', cgi_form):
                form['PRANK_F'] = cgi_form['F']
            
            if dict_key_defined_not_empty('maxiterate', cgi_form):
                form['MAFFT_maxiterate'] = cgi_form['maxiterate']
            
            if dict_key_defined_not_empty('pair', cgi_form):
                form['MAFFT_refinement'] = cgi_form['pair'] 
            else:
                form['MAFFT_refinement'] =''
                
            if dict_key_defined_not_empty('RERUN_SAME_SEQ', cgi_form):
                form['RERUN_SEQ_ONLY'] = cgi_form['RERUN_SAME_SEQ']
            else:
                form['RERUN_SEQ_ONLY'] = None
            
            #Handle and validate rerun for same seq
            if dict_key_defined_not_empty('RERUN_SAME_SEQ', cgi_form):
                old_run_dir = os.path.join(CONSTS.WEBSERVER_RESULTS_DIR, form['RERUN_SEQ_ONLY'])
                if os.path.isdir(old_run_dir):
                    form['RERUN_NO_CAPTCHA']="YES"
                else:
                    form['RERUN_NO_CAPTCHA']="NO"
            else:
                form['RERUN_NO_CAPTCHA']="NO"
            
            if dict_file_defined_not_empty('usrSeq_File', files):
                form['usrSeq_File'] = files['usrSeq_File'].filename
            else:
                form['usrSeq_File'] = None
            
            #    form['userMSA_File'] = files['MSAFile'].filename
            #else:
                
                
            form_path = os.path.join( var['WorkingDir'], "FORM.json")
            with open(form_path, 'w') as fp:
                json.dump(form, fp)
            fp.close()
        
        except:
            return f'GuidanceState.store elements: storing form failed'
            
        try:
            #VARS
            
            # defaults
            var['SP_SEQ_CUTOFF']="0.6"
            var['SP_COL_CUTOFF']="0.93"
            var['warning_messages'] = warning_messages
            ### Sending e-mails
            var['send_email_dir'] = CONSTS.SEND_EMAIL_DIR_IBIS
            var['smtp_server'] = CONSTS.SMTP_SERVER
            var['userName'] = CONSTS.ADMIN_USER_NAME
            var['userPass'] = CONSTS.ADMIN_PASSWORD
            var['IsSPAM']=0
            var['proc_num']=1
            if dict_key_defined_not_empty( 'LongestSeq', var):
                var['LongestSeq']=str(var['LongestSeq'])

            # Alignment Parameters
            if 'align_param' not in var:
                var['align_param'] = ''
            if form['MSA_Program'] == 'MAFFT':
                if dict_key_defined_not_empty('MAFFT_maxiterate', form):
                    if int(form['MAFFT_maxiterate']) > 0:
                        var['align_param'] = var['align_param'] + f" --maxiterate {form['MAFFT_maxiterate']}"
                if dict_key_defined_not_empty('MAFFT_refinement', form):
                    var['align_param'] = var['align_param'] + f" --{form['MAFFT_refinement']}"
                if dict_key_defined_not_empty('Align_Order', form):
                    if form['Align_Order'] == 'aligned':
                        var['align_param'] = var['align_param'] + " --reorder"
            elif form['MSA_Program'] == 'PRANK':
                 var['align_param'] = var['align_param'] + form['PRANK_F']
            
            if form['JOB_TITLE'] == 'daily_test':
                var['write_daily_test_cmd'] = f"python {CONSTS.WRITE_DAILY_TEST_SCRIPT} {CONSTS.DAILY_TEST_DIR} {var['run_number']}"
                
            var['OutLogFile'] = os.path.join( CONSTS.WEBSERVER_LOGS_DIR, f"{var['run_number']}.log") 
            var['WorkingDir'] = var['WorkingDir'] + '/'
            var['state'] = 'SUBMIT'
            
            if crash_flag:
                var['crashFlag'] = '1'
                
            vars_path = os.path.join( var['WorkingDir'], "VARS.json")
            with open(vars_path, 'w') as fp:
                json.dump(var, fp)
            fp.close()
            
            # convert to perl hash files
            # convertScript = '/home/josefspr/bioseq/guidance/guidance.v2.02/www/Guidance/json2hash.pl'
            # convertScript = '/Users/kpolonsky/Documents/GUIDANCE-guidance.v2.02/www/Guidance/json2hash.pl'
            convertScript = os.path.join(BIN_DIR, "script", 'json2hash.pl')
            form_data_path = os.path.join( var['WorkingDir'], 'form.data')
            var_data_path = os.path.join( var['WorkingDir'], 'input.data')
            cmd = f'perl {convertScript} {form_path} {form_data_path}'
            subprocess.run(cmd, shell=True)
            cmd = f'perl {convertScript} {vars_path} {var_data_path}'
            subprocess.run(cmd, shell=True)
            
        except:
            return f'GuidanceState.store elements: storing var failed'
        
        return "OK"
    
    def upload_files(self):
        
        form = self.form
        files = self.files
        var = self.var
        
        job_logger = logging.getLogger(var['run_number'])
        
        if not dict_key_value('Redirect_From_MAFFT','1', form): # regular run
        
            # upload data part
            job_logger.info(f'\n{"#" * 80}\nUploading data regular run\n')
            
            # upload sequence 
            seqsFile = os.path.join( var['WorkingDir'], var['SeqsFile'])
            alignment_file = os.path.join( var['WorkingDir'], var['Alignment_File'])
            alignment_file_not_empty = False
            if os.path.exists(alignment_file): 
                if os.path.getsize(alignment_file) > 0:
                    alignment_file_not_empty = True
            
            if form['FASTA_txt']:
                job_logger.info(f'writing to {seqsFile}\n')
                try: 
                    with open(seqsFile, "w") as f_seq:
                        f_seq.write(form['FASTA_txt'])
                    f_seq.close()
                except: 
                    error= f'GuidanceState.upload_files: can\'t open {seqsFile} for writing.'
                    raise Exception(error, "system")
                var['IsSPAM'] = checkForSpam(seqsFile, var['WorkingDir']);
            else:
                # upload user seqs file and alignment file
                try: 
                    if not os.path.exists(seqsFile) and files['usrSeq_File'].filename:
                        uploading_file = files['usrSeq_File'].filename
                        save_file_to_disk(var['WorkingDir'], job_logger, files['usrSeq_File'], var['SeqsFile'])
                    if not os.path.exists(alignment_file) and dict_key_defined_not_empty('userMSA_File', form):
                        uploading_file = files['userMSA_File'].filename
                        save_file_to_disk(var['WorkingDir'], job_logger, files['userMSA_File'], var['Alignment_File'])
                except: 
                    error= f'GuidanceState.upload_files: uploading file {uploading_file} failed.'
                    raise Exception(error, "system")

            # rename to codons file
            if form['Seq_Type'] == 'Codons': # and ('userMSA_File','') in form.items(): 
                if os.path.exists (seqsFile):
                    shutil.move( seqsFile, os.path.join( var['WorkingDir'], var['SeqsFile_Codons']))
                else:
                    error= f'GuidanceState.upload_files: {seqsFile} does not exist.'
                    raise Exception(error, "system")
                
            job_logger.info(f"ls of {var['WorkingDir']} fields:\n{os.listdir(var['WorkingDir'])}\n")
            
        else: # from mafft
            
            var['align_param'] = form['MAFFT_PARAM']
            form['MAFFT_RUN_OUTPUT_NAME'] = form['MAFFT_OUT_RUN_NEMBER']

            job_logger.info (f"====== THIS RUN WAS REDIRECTED FROM MAFFT SERVER:  http://mafft.cbrc.jp/alignment/server/spool/{form['MAFFT_RUN_OUTPUT_NAME']}.html ========")
            job_logger.info (f"====== MAFFT ORIGINAL RUNNING ARGUMENTS: {var['align_param']} =====")

            m = re.search('(.*)\s+input', var['align_param'])
            if m:
                 var['align_param'] = m.group(1)
                
            form['MAFFT_ALIGN'] = form['FASTA_txt'].replace('\r\n', '\n')
            form['userMSA_File']="Alignment from MAFFT"
            lines = re.split( r'>.*?\n', form['MAFFT_ALIGN']) 
            AA_Lines=0
            Nuc_Lines=0
            for line in lines: 
                line = line.rstrip()
                if line == '':
                    continue
                if line[0] != '>': 
                    line= line.replace ('\n', '')
                    if re.search( '^([ACTGUNRYSWKMBDHV\-actgunryswkmbdhv]+)$', line):  # including IUPAC: http://www.bioinformatics.org/sms/iupac.html
                        Nuc_Lines =+ 1
                    else:
                        AA_Lines += 1
                    
            if AA_Lines>0 and AA_Lines>Nuc_Lines:
                form['Seq_Type'] = 'AminoAcids'
            else:
                form['Seq_Type'] = 'Nucleotides'
            
            job_logger.info( f"Found {AA_Lines} AA lines and {Nuc_Lines} Nucleotides -  Decision:{form['Seq_Type']}")
            
        self.form = form
        self.var = var
            
    def update_state(self, state, error_msg = '', error_type = ''):
        # currenly there are 3 states: Init, Finished, and Error/Crashed
        form = self.form
        var = self.var
        warning_msg = ''
        
        job_logger = logging.getLogger(var['run_number'])
        job_logger.info(f'update_state({state.str()})')
        
        if 'state' in var.keys():
            if var['state'] == 'FINISHED' or var['state'] == 'ERROR':
                return
        var['state'] = state.str()
            
        if state == State.Init:
            files = self.files
             
            if form['Redirect_From_MAFFT'] == '1': 
    
                # align param
                if 'anysymbol' in var['align_param']: 
                    var['align_param'] = var['align_param'].replace('--anysymbol', '')
                    warning_msg = "<br><b><font color=\"red\" size=\"3\">Warnning:</b></font><font size=\"3\"> --anysymbol is not allowed in GUIDANCE, therefore if there are non-standard characters your run is expected to fail, relevant message will issue, then please fix your input and resubmit, otherwise GUIDANCE run continues</font></br>\n"

                if 'seed' in var['align_param']:

                    MAFFT_RunNumber = ''
                    job_logger.info(f"MAFFT_RunNumber = {form['MAFFT_RUN_OUTPUT_NAME']}")
                    
                    m = re.search( '_out([0-9]+)', form['MAFFT_RUN_OUTPUT_NAME'])
                    if m:
                       MAFFT_RunNumber = m.group(1)
                       
                    seeds = re.findall( '\-\-seed (str[0-9]+)', var['align_param'])
                    SeedsCount = len(seeds)
                    for i in range (0, SeedsCount-1):

                        seed_on_MAFFT = f"http://mafft.cbrc.jp/alignment/server/tmp/_seed{i}.var['MAFFT_RunNumber']"
                        job_logger.info( f"SEED:{seed_on_MAFFT}")
                        wget_cmd = f"links -source {seed_on_MAFFT} > {var['WorkingDir']}seedFile{i}"
                        job_logger.info(wget_cmd)
                        subprocess.run(wget_cmd, shell=True)
                        var['align_param'] = var['align_param'].replace( f'--seed {seeds[i]}', f"--seed {var['WorkingDir']}seedFile{seeds[i]}", var['align_param'])

                    job_logger.info( f'Total: {SeedsCount}')
                    if '--reorder' in var['align_param']:
                        var['align_param'] = var['align_param'].replace('--reorder', '') # if seed is provided reorder must be removed so the seeds will be first
                        job_logger.info ("WARNNING: --reorder is not allowed if seed alignment is provided, therefore the --reorder argument will be ignored and the output order will be same as input (with seeds first)")
                        warning_msg = "<br><b><font color=\"red\" size=\"3\">Warnning:</b></font><font size=\"3\"> --reorder is not allowed if seed alignment is provided, therefore the --reorder argument will be ignored and the output order will be same as input (with seed(s) first)</font></br>\n"
                        

                    wget_cmd = f"links -source http://mafft.cbrc.jp/alignment/server/tmp/_inx$VARS{var['MAFFT_RunNumber']} > {os.path.join (var['WorkingDir'], var['SeqsFile'])}" # we need the 'core' sequences without seed for the alignment
                    subprocess.run(wget_cmd, shell=True)
                    job_logger.info (wget_cmd)
                    try:
                        with open ( os.path.join( var['WorkingDir'], 'ALIGNMENT_FROM_MAFFT'), 'w') as f:
                            f.write(form['MAFFT_ALIGN'])
                        f.close()
                    except: 
                        error = f'GuidanceState.update_state(Init): failed to write to file ALIGNMENT_FROM_MAFFT'
                        raise Exception ( error, "system")

                else: # Really take the alignment only if seed not provided, otherwise we build the alignment again to make it simpler
                    try:
                        with open ( os.path.join(var['WorkingDir'],var['Alignment_File']), 'w') as f:
                            f.write(form['MAFFT_ALIGN'])
                        f.close()
                        job_logger.info (f"update_state(INIT): write input to {os.path.join(var['WorkingDir'],var['Alignment_File'])}")
                    except:
                        error = f"GuidanceState.update_state(Init): failed to write to file {var['Alignment_File']}"
                        raise Exception ( error, "system")
                    
                # KSENIA !!! DELETE THIS LINE
                # shutil.copy('/Users/kpolonsky/PycharmProjects/guidance_server_python/test_MAFFT_server_MSA_input4_redirectedFromMafft.txt', os.path.join(var['WorkingDir'],'ALIGNMENT_FROM_MAFFT'))
                # form['userMSA_File'] = 'ALIGNMENT_FROM_MAFFT'

                form['userMSA_File'] = 'Alignment from MAFFT'
                
                if 'qinsi' in var['align_param']:
                    error =  "GuidanceState.update_state(Init): We are sorry but Q-INS-i methodology is currently not supported by GUIDANCE2. Sorry for the inconvenience.<br>"
                    raise Exception ( error, "user")

                matches = re.findall('(allowshift)|(unalignlevel)|(leavegappyregion)|(regtable2seq)|(adjustdirection)|(dash)', var['align_param'])
                if len(matches) > 0: 
                    args = ''
                    for arg in matches:
                        args = args + ' ' + arg
                    job_logger.info( f"We are sorry but '$arg' option of MAFFT is currently not supported by GUIDANCE2. Sorry for the inconvenience.\nPARAMS: {var['align_param']}")
                    error =  f"GuidanceState.update_state(Init): We are sorry but $arg option of MAFFT is currently not supported by GUIDANCE2. Sorry for the inconvenience. Please feel free to <a href=\"mailto:evolseq\@taux.tau.ac.il\?subject=GUIDANCE2\%20Run\%20Number\%20{var['run_number']}\">contact us</a> about it<br>"
                    raise Exception ( error, "user")
        
                # unify case of input file
                if form['PROGRAM'] == "GUIDANCE":

                    if form['Seq_Type'] == "Nucleotides":
                        ans = convert_fs_to_lower_case(os.path.join(var["WorkingDir"],var["Alignment_File"])) # MAFFT ALWAYS OUTPUT NUC MSA IN LOWER CASE
                        if ans != 'OK':
                            raise Exception ( ans, "system")
                    if form['Seq_Type'] == "AminoAcids":
                        ans = convert_fs_to_upper_case(os.path.join(var["WorkingDir"],var["Alignment_File"])) # MAFFT ALWAYS OUTPUT AA  MSA IN UPPER CASE
                        if ans != 'OK':
                            raise Exception ( ans, "system")
                            
                else:# GUIDANCE2 or HoT
                    if  form['Seq_Type'] == "Nucleotides":
                        ans = convert_fs_to_upper_case(os.path.join(var["WorkingDir"],var["Alignment_File"])) # HoT and GUIDANCE2 ALWAYS OUTPUT NUC MSA IN UPPER CASE
                        if ans != 'OK':
                            raise Exception ( ans, "system")
                            
            # number of sequences analized
            if not dict_key_value('Redirect_From_MAFFT','1', form): #case regular run
                if dict_key_defined_not_empty( 'userMSA_File', form):
                    var['NumOfSeq'] = InputValidator.countSeq(os.path.join(var['WorkingDir'],form['userMSA_File']))
                else:
                    if form['Seq_Type'] != 'Codons':
                        var['NumOfSeq'] = InputValidator.countSeq(os.path.join(var['WorkingDir'],var['SeqsFile']))
                    else:
                        var['NumOfSeq'] = InputValidator.countSeq(os.path.join(var['WorkingDir'],var['SeqsFile_Codons']))
            else:
                var['NumOfSeq'] = '' 

            # sequences/msa link
            results_url = os.path.join( CONSTS.WEBSERVER_RESULTS_URL, var['run_number'])
            
            if dict_key_defined_not_empty( 'userMSA_File', form):
                if form['Redirect_From_MAFFT'] == '1' and 'seed' in var['align_param']: 
                    var['sequences_link'] = f"MSA File = <A HREF='{results_url}/ALIGNMENT_FROM_MAFFT' TARGET=_blank>{form['userMSA_File']}</A><br>"
                else:
                    var['sequences_link'] = f"MSA File = <A HREF='{results_url}/UserMSA.FIXED.ORIG' TARGET=_blank>{form['userMSA_File']}</A><br>"
            else:
                if form['Seq_Type'] == 'Codons':
                    if dict_file_defined_not_empty('usrSeq_File', files):
                        var['sequences_link'] = f"Sequences File = <A HREF='{results_url}/{var['SeqsFile_Codons']}' TARGET=USER_Seqs>{files['usrSeq_File'].filename}</A><br>"
                    elif form['JOB_TITLE']:
                        var['sequences_link'] = f"Sequences = <A HREF='{results_url}/{var['SeqsFile_Codons']}' TARGET=USER_Seqs>{form['JOB_TITLE']}</A><br>"
                    else:
                        var['sequences_link'] = f"Sequences = <A HREF='{results_url}/{var['SeqsFile_Codons']}' TARGET=USER_Seqs>Fasta Sequences</A><br>"
                    
                else: #NOT CODONS
                    if dict_file_defined_not_empty('usrSeq_File', files):
                        var['sequences_link'] = f"Sequences File = <A HREF='{results_url}/{var['SeqsFile']}' TARGET=USER_Seqs>{files['usrSeq_File'].filename}</A><br>"
                    elif form['JOB_TITLE']:
                        var['sequences_link'] = f"Sequences = <A HREF='{results_url}/{var['SeqsFile']}' TARGET=USER_Seqs>{form['JOB_TITLE']}</A><br>"
                    else:
                        var['sequences_link'] = f"Sequences = <A HREF='{results_url}/{var['SeqsFile']}' TARGET=USER_Seqs>Fasta Sequences</A><br>"
                                
            var['dataset']='MSA'
            var['res_pair_res_html_file']=f"{var['run_number']}/{var['dataset']}.{form['MSA_Program']}.Guidance_res_pair_res.html"
            var['Alignment_File_With_Names']= f"{var['run_number']}/{var['Alignment_File']}.With_Names"
            if form['PROGRAM'] == "GUIDANCE" or form['PROGRAM'] == "HoT":
                var['Output_Prefix']=f"{var['dataset']}.{form['MSA_Program']}.Guidance"
            elif form['PROGRAM'] == "GUIDANCE2":
                var['Output_Prefix']=f"{var['dataset']}.{form['MSA_Program']}.Guidance2"
            var['Seq_Scores']= f"{var['run_number']}/{var['Output_Prefix']}_res_pair_seq.scr_with_Names"
        
            if form['Seq_Type'] == 'AminoAcids' or form['Seq_Type'] == 'Codons':
                var['type_a'] = 'aa'
            elif form['Seq_Type'] ==  'Nucleotides': 
                var['type_a'] = 'nuc'
            
        elif state == State.Finished:
        
            var['Alignment_File_without_low_SP_Col'] = f"{var['dataset']}.{form['MSA_Program']}.Without_low_SP_Col"
            var['Alignment_File_without_low_SP_Col_with_Names'] = f"{var['Alignment_File_without_low_SP_Col']}.With_Names"
            var['removed_low_SP_SITE'] = f"{var['SeqsFile']}.{var['dataset']}.{form['MSA_Program']}.Removed_Col"
            
            var['Seq_File_without_low_SP_SEQ']=f"{var['SeqsFile']}.Without_low_SP_Seq"
            var['removed_low_SP_SEQ']=f"{var['SeqsFile']}.Removed_Seq"
            var['Seq_File_without_low_SP_SEQ_with_Names']=f"{var['Seq_File_without_low_SP_SEQ']}.With_Names"
            var['removed_low_SP_SEQ_With_Names']=f"{var['removed_low_SP_SEQ']}.With_Names"
            
            ''' # update sequences link
            if dict_key_defined_not_empty( 'userMSA_File', form):
                if form['Redirect_From_MAFFT'] == '1' and 'seed' not in var['align_param']: 
                    var['sequences_link'] = f"MSA File = <A HREF='{results_url}/UserMSA.FIXED.ORIG' TARGET=_blank>{form['userMSA_File']}</A><br>"
            '''
            
            wd = os.path.join(CONSTS.WEBSERVER_RESULTS_DIR, var['run_number'])
            residue_scores_path = os.path.join(wd, f"{var['Output_Prefix']}_res_pair_res.scr")
            if os.path.exists(residue_scores_path):
                var['residue_scores'] = f"{var['run_number']}/{var['Output_Prefix']}_res_pair_res.scr"
            else:
                var['residue_scores'] = f"{var['run_number']}/{var['Output_Prefix']}_res_pair_res.zip"
                
            residue_pair_scores_path = os.path.join(wd, f"{var['Output_Prefix']}_res_pair.scr")
            if os.path.exists(residue_scores_path):
                var['residue_pair_scores'] = f"{var['run_number']}/{var['Output_Prefix']}_res_pair.scr"
            else:
                var['residue_pair_scores'] = f"{var['run_number']}/{var['Output_Prefix']}_res_pair.zip"
            
            MSA_score_file = os.path.join (var['WorkingDir'] , f"{var['Output_Prefix']}_msa.scr")
            var['mean_res_pair_score'] = 0
            with open (MSA_score_file, 'r') as f: 
                for line in f: 
                    m = re.search('^\#MEAN_RES_PAIR_SCORE ([0-9.]+)', line)
                    if m:
                        var['mean_res_pair_score'] = m.group(1)
            f.close()
            
            var['select_remove_site_selection_mask'] = print_remove_site_selection_mask(residue_scores_path)
            
            with open ( os.path.join(var['WorkingDir'], 'MSA_LENGTH'), 'r') as f:
                msa_length = f.read()
            f.close()
            column_pair_scores_path = os.path.join(wd, f"{var['Output_Prefix']}_res_pair_col.scr" )
            var['select_remove_site_selection_box'] = print_remove_site_selection_box(column_pair_scores_path, int(msa_length))
            
            if not file_exists_not_empty( os.path.join( var['WorkingDir'], var['removed_low_SP_SITE']) ): 
                var['remove_columns_default_line'] = 'All positions had score higher than ' + "{:.3f}".format(float(var['SP_COL_CUTOFF'])) + '<br>'
            elif not os.path.exists( os.path.join( var['WorkingDir'],  var['Alignment_File_without_low_SP_Col_with_Names']) ):
                # filtered is empty file -> all positions were removed...
                var['remove_columns_default_line'] = '<font color=\'red\'><B>ATTENTION:</font></B> All positions had score below ' + "{:.3f}".format(float(var['SP_COL_CUTOFF'])) + '<br>'
            else: # give the filtered file done with default values
                var['remove_columns_default_line'] = f"<li><A HREF=\'{os.path.join( var['run_number'],  var['Alignment_File_without_low_SP_Col_with_Names'])}\' TARGET=_blank>The MSA after the removal of unreliable columns (below {var['SP_COL_CUTOFF']})</A><font size=-1> (see list of removed columns <A HREF=\'{os.path.join( var['run_number'], var['removed_low_SP_SITE'])}\' TARGET=_blank>here</A>)</font></li>"
            
            msa_depth_path = os.path.join(var['WorkingDir'], 'MSA_DEPTH')
            if not os.path.exists(msa_depth_path):
                sleep (3)
            
            if os.path.exists(msa_depth_path):
                with open ( msa_depth_path, 'r') as f:
                    msa_depth = f.read()
                f.close()
              
                # only create if msa depth file exists
                seq_pair_scores_path = os.path.join(wd, f"{var['Output_Prefix']}_res_pair_seq.scr" )
                var['select_remove_seq_selection_box'] = print_remove_seq_selection_box(seq_pair_scores_path, int(msa_depth))
            
            if not os.path.exists( os.path.join( var['WorkingDir'], var['Seq_File_without_low_SP_SEQ_with_Names']) ): 
                var['remove_seq_default_line'] = f"<font color='red'><B>ATTENTION:</font></B> All sequences had score below {var['SP_SEQ_CUTOFF']}<br>"
            elif not os.path.exists( os.path.join( var['WorkingDir'],  var['removed_low_SP_SEQ_With_Names']) ):
                var['remove_seq_default_line'] = f"All sequences had score higher than {var['SP_SEQ_CUTOFF']}<br>"
            else: 
                if form['Redirect_From_MAFFT'] == '1':
                    var['remove_seq_default_line'] = f'''
                        <A HREF="{var['Seq_File_without_low_SP_SEQ_with_Names']}" TARGET=_blank>The input sequences after the removal of unreliable sequences (with confidence score below {var['SP_SEQ_CUTOFF']})</A>
                        <font size=-1> (see list of removed sequences <A HREF="{var['removed_low_SP_SEQ_With_Names']}"
                        TARGET=_blank>here</A></font>)&nbsp;&nbsp;&nbsp;
                        <INPUT TYPE=\"BUTTON\" VALUE=\"run GUIDANCE on the confidently-aligned sequences only\" 
                        ONCLICK=\"var answer = confirm('ATTENTION: Running GUIDANCE on the confidently-aligned sequences only, ignores the parameters used for the original run on MAFFT server. It is therefore recommended to adjust these parameters or aligning the confidently-aligned sequences on MAFFT server and run GUIDANCE again from there');
                        if (answer) {{window.open(\'{CONSTS.WEBSERVER_URL}/rerun/{var['run_number']}/{var['Seq_File_without_low_SP_SEQ_with_Names']}\')}}\">
                        <br>'''
                else:
                    var['remove_seq_default_line'] = f'''
                        <A HREF="{var['Seq_File_without_low_SP_SEQ_with_Names']}" TARGET=_blank>The input sequences after the removal of unreliable sequences (with confidence score below {var['SP_SEQ_CUTOFF']})</A>
                        <font size=-1> (see list of removed sequences <A HREF="{var['removed_low_SP_SEQ_With_Names']}"
                        TARGET=_blank>here</A></font>)&nbsp;&nbsp;&nbsp;
                        <INPUT TYPE=\"BUTTON\" VALUE=\"run GUIDANCE on the confidently-aligned sequences only\" 
                        ONCLICK=\"window.open(\'{CONSTS.WEBSERVER_URL}/rerun/{var['run_number']}/{var['Seq_File_without_low_SP_SEQ_with_Names']}\')\">
                        <br>'''
                        
            var['List_Of_Alternative_MSAs'] = os.path.join( var['WorkingDir'], 'List_Of_Default_and_AltMSAs.txt')
            if form['PROGRAM'] == 'GUIDANCE2':
                var['SuperMSA_selection'] = print_SuperMSA_selection(var['List_Of_Alternative_MSAs'])
            
        elif state == State.Error or state == State.Crashed:
        
            var['error_type'] = error_type
            if error_type == 'system':
                var['error_msg'] = 'System Error'
            else:
                var['error_msg'] = error_msg
        
        # save vars
        vars_path = os.path.join( var['WorkingDir'], "VARS.json")
        try: 
            with open(vars_path, 'w') as fp:
                json.dump(var, fp)
            fp.close()
        except:
            error = f'GuidanceState({state.str()}): saving var failed'
            raise Exception ( error, "system")
        
        return warning_msg
        
    def validateInput(self):
    
        form = self.form
        var = self.var
        warning_msg = ''
        # KSENIA
        var['errors_file'] = f"results/{var['run_number']}/errors.txt"
        
        # validate seqs
        job_logger = logging.getLogger(var['run_number'])
        job_logger.info(f'\n{"#" * 80}\nvalidating seqs\n')
        
        if form['Redirect_From_MAFFT'] == '0': # regular run
        
            if not dict_file_defined_not_empty('userMSA_File', form): # Seq file provided
            
                if form['Seq_Type'] != 'Codons':
                    job_logger.info (f'validate_Seqs({var["WorkingDir"]},{var["SeqsFile"]},{form["Seq_Type"]}, False):\n')
                    ans = InputValidator.validate_Seqs( var['WorkingDir'], var['SeqsFile'] , form['Seq_Type'] , False)
                    var['LongestSeq'] = InputValidator.get_max_seq_length(os.path.join( var['WorkingDir'], var['SeqsFile']))
                else: 
                    job_logger.info (f'validate_Seqs({var["WorkingDir"]},{var["SeqsFile_Codons"]},{form["Seq_Type"]},False ):\n')
                    ans = InputValidator.validate_Seqs( var["WorkingDir"], var['SeqsFile_Codons'] , form['Seq_Type'] , False)
                    var['LongestSeq'] = InputValidator.get_max_seq_length(os.path.join( var["WorkingDir"], var['SeqsFile_Codons']))
                    
                job_logger.info(f'return: {join_list(ans)}\n')
                if ans[0] == "sys_error":
                    raise Exception ( ans[1], "system")
                elif ans[0] != 'OK':
                    #raise Exception ( join_list(ans), "user")
                    raise Exception ( ans, "user")
                var['SeqsFile'] = ans[2]
                var['NumOfSeq'] = ans[3]
                
                if ans[0] == 'OK' and ans[1]:
                    job_logger.warning(f'Warning: {ans[1]}. Nevertheless calculation is continued\n')
                    
                    warning_msg = f"<br><b><font color=\"red\" size=\"3\">Warning:</b></font><font size=\"3\"> {ans[1]}; The calculation continues.</font>\n"

                # uncommented this
                if form['Seq_Type'] == 'Codons':
                   shutil.move(os.path.join( var['WorkingDir'], var['SeqsFile']), os.path.join( var['WorkingDir'], var['SeqsFile_Codons']))
                    
            else:

                alignment_file = os.path.join( var['WorkingDir'], var['Alignment_File'])
                alignment_file_not_empty = False
                if os.path.exists(alignment_file): 
                    if os.path.getsize(alignment_file) > 0:
                        alignment_file_not_empty = True
                        
                if alignment_file_not_empty:
                    # Alignment provided
                    if form['Seq_Type'] != 'Codons':
                        job_logger.info (f"validate_Seqs({var['WorkingDir']},{var['Alignment_File']},{form['Seq_Type']},True):\n")
                        ans = InputValidator.validate_Seqs( var['WorkingDir'], var['Alignment_File'] , form['Seq_Type'] , True)
                        var['LongestSeq'] = InputValidator.get_max_seq_length(alignment_file)
                        
                    else: # Codon Alignment
                        job_logger.info (f"validate_Seqs({var['WorkingDir']},{var['Alignment_File']},{form['Seq_Type']}, True, {form['CodonTable']}):\n")
                        ans = InputValidator.validate_Seqs( var['WorkingDir'], var['Alignment_File'] , form['Seq_Type'] , True, form['CodonTable'])
                        var['LongestSeq'] = InputValidator.get_max_seq_length(alignment_file)
                        
                    job_logger.info(f'return: {join_list(ans)}\n')
                    if ans[0] == "sys_error":
                        raise Exception ( ans[1], "system")
                    elif ans[0] != 'OK':
                        #raise Exception (join_list(ans), "user")
                        raise Exception (ans, "user")
                    var['Alignment_File'] = ans[2]
                    var['NumOfSeq'] = ans[3]

                    if ans[0] == 'OK' and ans[1]:
                        job_logger.warning(f'Warning: {ans[1]}. Nevertheless calculation is continued\n')
                        
                        warning_msg = f"<br><b><font color=\"red\" size=\"3\">Warning:</b></font><font size=\"3\"> {ans[1]}; The calculation continues.</font>\n"

            if var['NumOfSeq'] < 4  and form['PROGRAM'] == "GUIDANCE":
                error = f"Only {var['NumOfSeq']} sequences were provided, however at least 4 sequences are requiered for GUIDANCE.<br>You can run HoT algorithm instead."
                raise Exception ( error, "user")
            
            if var['NumOfSeq'] < 3  and form['PROGRAM'] == "GUIDANCE2":
                error = f"Only {var['NumOfSeq']} sequences were provided, however at least 3 sequences are requiered for GUIDANCE2.<br>You can run HoT algorithm instead."
                raise Exception ( error, "user")
                
            if var['NumOfSeq'] >=300: 
                error = f"Due to limited computational resources, the web-server support analysis of up to 300 sequences. You can <a href=\"/source\"  target=\"_blank\">install and use the command-line version</a> or reduce the number of sequences and submit the job again.<br>Feel free to <a href=\"mailto:haimash\@tau.ac.il\" target=\"_blank\">contact us</a> for additional help. Sorry for the inconvenience."
                raise Exception ( error, "user")
            
            if var['LongestSeq'] > 6000 and form['Seq_Type'] == "Codons":
                error = f"Due to limited computational resources, the web-server support analysis of sequences no longer than 6,000 bp. You can <a href=\"/source\"  target=\"_blank\">install and use the command-line version</a> or reduce the size of sequences and submit the job again.<br>Feel free to <a href=\"mailto:haimash\@tau.ac.il\" target=\"_blank\">contact us</a> for additional help. Sorry for the inconvenience."
                raise Exception ( error, "user")
            
            if var['LongestSeq'] > 6000 and form['Seq_Type'] == "Nucleotides":
                error = f"Due to limited computational resources, the web-server support analysis of sequences no longer than 6,000 bp. You can <a href=\"/source\"  target=\"_blank\">install and  use the command-line version</a> or reduce the size of sequences and submit the job again.<br>Feel free to <a href=\"mailto:haimash\@tau.ac.il\" target=\"_blank\">contact us</a> for additional help. Sorry for the inconvenience."
                raise Exception ( error, "user")
                
            if var['LongestSeq'] > 2000 and form['Seq_Type'] == "AminoAcids":
                error = f"Due to limited computational resources, the web-server support analysis of sequences no longer than 2,000 AA. You can <a href=\"/source\"  target=\"_blank\">install and use the command-line version</a> or reduce the size of sequences and submit the job again.<br>Feel free to <a href=\"mailto:haimash\@tau.ac.il\" target=\"_blank\">contact us</a> for additional help. Sorry for the inconvenience."
                raise Exception ( error, "user")
        
        else:  # mafft case
        
            job_logger.info (f'validate_Seqs({var["WorkingDir"]},{var["Alignment_File"]},{form["Seq_Type"]}, True):')
            ans = InputValidator.validate_Seqs( var['WorkingDir'], var['Alignment_File'] , form['Seq_Type'] , True)
            
            job_logger.info(f'return: {join_list(ans)}\n')
            if ans[0] == "sys_error":
                raise Exception ( ans[1], "system")
            elif ans[0] != 'OK':
                #raise Exception ( join_list(ans), "user") JS
                raise Exception ( ans, "user")
            var['Alignment_File'] = ans[2]
            var['NumOfSeq'] = ans[3]
    
            if ans[0] == 'OK' and ans[1]:
                job_logger.warning(f'Warning: {ans[1]}. Nevertheless calculation is continued\n')
                warning_msg = f"<br><b><font color=\"red\" size=\"3\">Warning:</b></font><font size=\"3\"> {ans[1]}; The calculation continues.</font>\n"
                
        job_logger.info(f"ls of {var['WorkingDir']} fields:\n{os.listdir(var['WorkingDir'])}\n")
        
        return 'OK', warning_msg
        
    def submit_job(self):
    
        job_logger = logging.getLogger(self.var['run_number'])
        job_logger.info('going to store elements')
        status = self.save_state()
        job_logger.info(f'returned from save_state: status = {status}')
        var = self.var
        
        # Prepare and fqsub script
        parameters = f"{os.path.join(var['WorkingDir'],'input.data')} {os.path.join(var['WorkingDir'],'form.data')}"
        cmds_file = os.path.join(var['WorkingDir'], 'qsub.cmds')
        GuidanceState.write_cmds_file(cmds_file, parameters, var['WorkingDir'], var['run_number'])

        job_id_file = os.path.join(var['WorkingDir'], 'QSTAT_NO') 

        # a simple command when using shebang header (#!) in q_submitter_power.py
        #submission_cmd = f"{CONSTS.Q_SUBMITTER_SCRIPT} {cmds_file} {var['WorkingDir']} -q 'pupkowebr' --verbose > {job_id_file}"
        submission_cmd = f"{CONSTS.Q_SUBMITTER_SCRIPT} {cmds_file} {var['WorkingDir']} -q 'pupkowebr' --verbose"

        # subprocess call
        job_logger.info(f'\nSUBMITTING JOB TO QUEUE:\n{submission_cmd}\n')
        #subprocess.call(submission_cmd, shell=True)
        
        job_run_output = subprocess.run(submission_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        jobId = job_run_output.stdout.decode('utf-8').split('.')[0]
        
        job_logger.info(f'returned jobId = {jobId}\n')
        '''
        try: 
            with open (job_id_file, "r") as f:
                jobId = f.read()
                jobId.rstrip()
            f.close()
            
        except:
            
            raise Exception (f"submit job: openning {job_id_file} for reading failed", "sys")
        '''
        
        with open (job_id_file, 'w') as f:
            f.write(jobId)
        f.close()
        
        return jobId
    
    '''
    def send_email(self):
    
        form = self.form
        email = form['email_add'].strip()
        job_title = form['JOB_TITLE'].strip()
        run_number = self.var['run_number']
        
        notification_content = f'Your submission details are:\n{job_title}\n'

        notification_content += f'Once the analysis will be ready, we will let you know! Meanwhile, you can track the ' \
            f'progress of your job at:\n{CONSTS.WEBSERVER_URL}/results/{run_number}\n\n'

        send_email(smtp_server=CONSTS.SMTP_SERVER,
                   sender=CONSTS.ADMIN_EMAIL,
                   receiver=f'{email}',
                   subject=f'{CONSTS.WEBSERVER_NAME.upper()} - your job has been submitted! (Run number: {run_number})',
                   content=notification_content)
    '''
              
    def send_system_error_email(self):
    
        email = self.form['email_add'].strip()
        run_number = self.var['run_number']
        
        # get log file
        job_logger = logging.getLogger(run_number)
        handler = job_logger.handlers[0]
        logFile = os.path.basename(handler.baseFilename)
        
        if email != '': 
            send_email(smtp_server=CONSTS.SMTP_SERVER,
                       sender=CONSTS.ADMIN_EMAIL,
                       receiver=f'{CONSTS.ADMIN_EMAIL}',
                       subject=f'{CONSTS.WEBSERVER_NAME.upper()} job {run_number} by {email} has  failed!',
                       content=f"{email}\n\n{os.path.join(CONSTS.WEBSERVER_URL, 'results', run_number)}\n"
                       f"\n{os.path.join(CONSTS.WEBSERVER_URL, 'results', run_number, logFile)}")
        
        
    def get_process_id(self):
        
        return self.var['run_number']
    
    # to remove
    def handle_error(self, msg, error_type):
    
        email = form['email_add'].strip()
        run_number = self.var['run_number']
        
        # get log file
        job_logger = logging.getLogger(run_number)
        handler = job_logger.handlers[0]
        logFile = os.path.basename(handler.baseFilename)
        
        output_path = os.path.join(self.var['WorkingDir'], CONSTS.RESULT_WEBPAGE_NAME)
        
        job_logger.info(f'output_path = {output_path}')
        
        
        if email != '' and error_type == "system": 
            send_email(smtp_server=CONSTS.SMTP_SERVER,
                       sender=CONSTS.ADMIN_EMAIL,
                       receiver=f'{CONSTS.ADMIN_EMAIL}',
                       subject=f'{CONSTS.WEBSERVER_NAME.upper()} job {run_number} by {email} has  failed!',
                       content=f"{email}\n\n{os.path.join(CONSTS.WEBSERVER_URL, 'results', run_number)}\n"
                       f"\n{os.path.join(CONSTS.WEBSERVER_URL, 'results', run_number, logFile)}")

        job_logger.info('send_email')
        
        # Must be after flushing all previous data. Otherwise it might refresh during the writing.. :(
        sleep(2)
        
        job_logger.info('end of method')
    
    def log_job (self):
        
        form = self.form
        if 'email_add' in form.keys() and form['email_add'] != '':
            email = form['email_add'].strip()
        else:
            email = ''
        
        form = self.form
        if form['Redirect_From_MAFFT'] == '1':
            caller = 'FROM_MAFFT_SERVER'
        else:
            caller = 'GUIDANCE2'
        
        hostname = socket.gethostname()    
        user_ip = socket.gethostbyname(hostname)
    
        with open(CONSTS.SUBMISSIONS_LOG, 'a') as f:
            f.write(f"{currentTime()} {self.var['run_number']} {user_ip} {email}\t{caller}\n")
        f.close()
        
    # global methods
    
    def write_cmds_file(cmds_file, parameters, working_dir, run_number):
        # the queue does not like very long commands so I use a dummy delimiter (!@#) to break the commands for q_submitter
        new_line_delimiter = '!@#'

        required_modules_as_str = ' '.join(CONSTS.REQUIRED_MODULES)
        with open(cmds_file, 'w') as f:
            f.write(f'module load {required_modules_as_str};')
            f.write('bash activate /home/josefspr/bioseq/guidance/miniconda_guidance_env;')
            f.write(f'cd {working_dir};')
            f.write(new_line_delimiter)
            f.write(f'perl {CONSTS.MAIN_SCRIPT} {parameters} > {os.path.join(working_dir, "std.out")}\t{CONSTS.WEBSERVER_NAME}_{run_number}')
            
            f.write('\n')
        f.close()

    def get_state(process_id):
    
        job_logger = logging.getLogger(process_id)
        job_logger.info(f'entered method')
        
        running = False
        
        working_dir = os.path.join( CONSTS.WEBSERVER_RESULTS_DIR, process_id)
        print (f'get_state: working_dir = {working_dir}')
        if not os.path.exists(working_dir):
            return None

        var_path = os.path.join( CONSTS.WEBSERVER_RESULTS_DIR, process_id, 'VARS.json')
        with open (var_path, 'r') as f:
            var = json.load (f)
        f.close()
        
        if var['state'] == 'FINISHED':
            return State.Finished
        if var['state'] == 'ERROR':
            return State.Error
        elif var['state'] == 'FAILED': 
            return State.Crashed
        elif var['state'] == 'SUBMIT':
            return State.Running
        else:
            return State.Init
        
        
    def job_ended(process_id):
    
        end_page = os.path.join( CONSTS.WEBSERVER_RESULTS_DIR, process_id, f'GUIDANCE_{process_id}.END_OK')
        job_logger = logging.getLogger(process_id)
        job_logger.info(f'end_page = {end_page}')
        
        if os.path.exists (end_page):
            return True
        else:
            return False
        
    def get_results_html(process_id): 
    
        # check process_id exists
        working_dir = os.path.join( CONSTS.WEBSERVER_RESULTS_DIR, process_id)
        result_page = os.path.join(working_dir, CONSTS.RESULT_WEBPAGE_NAME)
        if not os.path.exists(result_page):
            return None
            
        results_url = os.path.join(CONSTS.WEBSERVER_RESULTS_URL, process_id)
        output_url = os.path.join(results_url, CONSTS.RESULT_WEBPAGE_NAME)
        
        return output_url
    
    
