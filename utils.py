from enum import Enum
import datetime
import logging
import os
import sys
import shutil
import SharedConsts as CONSTS
import pickle
import string
import re
import uuid
import json
from time import time
from random import randint


LOGGER_LEVEL_JOB_MANAGE_THREAD_SAFE = logging.DEBUG
LOGGER_LEVEL_JOB_MANAGE_API = logging.DEBUG
Bin = os.path.dirname(os.path.abspath(sys.argv[0]))
BIN_DIR = os.path.dirname(Bin)
SERVER_DIR = BIN_DIR

def init_dir_path():
    path2change = SERVER_DIR
    os.chdir(path2change)
    
init_dir_path()
logging_file_name = os.path.join('logs/', datetime.datetime.now().strftime('%Y_%m_%d_%H_%M.log'))
FORMAT = '%(asctime)s[%(levelname)s][%(filename)s][%(funcName)s]: %(message)s'
formatter = logging.Formatter('%(asctime)s[%(levelname)s][%(filename)s][%(funcName)s]: %(message)s')
logging.basicConfig(filename = logging_file_name, level=logging.WARNING, format=FORMAT)
logger = logging.getLogger('main')


class State(Enum):
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
            State.NotExists: 'NOT EXISTS'
        }[self]


def send_email_old(smtp_server, sender, receiver, subject='', content=''):
    from email.mime.text import MIMEText
    from smtplib import SMTP
    msg = MIMEText(content)
    msg['Subject'] = subject
    msg['From'] = sender
    msg['To'] = receiver
    s = SMTP(smtp_server)
    s.send_message(msg)
    s.quit()

def send_email ( smtp_server, sender, receiver, subject='', content=''):
    cmd = f'sendemail -f {sender} -t {receiver} -u \'{subject}\' -m \'{content}\' -s {CONSTS.SMTP_SERVER} -o tls=yes -xu {sender} -xp {CONSTS.ADMIN_PASSWORD}'
    os.system(cmd)

def currentTime(): 
    now = datetime.datetime.now()
    return now.strftime("%H:%M:%S %Y-%m-%d")

def remove_job_from_running_log(run_number):
    
    with open(CONSTS.PEPITOPE_RUNNING_JOBS, "r+") as f:
        d = f.readlines()
        f.seek(0)
        for i in d:
            if not run_number in i:
                f.write(i)
        f.truncate()
    f.close()

def write_daily_test(run_number, status):
    date = datetime.datetime.today().strftime('%d%m%Y')
    results_url = f'http://pepitope.tau.ac.il/results/{run_number}/output.html'
    with open(os.path.join(CONSTS.DAILY_TEST_DIR, f'pepitope_{date}.txt'), "w") as f:
        f.write(f'{status},{results_url}')
    f.close()
    
def save_obj( obj, path):
    with open(path, "wb") as f: 
        pickle.dump( obj, f)
    f.close()
    
def load_obj(path):
    with open(path, "rb") as f: 
        obj = pickle.load(f)
    f.close()
    return obj
    
def join_list (list):
    join_str = ''
    comma = ''
    for n in list: 
        join_str= join_str + comma + str(n)
        comma = ','
    return join_str

# used by Job_Manager
def dict_key_defined_not_empty(key, form):
    if key not in form.keys():
        return False
    return form[key] != ''

def dict_file_defined_not_empty(key, files):
    if key not in files:
        return False
    return files[key].filename != ''
    # return files[key] != ''

def dict_key_value(key, value, form):
    if key not in form.keys():
        return False
    return form[key] == value

def str_defined_and_not_equal_char(line, ch):
    if not line or line == '':
        return False
    if line[0] != ch: 
        return True
    else:
        return False

def checkForSpam (inFile, working_dir): # check the input for the possible of spammer
    try: 
        with open (inFile, "r") as f: 
            isSpam = 0
            for line in f: 
                if re.search('href\s*\=\s*', line):
                    isSpam = 1
                    break
                if re.search('http.*?\:\/\/', line):
                    isSpam = 1
                    break
        f.close()
        if isSpam: 
            with open (os.path.join( working_dir, inFile)) as f_out:
                pass
            f_out.close()
        return isSpam
        
    except:
        raise Exception (f"checkForSpam can\'t open file: '{inFile}'", "sys")
        
def save_file_to_disk(wd, logger, file, file_name_on_disk=None):

    logger.info('Saving FILE to disk')
 
    logger.info(f'Uploaded file name is:{file.filename}')

    data_path = os.path.join(f'{wd}/{file_name_on_disk}')
    file.save(data_path)
    
    '''
    data = form[form_field_name]
    logger.info(f'{uploaded_file_name} first 100 chars are: {data[:100]}\n')

    if file_name_on_disk is None:
        file_name_on_disk = uploaded_file_name
    data_path = os.path.join(f'{wd}/{file_name_on_disk}')

    with open(data_path, 'wb') as data_f:
        data_f.write(data)
    data_f.close()
    '''
    
    logger.info(f'Uploaded data was saved to {data_path} successfully\n')

    return data_path, file_name_on_disk

def save_text_to_disk(form, wd, logger, form_field_name, file_name_on_disk):

    logger.info(f'Saving TEXT to disk')

    data = form[form_field_name].rstrip()
    logger.info(f'first 100 chars are: {data[:100]}\n')

    data_path = os.path.join(f'{wd}/{file_name_on_disk}')
    with open(data_path, 'w') as data_f:
        data_f.write(data)
    data_f.close()
    
    logger.info(f'Uploaded data was saved to {data_path} successfully\n')

    return data_path
    
def peek_form(form, files, logger):

    #logger = logging.getLogger('cgi')
    # for debugging
    sorted_form_keys = sorted(form.keys())
    logger.info(f'These are the keys that the CGI received:\n{"; ".join(sorted_form_keys)}\n\n')
    logger.info('Form values are:\n')
    for key in sorted_form_keys:
        logger.info(f'{key} = {form[key]}')
    if files:
        if 'usrSeq_File' in files:
            file = files['usrSeq_File']
            if file and file.filename != '':
                logger.info(f'usrSeq_File = {file.filename}')
            '''
            if form[key].filename:
                logger.info(f'100 first characters of {key} = {form[key][:100]}\n')
            else:
                logger.info(f'{key} = {form[key]}\n')
            '''
        
def is_float(str):
    try:
        val = float(str)
        return True
    except ValueError:
        return False

def print_remove_site_selection_mask (sp_res_file): 

    return_str=''
    Cutoff_Removed_Pos_Hash = {}
    count_lines = 0
    with open (sp_res_file, 'r') as f: 
        for line in f:
            count_lines = count_lines + 1
            line = line.strip()
            if re.search( '^[0-9]+', line): 
                line_array = line.split()
                if is_float(line_array[2]):
                    score = "{:.3f}".format(float(line_array[2]))
                    #if score == 'nan':
                    #   continue
                    if score in Cutoff_Removed_Pos_Hash.keys(): 
                        Cutoff_Removed_Pos_Hash[score] = Cutoff_Removed_Pos_Hash[score] + 1
                    else:
                        Cutoff_Removed_Pos_Hash[score] = 1
        
    f.close()
    
    return_str = '<select name="cutoff" id="cutoff">'
    for cutoff in sorted(Cutoff_Removed_Pos_Hash.keys()): 
        res_bellow_Cutoof=0
        for score in Cutoff_Removed_Pos_Hash.keys():
            if score < cutoff:
                res_bellow_Cutoof += int(Cutoff_Removed_Pos_Hash[score])
        res_remain_precent = "{:.1f}".format( (1-(res_bellow_Cutoof/count_lines)) * 100) 
        if res_remain_precent != '100.0' and res_remain_precent != '0.0' and cutoff != 'nan':
            return_str += f'<option value="{cutoff}">{cutoff} ({res_remain_precent}% of residues remain)</option>'
    
    return return_str

def print_remove_site_selection_box (sp_res_file, msa_length): 

    return_str=''
    Cutoff_Removed_Pos_Hash = {}
    with open (sp_res_file, 'r') as f: 
        for line in f:
            line = line.strip()
            if re.search( '^[0-9]+', line): 
                line_array = line.split()
                if is_float(line_array[1]):
                    score = "{:.3f}".format(float(line_array[1]))
                    #if score == 'nan':
                    #   continue
                    if score in Cutoff_Removed_Pos_Hash.keys(): 
                        Cutoff_Removed_Pos_Hash[score] = Cutoff_Removed_Pos_Hash[score] + 1
                    else:
                        Cutoff_Removed_Pos_Hash[score] = 1
        
    f.close()
    # KSENIA
    return_str = '<select name="Col_Cutoff" id="Col_Cutoff">'
    for cutoff in sorted(Cutoff_Removed_Pos_Hash.keys()): 
        col_bellow_Cutoof=0
        for score in Cutoff_Removed_Pos_Hash.keys():
            if score < cutoff:
                col_bellow_Cutoof += int(Cutoff_Removed_Pos_Hash[score])
        col_remain_precent = "{:.1f}".format( (1-(col_bellow_Cutoof/msa_length)) * 100) 
        if col_remain_precent != '100.0' and col_remain_precent != '0.0' and cutoff != 'nan':
            return_str += f'<option value="{cutoff}">{cutoff} ({col_remain_precent}% of columns remain)</option>'
    
    return return_str

def print_remove_seq_selection_box (sp_res_file, msa_depth): 

    return_str=''
    Cutoff_Removed_Pos_Hash = {}
    with open (sp_res_file, 'r') as f: 
        for line in f:
            line = line.strip()
            if re.search( '^[0-9]+', line): 
                line_array = line.split()
                if is_float(line_array[1]):
                    score = "{:.3f}".format(float(line_array[1]))
                    if score in Cutoff_Removed_Pos_Hash.keys(): 
                        Cutoff_Removed_Pos_Hash[score] = Cutoff_Removed_Pos_Hash[score] + 1
                    else:
                        Cutoff_Removed_Pos_Hash[score] = 1
        
    f.close()

    return_str = '<select name="Seq_Cutoff" id="Seq_Cutoff">'
    for cutoff in sorted(Cutoff_Removed_Pos_Hash.keys()): 
        seq_bellow_Cutoof=0
        for score in Cutoff_Removed_Pos_Hash.keys():
            if score < cutoff:
                seq_bellow_Cutoof += int(Cutoff_Removed_Pos_Hash[score])
        seq_remain_precent = "{:.1f}".format( (1-(seq_bellow_Cutoof/msa_depth)) * 100) 
        if seq_remain_precent != '100.0' and seq_remain_precent != '0.0' and cutoff != 'nan':
            return_str += f'<option value="{cutoff}">{cutoff} ({seq_remain_precent}% of sequences remain)</option>'
    
    return return_str

def print_SuperMSA_selection(AltMSA_list): 

    Default_number_of_alternative_MSA_to_concat=20
    
    try: 
        with open( AltMSA_list, 'r') as f:
        
            return_str = "<select name=\"NumOfMSAs\" id=\"NumOfMSAs\">"
            AltMSA_to_Concat_including_baseMSA=0
            
            for line in f:
            
                AltMSA_to_Concat_including_baseMSA += 1
                AltMSA_to_Concat = AltMSA_to_Concat_including_baseMSA - 1
                if AltMSA_to_Concat > 0: 
                    if AltMSA_to_Concat == Default_number_of_alternative_MSA_to_concat:
                        return_str += f"<option selected=\"selected\" value=\"{AltMSA_to_Concat_including_baseMSA}\">{AltMSA_to_Concat}</option> "
                    else:
                        return_str += f"<option value=\"{AltMSA_to_Concat_including_baseMSA}\">{AltMSA_to_Concat}</option> "

            return_str += "</select>"
            
            f.close()
            return return_str
            
    except: 
        raise Exception (f"print_SuperMSA_selection: Can't open the list with AltMSA file: '{AltMSA_list}'", "sys")
        
def file_exists_not_empty(filename):
    if not os.path.exists(filename):
        return False
    else:
        if os.path.getsize(filename) > 0:
            return True
        else:
            return False

def convert_fs_to_lower_case (file):

    try:
        with open (file, 'r') as f:
            content = f.read()
        f.close()
    except: 
        return f'convert_fs_to_lower_case: Fail to open {file} for reading'
        
    try: 
        with open (file, 'r') as f:
            for line in f:
                if line[0] == '>':
                    f.write (line)
                else:
                    f.write(line.lower())
        f.close()
    except:
        f'convert_fs_to_lower_case: Fail to open {file} for writing'
        
    return 'OK'
    
def convert_fs_to_upper_case (file):

    try:
        with open (file, 'r') as f:
            content = f.read()
        f.close()
    except: 
        return f'convert_fs_to_upper_case: Fail to open {file} for reading'
        
    try: 
        with open (file, 'r') as f:
            for line in f:
                if line[0] == '>':
                    f.write(line)
                else:
                    f.write(line.upper())
        f.close()
    except:
        f'convert_fs_to_upper_case: Fail to open {file} for writing'
        
    return 'OK'
    
def get_job_logger(jobId):
    wd = os.path.join(CONSTS.WEBSERVER_RESULTS_DIR, jobId)
    if not os.path.exists(wd):
        return None
    logger = logging.getLogger(jobId)
    logger.setLevel(logging.INFO) # or whatever
    if len(logger.handlers) == 0:
        logFile = os.path.join(wd, f'{jobId}.log')
        handler = logging.FileHandler(logFile, 'w', 'utf-8') # or whatever
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.propagate = False
    return logger
    
def get_new_process_id():
    time_str = str(round(time()))
    rand_str = str(randint(1000,9999))
    return f'{time_str}{rand_str}'
        
