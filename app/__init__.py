#module load python/python-3.6.7

from flask import Flask, flash, request, redirect, url_for, render_template, Response, make_response, jsonify, send_file, send_from_directory
from werkzeug.middleware.proxy_fix import ProxyFix
from Job_Manager_API import Job_Manager_API
from GuidanceState import GuidanceState
from SharedConsts import UI_CONSTS
import SharedConsts as CONSTS
import logging
from utils import State, logger, get_job_logger, LOGGER_LEVEL_JOB_MANAGE_API, get_new_process_id #JS added
import os, sys
import warnings
from time import time, sleep
import re
import subprocess
# from flask_recaptcha import ReCaptcha
from google_recaptcha_flask import ReCaptcha
from GuidanceJobSubmitter import GuidanceJobSubmitter

# necessary in order to load the secret keys from .env
from dotenv import load_dotenv
load_dotenv()

#TODO think about it
warnings.filterwarnings("ignore")

USER_FILE_NAME = '' # not used by guidance
MAX_NUMBER_PROCESS = 10
TIME_OF_STREAMING_UPDATE_REQUEST_BEFORE_DELETING_IT_SEC = 1200

app = Flask(__name__, static_url_path='/guidance/static')
app.wsgi_app = ProxyFix(app.wsgi_app, x_for=1, x_proto=1, x_host=1, x_port=1)

# Ksenia
# app.config['APPLICATION_ROOT'] = '/'
app.config['APPLICATION_ROOT'] = '/guidance'
PREFIX = "/guidance"
# Ksenia

# all keys should be located in file .env due to security considerations
app.config['PREFERRED_URL_SCHEME'] = 'https'
app.config['SECRET_KEY'] = os.getenv("SECRET_KEY")
app.config['UPLOAD_FOLDERS_ROOT_PATH'] = CONSTS.WEBSERVER_RESULTS_DIR # path to folder
app.config['MAX_CONTENT_LENGTH'] = 16 * 1000 * 1000 * 1000 # MAX file size to upload
app.config['RECAPTCHA_SITE_KEY'] = os.getenv('RECAPTCHA_SITE_KEY')
app.config['RECAPTCHA_SECRET_KEY'] = os.getenv('RECAPTCHA_SECRET_KEY')
recaptcha = ReCaptcha(app) # Create a ReCaptcha object by passing in 'app' as parameter
process_id2update = []


@app.route(PREFIX + '/ConcatMSAs/<process_id>', methods=['GET', 'POST'])
def ConcatMSAs(process_id):

    if request.method == 'POST':
    
        working_dir = os.path.join( CONSTS.WEBSERVER_RESULTS_DIR, process_id)
        MSA_List = os.path.join(working_dir, 'List_Of_Default_and_AltMSAs.txt')
        NumOfMSAs = request.form['NumOfMSAs']
        Num_of_Alt = int(NumOfMSAs) - 1 # the base MSA is already included in the NumOfMSAs by the form
        OutPath = os.path.join( working_dir, f'SuperMSA_DefaultMSA_and_{Num_of_Alt}_Alt.fas')
        OutHTML = os.path.join(working_dir, CONSTS.RESULT_WEBPAGE_NAME)

        cmd = f"python3 {CONSTS.CONCAT_SCRIPT} {MSA_List} {OutPath} {str(NumOfMSAs)} NO YES {OutHTML}"
        try:
            subprocess.run(cmd, shell=True, timeout=300)
        except subprocess.TimeoutExpired:
            raise RuntimeError(f"Command timed out after 300s: {cmd}")
        return redirect(url_for('results', _anchor='remove_seq', process_id = process_id))

@app.route(PREFIX + '/remove_seq/<process_id>', methods=['GET', 'POST'])
def remove_seq(process_id):

    if request.method == 'POST':

        cmd = f"python3 {CONSTS.REMOVE_SEQ_SCRIPT} {request.form['VARS_json']} {request.form['FORM_json']} {request.form['Seq_Cutoff']}"
        try:
            subprocess.run(cmd, shell=True, timeout=300)
        except subprocess.TimeoutExpired:
            raise RuntimeError(f"Command timed out after 300s: {cmd}")
        return redirect(url_for('results', _anchor='remove_seq', process_id = process_id))

@app.route(PREFIX + '/remove_pos/<process_id>', methods=['GET', 'POST'])
def remove_pos(process_id):

    if request.method == 'POST':

        cmd = f"python3 {CONSTS.REMOVE_POS_SCRIPT} {request.form['VARS_json']} {request.form['Col_Cutoff']}"
        try:
            subprocess.run(cmd, shell=True, timeout=300)
        except subprocess.TimeoutExpired:
            raise RuntimeError(f"Command timed out after 300s: {cmd}")
        return redirect(url_for('results', _anchor='remove_pos', process_id = process_id))

@app.route(PREFIX + '/mask/<process_id>', methods=['GET', 'POST'])
def mask(process_id):

    if request.method == 'POST':

        cmd = f"python3 {CONSTS.MASK_SCRIPT} {request.form['VARS_json']} {request.form['type_a']} {request.form['cutoff']}"
        try:
            subprocess.run(cmd, shell=True, timeout=300)
        except subprocess.TimeoutExpired:
            raise RuntimeError(f"Command timed out after 300s: {cmd}")
        return redirect(url_for('results', _anchor='masks', process_id = process_id))
        
@app.route(PREFIX + '/update_html_ext/<process_id>')
def update_html_ext(process_id):

    job_logger = get_job_logger(process_id)
    log_msg = f'process_id = {process_id}'
    logger.info(log_msg)
    job_logger.info(log_msg)
    
    '''
    if process_id:
        process_id2update.append(process_id)
    '''
    return jsonify('data')
    
def update_html(process_id, state):

    job_logger = get_job_logger(process_id)
    log_msg = f'process_id = {process_id} state = {state}'
    logger.info(log_msg)
    job_logger.info(log_msg)
    
    if process_id:
        process_id2update.append(process_id)
    
@app.route(PREFIX + "/process_page_update/<process_id>")
def update_process_page(process_id):
    if process_id in process_id2update:
        process_id2update.remove(process_id)
        return 'update' #UI_CONSTS.TEXT_TO_RELOAD_HTML

#JS commented out
manager = Job_Manager_API(MAX_NUMBER_PROCESS, app.config['UPLOAD_FOLDERS_ROOT_PATH'], USER_FILE_NAME, update_html)
       
@app.route(PREFIX + '/process_state/<process_id>')
def process_state(process_id):

    job_state = manager.get_guidance_job_state(process_id)
    if job_state:
        job_logger = get_job_logger(process_id)
    else: 
        job_logger = None
    log_msg = f'process_id = {process_id}, job_state = {job_state}'
    logger.info (log_msg)
    if job_logger:
        job_logger.info (log_msg)
    
    if job_state == None:
        return redirect(url_for('error', error_type=UI_CONSTS.UI_Errors.UNKNOWN_PROCESS_ID.name))
        
    guidance_state = None
    try:
        guidance_state = GuidanceState(jobId = process_id)
    except:
        log_msg = f'GuidanceState.__init__({process_id}) failed'
        logger.info (log_msg)
        if job_logger:
            job_logger.info(log_msg)
        kwargs = {
            "var": guidance_state.var if guidance_state else {},
        }
        return render_template('error_page.html', error_text=CONSTS.SYS_ERROR_MSG, **kwargs)
    
    # get runtime warnings
    runtime_warnings_path = os.path.join (CONSTS.WEBSERVER_RESULTS_DIR, process_id, 'Warnings.txt')
    runtime_warnings = ''
    if os.path.exists( runtime_warnings_path):
        with open (runtime_warnings_path, 'r') as f:
            runtime_warnings = f.read()
        f.close()
    else:
        runtime_warnings = ''
        
    # update progress report
    progress_report_html = os.path.join (CONSTS.WEBSERVER_RESULTS_DIR, process_id, 'ProgressReport.html')
    progress_report = ''
    if os.path.exists( progress_report_html):
        with open (progress_report_html, 'r') as f:
            progress_report = f.read()
        f.close()
    else:
        progress_report = ''
    MSA_status_file = os.path.join (CONSTS.WEBSERVER_RESULTS_DIR, process_id, 'MSA_STATUS.txt')
    if os.path.exists( MSA_status_file):
        with open (MSA_status_file, 'r') as f:
            MSA_status = f.read()
        f.close()
    else:
        MSA_status = ''
    if progress_report:
        progress_report = progress_report.replace('REPLACE', MSA_status)
    
    # update state only if crashed or finished
    if job_state == State.Crashed:
        errors_file = os.path.join(CONSTS.WEBSERVER_RESULTS_DIR, process_id, 'errors.txt')
        if os.path.exists(errors_file):
            with open(errors_file) as _ef:
                error_content = _ef.read().strip()
            guidance_state.update_state(state=State.Crashed, error_msg=error_content, error_type='user')
        else:
            guidance_state.update_state(state=job_state, error_msg='System error', error_type='system')
    elif job_state == State.Finished:
        
        # check if GuidanceState should be updated
        if guidance_state.var['state'] != "FINISHED":
            guidance_state.update_state(state = State.Finished)
            sleep(5)
        
    if job_state != State.Finished and job_state != State.Crashed and job_state != State.Error:
    
        # job is still running, state is not updated
        
        kwargs = {
            "reload_interval": CONSTS.RELOAD_INTERVAL,
            "var": guidance_state.var, 
            "form": guidance_state.form,
            "job_state": job_state.str(),
            "runtime_warnings": runtime_warnings,
            "progress_report": progress_report,
            "mask_residues_list": [], 
            "remove_columns_list": [], 
            "remove_sequences_list": [],
            "super_msa_list": []
        }
        return render_template('running.html', **kwargs)
        
    else:
    
        # redirect to results
        return redirect(url_for('results', process_id = process_id))

@app.route(PREFIX + '/results/<int:process_id>', methods=['GET', 'POST'])
def results(process_id):
    
    process_id = str(process_id) 
    job_state = GuidanceState.get_state(process_id)
    if job_state:
        job_logger = get_job_logger(process_id)
    else: 
        job_logger = None
    log_msg = f'process_id = {process_id}, job_state = {job_state}'
    logger.info (log_msg)
    if job_logger:
        job_logger.info (log_msg)
    
    if job_state == None:
        guidance_state = GuidanceState(jobId=process_id)
        kwargs = {
            "var": guidance_state.var,
        }
        return render_template('error_page.html', error_text=f"Job does not exist", **kwargs)
    
    guidance_state = None
    try:
        guidance_state = GuidanceState(jobId = process_id)
    except:
        log_msg = f'GuidanceState.__init__({process_id}) failed'
        logger.info (log_msg)
        if job_logger:
            job_logger.info (log_msg)

        kwargs = {
            "var": guidance_state.var if guidance_state else {},
        }

        return render_template('error_page.html', error_text=CONSTS.SYS_ERROR_MSG, **kwargs)

    # if running check if job ended
    if job_state == State.Running:
    
        if GuidanceState.job_ended(process_id):
            # errors.txt is written by exit_on_error() when the job aborts on a user error.
            # If it exists the process exited via error, not normal completion.
            errors_file = os.path.join(CONSTS.WEBSERVER_RESULTS_DIR, process_id, 'errors.txt')
            if os.path.exists(errors_file):
                with open(errors_file) as _ef:
                    error_content = _ef.read().strip()
                guidance_state.update_state(state=State.Crashed, error_msg=error_content, error_type='user')
                job_state = State.Crashed
            else:
                guidance_state.update_state(State.Finished)
                job_state = State.Finished
        else:
            job_state_man = manager.get_guidance_job_state(process_id)
            if job_state_man == None or job_state_man == State.Crashed:
                # Check if errors.txt was written by exit_on_error() for a user error
                errors_file = os.path.join(CONSTS.WEBSERVER_RESULTS_DIR, process_id, 'errors.txt')
                if os.path.exists(errors_file):
                    with open(errors_file) as _ef:
                        error_content = _ef.read().strip()
                    guidance_state.update_state(state=State.Crashed, error_msg=error_content, error_type='user')
                else:
                    guidance_state.update_state(State.Crashed)
                job_state = State.Crashed
            
    # if jobs did not finish or is not an error redirect to process_state
    if job_state != State.Finished and job_state != State.Crashed and job_state != State.Error:
        return redirect(url_for('process_state', process_id = process_id))
    
    # get runtime warnings
    runtime_warnings_path = os.path.join (CONSTS.WEBSERVER_RESULTS_DIR, process_id, 'Warnings.txt')
    runtime_warnings = ''
    if os.path.exists( runtime_warnings_path):
        with open (runtime_warnings_path, 'r') as f:
            runtime_warnings = f.read()
        f.close()
    else:
        runtime_warnings = ''
        
    # update progress report
    progress_report_html = os.path.join (CONSTS.WEBSERVER_RESULTS_DIR, process_id, 'ProgressReport.html')
    progress_report = ''
    if os.path.exists( progress_report_html):
        with open (progress_report_html, 'r') as f:
            progress_report = f.read()
        f.close()
    else:
        progress_report = ''
    MSA_status_file = os.path.join (CONSTS.WEBSERVER_RESULTS_DIR, process_id, 'MSA_STATUS.txt')
    if os.path.exists( MSA_status_file):
        with open (MSA_status_file, 'r') as f:
            MSA_status = f.read()
        f.close()
    else:
        MSA_status = ''
    if progress_report:
        progress_report = progress_report.replace('REPLACE', MSA_status)
    
    # mask specific residues
    mask_residues_list = [] 
    for file in os.listdir(os.path.join (CONSTS.WEBSERVER_RESULTS_DIR, process_id)):
        m = re.search('Mask_Residues_Res_([0-9]*\.[0-9]+).aln',file)
        if m:
            if m.group(1) not in mask_residues_list:
                mask_residues_list.append(m.group(1))
                
    # remove unreliable columns
    remove_columns_list = [] 
    for file in os.listdir(os.path.join (CONSTS.WEBSERVER_RESULTS_DIR, process_id)):
        m = re.search('Without_low_SP_Col.([0-9]*\.[0-9]+)',file)
        if m:
            if m.group(1) not in remove_columns_list:
                remove_columns_list.append(m.group(1))
                
    # remove unreliable sequences
    remove_sequences_list = [] 
    for file in os.listdir(os.path.join (CONSTS.WEBSERVER_RESULTS_DIR, process_id)):
        m = re.search('Without_low_SP_Seq.([0-9]*\.[0-9]+)',file)
        if m:
            if m.group(1) not in remove_sequences_list:
                remove_sequences_list.append(m.group(1))
    
    # superMSA's
    super_msa_list = []
    for file in os.listdir(os.path.join (CONSTS.WEBSERVER_RESULTS_DIR, process_id)):
        m = re.search('SuperMSA_DefaultMSA_and_([0-9]+)_Alt.fas', file)
        if m:
            if m.group(1) not in super_msa_list:
                super_msa_list.append(m.group(1))

    # best MSA selected by the pretrained DL model
    best_msa_file = ''
    dataset = guidance_state.var.get('dataset', 'MSA')
    msa_program = guidance_state.form.get('MSA_Program', 'MAFFT')
    candidate_best_msa = f"{dataset}.{msa_program}.Guidance2_BestMSA.fasta"
    if os.path.exists(os.path.join(CONSTS.WEBSERVER_RESULTS_DIR, process_id, candidate_best_msa)):
        best_msa_file = candidate_best_msa

    sleep(3)

    errors_file = os.path.join(CONSTS.WEBSERVER_RESULTS_DIR, process_id, 'errors.txt')
    kwargs = {
        "reload_interval": CONSTS.RELOAD_INTERVAL,
        "var": guidance_state.var,
        "form": guidance_state.form,
        "job_state": job_state.str(),
        "runtime_warnings": runtime_warnings,
        "progress_report": progress_report,
        "mask_residues_list": sorted( mask_residues_list, reverse=True),
        "remove_columns_list": sorted( remove_columns_list, reverse=True),
        "remove_sequences_list": sorted( remove_sequences_list, reverse=True),
        "super_msa_list": super_msa_list,
        "best_msa_file": best_msa_file,
        "errors_file_exists": os.path.exists(errors_file),
        "admin_email": CONSTS.ADMIN_EMAIL,
    }

    #
    return render_template('output.html', **kwargs) #TODO - testing removing cache

    # response = make_response(render_template('output.html', **kwargs))
    # response.headers['Cache-Control'] = 'no-store, no-cache, must-revalidate, max-age=0'
    # response.headers['Pragma'] = 'no-cache'
    # response.headers['Expires'] = '0'
    # return response
            
@app.route(PREFIX + '/results/<int:_process_id>/<string:show_file>', methods=['GET', 'POST'])
def show_file( _process_id, show_file):

    process_id = str(_process_id)
    if '.html' in show_file:
        # # ksenia's code
        path = os.path.join(CONSTS.WEBSERVER_RESULTS_DIR, process_id, show_file)
        resp = make_response(open(path, 'r', encoding="ISO-8859-2").read())
        resp.headers["Content-type"] = "text/html;charset=UTF-8"
        return resp

        # return render_template(f'results/{process_id}/{show_file}')
    elif '.tar.gz' in show_file or '.zip' in show_file:
        path = os.path.join (CONSTS.WEBSERVER_RESULTS_DIR, process_id, show_file)
        return send_file(path, mimetype='application/octet-stream')
    elif '.php' in show_file:
        old_results_path = os.path.join( CONSTS.WEBSERVER_RESULTS_OLD_URL, process_id, show_file)
        return redirect (old_results_path)
    else:
        path = os.path.join(CONSTS.WEBSERVER_RESULTS_DIR, process_id, show_file)
        if not os.path.exists(path):
            return make_response(f"File not found: {show_file}", 404)
        resp = make_response(open(path, 'r', encoding="ISO-8859-1").read())
        resp.headers["Content-type"] = "text/plain;charset=UTF-8"
        return resp

@app.route(PREFIX + '/logs/<string:show_file>', methods=['GET', 'POST'])

def show_log_file(show_file):
    path = os.path.join (CONSTS.WEBSERVER_LOGS_DIR, show_file)
    resp = make_response(open(path, 'r', encoding="ISO-8859-1").read())
    resp.headers["Content-type"]="text/plain;charset=UTF-8"
    return resp

@app.route(PREFIX + '/testpost', methods=['GET', 'POST'])
def testpost():
    if request.method == 'POST':
    	return render_template('posted.html')
    else:
    	return render_template('testpost.html')
    	
@app.route(PREFIX + '/error/<error_type>')
def error(error_type):
    # checking if error_type exists in error enum
    # KSENIa REPLACE var

    kwargs = {
        "var": error_type,
    }
    try:
        return render_template('error_page.html', error_text=UI_CONSTS.UI_Errors[error_type].value, **kwargs)
    except:
        return render_template('error_page.html', error_text=f'Unknown error, \"{error_type}\" is not a valid error code', **kwargs)

@app.route(PREFIX + '/', methods=['GET', 'POST'])
def home():
    if request.method == 'POST':
    	
        new_process_id = get_new_process_id() #manager.get_new_process_id() # JS commented out
        job_logger = None
        logger.info (f'new_process_id = {new_process_id}')
        
        # add to daily runs
        #daily_page = os.path.join( CONSTS.WEBSERVER_LOGS_DIR, 'daily_runs', new_process_id)
        #open(daily_page, 'a').close()
        
        # create working directory
        #wd = os.path.join(app.config['UPLOAD_FOLDERS_ROOT_PATH'], new_process_id)
        #os.mkdir(wd)
        
        # guidance state part
        try:
        
            warning_messages = ''
            guidance_state = None
            
            if 'no_captcha' not in request.form.keys(): 
                if not recaptcha.verify():
                    raise Exception ("Run failed to pass I'm not a robot test", "user")
                
            # create Job Manager instance
            guidance_state = GuidanceState(jobId = new_process_id, form = request.form, files = request.files)
            job_logger = get_job_logger(new_process_id)
                
            # upload files
            guidance_state.upload_files()
            msg = guidance_state.update_state(state = State.Init)
            if msg != '':
                warning_messages = warning_messages + msg
            
            # validate input
            status, msg = guidance_state.validateInput()
            job_logger.info(f'validated data, status = {status}')
            # KSENIA added **kwargs
            if status != 'OK':
                kwargs = {
                    "var": guidance_state.var,
                }
                return render_template('error_page.html', error_text=msg, **kwargs) # change
            if msg != '':
                warning_messages = warning_messages + msg
            job_logger.info(f'validated data')
            
            # store elements
            job_logger.info(f'storing state')
            if 'crash' in request.form['JOB_TITLE']:
                crash_flag = True
            else:
                crash_flag = False
            status = guidance_state.save_state(warning_messages, crash_flag=crash_flag)
            if status != 'OK':
                kwargs = {
                    "var": guidance_state.var,
                }
                render_template('error_page.html', error_text=msg, **kwargs)  # change
                raise Exception(status, "system")
            
            #return render_template('posted.html', msg = 'stored')
            
        except Exception as e:
        
            if len(e.args) > 0:
                logger.info(e.args[0])
                if job_logger: 
                    job_logger.info(e.args[0])
                
            # html FAIL message and footer
            if len(e.args) > 1 and e.args[1]=="user":
                msg = e.args[0]
                error_type = "user"
            else:
                msg = CONSTS.SYS_ERROR_MSG
                error_type = "system"
            
            if len(e.args) > 1:
                log_error_msg = e.args[0]
            else:
                log_error_msg = 'no error message'
            
            log_msg = f'Guidance State caught exception: msg = {log_error_msg}, error_type = {error_type}'
            logger.info(log_msg)
            if job_logger: 
                job_logger.info(log_msg)
            
            if guidance_state: 
                guidance_state.update_state(state = State.Error, error_msg = msg, error_type = error_type)
                if error_type == "system":
                    guidance_state.send_system_error_email()
                    
            fail_page = os.path.join(CONSTS.WEBSERVER_RESULTS_DIR, new_process_id, f'GUIDANCE_{new_process_id}.END_FAIL')
            open(fail_page, 'a').close()

            kwargs = {
                "var": guidance_state.var,
            }
            
            return render_template('error_page.html', error_text=msg, **kwargs)
                
        # here the job is submitted
        
        email_address =  request.form['email_add']
        
        #add to users log file
        guidance_state.log_job()
        
        #return render_template('posted.html', msg = 'log job')
        
        if 'norun' not in request.form['JOB_TITLE']:
            #man_results = manager.add_guidance_process(new_process_id, email_address, request.form['JOB_TITLE'])
            working_dir = os.path.join ( CONSTS.WEBSERVER_RESULTS_DIR, new_process_id)
            returnVal = GuidanceJobSubmitter.submit_job(working_dir, email_address)
            if returnVal == 0:
            	return            
            man_results = True
        else:
            man_results = False
        
        #return render_template('posted.html', msg = 'added process')
        
        if not man_results:
            msg = f'job_manager_api can\'t add process {new_process_id}'
            logger.error(msg)
            job_logger.info(msg)
            if 'norun' not in request.form['JOB_TITLE']:
                kwargs = {
                    "var": guidance_state.var,
                }
                return render_template('error_page.html', error_text=CONSTS.SYS_ERROR_MSG, **kwargs)
            else:
                kwargs = {
                    "var": guidance_state.var,
                }
                return render_template('error_page.html', error_text=new_process_id, **kwargs)
        log_msg = f'process added man_results = {man_results}, redirecting url'
        logger.info(log_msg)
        job_logger.info(log_msg)
        return redirect(url_for('process_state', process_id=new_process_id))
        
    else:
        return render_template('home.html', FASTA_txt='') # JS daily_test = 'yes'

@app.route(PREFIX + '/rerun/<process_id>/<seqFile>', methods=['GET'])
def rerun( process_id, seqFile):

    job_logger = get_job_logger(process_id)
    working_dir = os.path.join(CONSTS.WEBSERVER_RESULTS_DIR, process_id)
    if not os.path.exists(working_dir):
        guidance_state = GuidanceState(jobId=process_id)
        kwargs = {
            "var": guidance_state.var,
        }
        return render_template('error_page.html', error_text=f'Process ID does not exist', **kwargs)
        
    try: 
        rerun_params = {}
        # get contents of seqFile
        FASTA_txt = ''
        with open(os.path.join(working_dir, seqFile), 'r') as f:
            for line in f:
                FASTA_txt = FASTA_txt + line.rstrip() + '\n'
        f.close()
        
        # get parameters of run
        with open(os.path.join(working_dir, 'rerun_param'), 'r') as f:
            for line in f:
                m = re.search( '\$([A-Za-z_]+)=([0-9\.]+);', line)
                if m:
                    rerun_params[m.group(1)] = m.group(2)
        f.close()
        
        log_msg = f'process_id = {process_id}, seqFile = {seqFile}, FASTA_txt = {FASTA_txt}, rerun_params = {rerun_params}'
        logger.info(log_msg)
        job_logger.info(log_msg)
        
        return render_template('home.html', FASTA_txt=FASTA_txt, rerun_params=rerun_params) #daily_test = 'yes'
    
    except Exception as e:
    
        log_msg = f'caught exception: process_id = {process_id}, seqFile = {seqFile}, FASTA_txt = {FASTA_txt}'
        logger.info(log_msg)
        job_logger.info(log_msg)
        guidance_state = GuidanceState(jobId=process_id)
        kwargs = {
            "var": guidance_state.var,
        }
        return render_template('error_page.html', error_text='Rerun failed', **kwargs)

#@app.route(PREFIX + '/run_from_mafft', methods=['GET'])
@app.route(PREFIX + '/index_FromMAFFT.php', methods=['GET'])
def run_from_mafft():

    import requests
    
    args = request.args
    calling_run = args.get("run")
    calling_args= args.get("args")
    
    url= f'http://mafft.cbrc.jp/alignment/server/spool/{calling_run}.pir'
    r = requests.get(url, allow_redirects=True)
    
    return render_template('home.html', FASTA_txt = r.content.decode('ascii'), back_from_mafft = 1, calling_run = calling_run, calling_args = calling_args)
    #return render_template('home.html', FASTA_txt = 'AAAA')

@app.route(PREFIX + '/index_rerun_same_seq.php', methods=['GET'])
def rerun_same_seq_php():
    import requests
    
    args = request.args
    calling_run =  args.get('run')
    seqFile =  args.get('file')
    
    path = os.path.join( '/rerun', calling_run, seqFile)
    return redirect (path)
    
@app.route(PREFIX + '/index_rerun.php', methods=['GET'])
def rerun_php():
    import requests
    
    args = request.args
    calling_run =  args.get('run')
    seqFile =  args.get('file')
    
    path = os.path.join( '/rerun', calling_run, seqFile)
    return redirect (path)
        
@app.route(PREFIX + '/overview', methods=['GET'])
def overview():
    return render_template('overview.html')

@app.route(PREFIX + '/gallery', methods=['GET'])
def gallery():
    return render_template('gallery.html')

@app.route(PREFIX + '/source', methods=['GET'])
def source():
    return render_template('source.html')
    
@app.route(PREFIX + '/credits', methods=['GET'])
def credits():
    return render_template('credits.html')
    
@app.errorhandler(404)
def page_not_found(e):
    return redirect(url_for('error', error_type=UI_CONSTS.UI_Errors.PAGE_NOT_FOUND.name))

@app.route(PREFIX + "/about")
def about():
    return render_template('about.html')
    
@app.route(PREFIX + '/run/<int:index>', methods=['GET', 'POST'])
def run(index):
    
    logger.info('calling run')
            
    forms = []
    form = {"MSA_Program": "MAFFT", "MSA_Program_GUIDANCE" : "MAFFT", "JOB_TITLE": "test", "Bootstraps": "100", "email_add": "josefspr@gmail.com", "PROGRAM": "GUIDANCE3", "Seq_Type": "AminoAcids", "FASTA_txt": ">NC_001802.HXB2\r\nMQPIPIVAIVALVVAIIIAIVVWSIVIIEYRKILRQRKIDRLIDRLIERAEDSGNESEGEISALVEMGVEMGHHAPWDVDDL\r\n>EF637049.B\r\nMQSLQIVAIVALVVTAIIAIVVWSIVLIEYRKLLRQRKIDRLIDRIRERAEDSGNESEGDQEELAGLVERGHLAPWDVDDL\r\n>EF514700.B\r\nMQPLEILAIVALVVAIILAIVVWTIVFIEYKKILRQRKIDRLIDRIAERAEDSGNESEGDQEELSALVDMGHDAPWVVVDQ\r\n>DQ056417.C\r\nMLESIDYRLGVAALLLALIIAIIVWIIAYLEYRKLLRQRRIDKLIKRIRERAEDSGNESEGDIEELSTMVDVEHLRLLDVNNL\r\n>AY463217.C\r\nMVDLLAGVDYRVGVGALIIALIIAIIVWIWVYIEYRKLLRQRKIDWLIKRLREREEDSGNESEGDTEELATMVDMGHLRLLDDNNV\r\n>DQ011165.C\r\nMLNFLAGVDYRIGVGALIVGLIIAIVVWIIVYLEYRKLVKQRKIDWLIERIRERAEDSGNESEGDTEELATMVDMGHLRLLDAYDL\r\n>AB254142.C\r\nMINFAARVDYRVGVAAFTIALIIAIVVWIIVYLELVRQRKIDQLIIRIREREEDSGNESEGDIEELSTMVDMGQLRLLDGNGL\r\n>AY901969.C\r\nMVNLLEKVNLFEKVDYRLGVGALLIALVIAIIVWTIAYIEYRKLVRQRKIDWLVKRIRERAEDSGNESDGDTEELSTMVDLGHLRLLDVAEL\r\n>EU110088.A1\r\nMNQLQILAIXGLVVALILAIVVWTIVGIEYRKLLRQRRIDRLIKRISERAEDSGNESDGDTEELSQLVEMGNYNLGFDDNL\r\n>AB253428.A1\r\nMQLLEICAVVGLVVALIIAIVVWTIVGIEYKKLLKQRKIDRLVDRIRERAEDSGNESDGDREELSLLVDMGDYDLGDDNNL\r\n>AF457052.A1\r\nMLSALEICAIAGLVIALIIAIVVWTIVGIEYRRLLKQRKIDRLIERIRERAEDSGNESDGDTEELAALIEMGNYDLGDANDL\r\n>AF077336.F1\r\nMSYLLAIGIAALIVALIIAIVVWTIVYIEYKKLVRQRKINKLYKRIRERAEDSGNESEGDAEELAALGEMGPFIPGDINNL\r\n>DQ168575.G\r\nMKSLEISAIVGLIVAFIAAIVVWTIVLIEYRKIRKQKRIDKILDRIRERAEDSGNESEGDTEELATLVDMVDFEPWVGDNL\r\n>AY795907.D\r\nMQTLEILSIVALVIAAIIAIIVWTIVYIEYRKIRRQRKIDQLIDRIRERAEDSGNESEGDEEELSTLMEMGHAAPWNVADDL\r\n", "SP_COL_CUTOFF" : "0.93", "SP_SEQ_CUTOFF": "0.6", "maxiterate":"0", "outorder": "aligned", "GENCODE": "1", "F": "+F", "solved" : "false", "pair" : ""}
    forms.append(form)
    form = {"MSA_Program": "MAFFT", "MSA_Program_GUIDANCE" : "MAFFT", "JOB_TITLE": "test", "Bootstraps": "100", "email_add": "josefspr@gmail.com", "PROGRAM": "GUIDANCE", "Seq_Type": "AminoAcids", "FASTA_txt": ">NC_001802.HXB2\r\nMQPIPIVAIVALVVAIIIAIVVWSIVIIEYRKILRQRKIDRLIDRLIERAEDSGNESEGEISALVEMGVEMGHHAPWDVDDL\r\n>EF637049.B\r\nMQSLQIVAIVALVVTAIIAIVVWSIVLIEYRKLLRQRKIDRLIDRIRERAEDSGNESEGDQEELAGLVERGHLAPWDVDDL\r\n>EF514700.B\r\nMQPLEILAIVALVVAIILAIVVWTIVFIEYKKILRQRKIDRLIDRIAERAEDSGNESEGDQEELSALVDMGHDAPWVVVDQ\r\n>DQ056417.C\r\nMLESIDYRLGVAALLLALIIAIIVWIIAYLEYRKLLRQRRIDKLIKRIRERAEDSGNESEGDIEELSTMVDVEHLRLLDVNNL\r\n>AY463217.C\r\nMVDLLAGVDYRVGVGALIIALIIAIIVWIWVYIEYRKLLRQRKIDWLIKRLREREEDSGNESEGDTEELATMVDMGHLRLLDDNNV\r\n>DQ011165.C\r\nMLNFLAGVDYRIGVGALIVGLIIAIVVWIIVYLEYRKLVKQRKIDWLIERIRERAEDSGNESEGDTEELATMVDMGHLRLLDAYDL\r\n>AB254142.C\r\nMINFAARVDYRVGVAAFTIALIIAIVVWIIVYLELVRQRKIDQLIIRIREREEDSGNESEGDIEELSTMVDMGQLRLLDGNGL\r\n>AY901969.C\r\nMVNLLEKVNLFEKVDYRLGVGALLIALVIAIIVWTIAYIEYRKLVRQRKIDWLVKRIRERAEDSGNESDGDTEELSTMVDLGHLRLLDVAEL\r\n>EU110088.A1\r\nMNQLQILAIXGLVVALILAIVVWTIVGIEYRKLLRQRRIDRLIKRISERAEDSGNESDGDTEELSQLVEMGNYNLGFDDNL\r\n>AB253428.A1\r\nMQLLEICAVVGLVVALIIAIVVWTIVGIEYKKLLKQRKIDRLVDRIRERAEDSGNESDGDREELSLLVDMGDYDLGDDNNL\r\n>AF457052.A1\r\nMLSALEICAIAGLVIALIIAIVVWTIVGIEYRRLLKQRKIDRLIERIRERAEDSGNESDGDTEELAALIEMGNYDLGDANDL\r\n>AF077336.F1\r\nMSYLLAIGIAALIVALIIAIVVWTIVYIEYKKLVRQRKINKLYKRIRERAEDSGNESEGDAEELAALGEMGPFIPGDINNL\r\n>DQ168575.G\r\nMKSLEISAIVGLIVAFIAAIVVWTIVLIEYRKIRKQKRIDKILDRIRERAEDSGNESEGDTEELATLVDMVDFEPWVGDNL\r\n>AY795907.D\r\nMQTLEILSIVALVIAAIIAIIVWTIVYIEYRKIRRQRKIDQLIDRIRERAEDSGNESEGDEEELSTLMEMGHAAPWNVADDL\r\n", "SP_COL_CUTOFF" : "0.93", "SP_SEQ_CUTOFF": "0.6", "maxiterate":"0", "outorder": "aligned", "GENCODE": "1", "F": "+F", "solved" : "false", "pair" : ""}
    forms.append(form)
    
    new_process_id = manager.get_new_process_id()
    job_logger = None
    
    # guidance state part
    try:
    
        guidance_state = None
        
        '''
        if not recaptcha.verify():
            raise Exception ("Run failed to pass I'm not a robot test", "user")
        '''
        
        form = forms[index]
        warning_messages = ''
        # create Job Manager instance
        guidance_state = GuidanceState(jobId = new_process_id, form = form, isRequest = False)
        job_logger = get_job_logger(new_process_id)
        
        # upload files
        guidance_state.upload_files()
        msg = guidance_state.update_state(state = State.Init)
        if msg != '':
            warning_messages = warning_messages + msg
        
        # validate input
        status, msg = guidance_state.validateInput()
        job_logger.info(f'validated data, status = {status}')
        
        if status != 'OK':
            kwargs = {
                "var": guidance_state.var,
            }
            return render_template('error_page.html', error_text=msg, **kwargs) # change
        if msg != '':
            warning_messages = warning_messages + msg
        job_logger.info(f'validated data')
        
        
        # store elements
        job_logger.info(f'storing state')
        if 'crash' in form['JOB_TITLE']:
            crash_flag = True
        else:
            crash_flag = False
        status = guidance_state.save_state(warning_messages, crash_flag=crash_flag)
        if status != 'OK':
            raise Exception (status, "system")
        
    except Exception as e:
    
        if len(e.args) > 0:
            logger.info(e.args[0])
            if job_logger: 
                job_logger.info(e.args[0])
            
        # html FAIL message and footer
        if len(e.args) > 1 and e.args[1]=="user":
            msg = e.args[0]
            error_type = "user"
        else:
            msg = CONSTS.SYS_ERROR_MSG
            error_type = "system"
        
        if len(e.args) > 1:
            log_error_msg = e.args[0]
        else:
            log_error_msg = 'no error message'
            
        log_msg = f'Guidance State caught exception: msg = {log_error_msg}, error_type = {error_type}'
        logger.info(log_msg)
        if job_logger:
            job_logger.info(log_msg)
        
        if guidance_state: 
            guidance_state.update_state(state = State.Error, error_msg = msg, error_type = error_type)
            if error_type == "system":
                guidance_state.send_system_error_email()
            kwargs = {
                "var": guidance_state.var,
            }
        return render_template('error_page.html', error_text=msg, **kwargs)
            
    # here the job is submitted
    
    email_address =  form['email_add']
    if 'norun' not in form['JOB_TITLE']:
        log_msg = f'submitting job: {new_process_id} {email_address} {form["JOB_TITLE"]}'
        logger.info(log_msg)
        job_logger.info(log_msg)
        man_results = manager.add_guidance_process(new_process_id, email_address, form['JOB_TITLE'])
        logger.info('returned from add_guidance_process')
        job_logger.info('returned from add_guidance_process')
    else:
        man_results = False
        
    if not man_results:
        msg = f'job_manager_api can\'t add process {new_process_id}'
        logger.warning(msg)
        kwargs = {
            "var": guidance_state.var,
        }
        if 'norun' not in form['JOB_TITLE']:
            return render_template('error_page.html', error_text=CONSTS.SYS_ERROR_MSG, **kwargs)
        else:
            return render_template('error_page.html', error_text=new_process_id, **kwargs)
            
    log_msg = f'process added man_results = {man_results}, redirecting url'
    logger.info(log_msg)
    job_logger.info(log_msg)
    
    return redirect(url_for('process_state', process_id=new_process_id))
        
@app.route(PREFIX + "/monitor")
def monitor():
    from Job_Manager_Thread_Safe import Job_Manager_Thread_Safe
    
    processes_state_dict = manager.get_processes_state_dict()
    waiting_list = manager.get_waiting_list()
    
    running_jobs_list = []
    waiting_jobs_list = []
    count_finished = 0
    count_crashed = 0
    count_init = 0
    
    for process_id in list(processes_state_dict):
        state = processes_state_dict[process_id].get_job_state("guidance")
        if state == State.Finished:
            count_finished += 1
        elif state == State.Crashed:
            count_crashed += 1
        elif state == State.Init:
            count_init += 1
        else:
            running_jobs_list.append (f'{process_id}, {state}') 
            
    for process_id in waiting_list:
        waiting_jobs_list.append (f'{process_id[0]}, State.Waiting')
            
    return render_template('monitor.html', count_finished = count_finished, count_crashed = count_crashed, count_init = count_init, running_jobs_list = running_jobs_list, waiting_jobs_list = waiting_jobs_list)

@app.route(PREFIX + "/clear_waiting_list")
def clear_waiting_list():
    waiting_list = manager.get_waiting_list()
    n_waiting = len(waiting_list)
    for i in range(n_waiting):
        manager.add_process_from_waiting_list()
    return render_template('msg.html', len = n_waiting)
        
@app.route(PREFIX + "/daily_test", methods=['GET'])
def daily_test():

    FASTA_txt = CONSTS.DAILY_TEST_SEQUENCE
    return render_template('home.html', FASTA_txt = FASTA_txt, daily_test = 'yes')

if __name__ == "__main__":
    # to see in the browser use http://127.0.0.1:5000/guidance/
    # app.run(debug=True)
    app.config['APPLICATION_ROOT'] = '/guidance'  # Ksenia
    app.run(debug=True, port=3000)     # Ksenia
