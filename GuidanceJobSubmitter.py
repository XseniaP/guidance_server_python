import pathlib
import subprocess
import os
import sys
from SharedConsts import EMAIL_CONSTS
# from script.config import EMAIL_CONSTS
import SharedConsts as CONSTS
# import script.config as CONSTS
from utils import logger, get_job_logger

sys.path.append('/guidance/guidance_server_python/script/')
from submit_slurm import submit_job_to_Q
from email_sender import send_email

class GuidanceJobSubmitter:
    
    @staticmethod
    def submit_job(working_dir, run_parameters, is_daily_test=False):
        
        job_id = str(pathlib.Path(working_dir).stem)
        parameters = f"{os.path.join(working_dir, 'VARS.json')} {os.path.join(working_dir, 'FORM.json')}"
        to_email = run_parameters
        cmds_file = os.path.join(working_dir, 'qsub.cmds')
        GuidanceJobSubmitter._write_cmds_file(cmds_file, parameters, working_dir, job_id, to_email)
        
        job_id_file = os.path.join(working_dir, 'QSTAT_NO') 

        # a simple command when using shebang header (#!) in q_submitter_power.py
        # replace by JS
        submission_cmd = f'cd {working_dir}\n python3 {CONSTS.MAIN_SCRIPT} {parameters} > {os.path.join(working_dir, "std.out")}'
        if to_email:
            email_cmd = f"python3 {CONSTS.EMAIL_SCRIPT} {CONSTS.SMTP_SERVER} {CONSTS.ADMIN_EMAIL} {to_email} --subject '{EMAIL_CONSTS.FINISHED_TITLE}' --content '{EMAIL_CONSTS.FINISHED_CONTENT.format(results_url=CONSTS.WEBSERVER_RESULTS_URL_EXT, process_id=job_id)}'"
            submission_cmd = f'{submission_cmd}\n{email_cmd}'
        if is_daily_test:
            write_daily_test_cmd = f"python {CONSTS.WRITE_DAILY_TEST_SCRIPT} {CONSTS.DAILY_TEST_DIR} {job_id}"
            submission_cmd = f'{submission_cmd}\n {write_daily_test_cmd}'
        #terminal_cmd = f'/opt/pbs/bin/qsub {str(temp_script_path)}'
        
        #send_email(CONSTS.SMTP_SERVER, CONSTS.ADMIN_EMAIL,
        #                to_email, subject=EMAIL_CONSTS.INIT_TITLE,
        #                content= EMAIL_CONSTS.INIT_CONTENT.format(results_url=CONSTS.WEBSERVER_PROCESS_STATE_URL_EXT, process_id=job_id))
        # send_email(CONSTS.SMTP_SERVER, CONSTS.ADMIN_EMAIL,
        #                CONSTS.ADMIN_EMAIL, subject=f'{EMAIL_CONSTS.INIT_TITLE} by {to_email}',
        #                content= EMAIL_CONSTS.INIT_CONTENT.format(results_url=CONSTS.WEBSERVER_PROCESS_STATE_URL_EXT, process_id=job_id))
        pid = os.fork()
        if pid == 0:
            #print(submission_cmd)
            #logger.info(f'submission_cmd={submission_cmd}')
            job_logger = get_job_logger(job_id)
            if job_logger:
                job_logger.info(f'submission_cmd={submission_cmd}')
            #job_run_output = subprocess.run(submission_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            try:
                job_run_output = submit_job_to_Q( working_dir, submission_cmd)
            except:
                if job_logger:
                    job_logger.info('submit_job_to_Q failed')
            if job_logger:
                job_logger.info(f'process id ={job_run_output}')
            # print(job_run_output.stdout)
                      
            # check return code and raise exception if failed
            if not job_run_output:
                error = f'submit job {job_id} failed:{job_run_output.returncode}'
                logger.error(error)
                job_logger = get_job_logger(job_id)
                if job_logger:
                    job_logger.error(error)
                # create failure file
                fail_page = os.path.join( CONSTS.WEBSERVER_RESULTS_DIR, job_id, f'GUIDANCE_{job_id}.END_FAIL')
                open(fail_page, 'a').close()
                send_email(CONSTS.SMTP_SERVER, CONSTS.ADMIN_EMAIL,
                       to_email, subject=EMAIL_CONSTS.CRASHED_TITLE,
                       content= EMAIL_CONSTS.CRASHED_CONTENT.format(results_url=CONSTS.WEBSERVER_RESULTS_URL_EXT, process_id=job_id))
                send_email(CONSTS.SMTP_SERVER, CONSTS.ADMIN_EMAIL,
                       CONSTS.ADMIN_EMAIL, subject=EMAIL_CONSTS.CRASHED_TITLE,
                       content= EMAIL_CONSTS.CRASHED_CONTENT.format(results_url=CONSTS.WEBSERVER_RESULTS_URL_EXT, process_id=job_id))
            else:
                #send_email(CONSTS.SMTP_SERVER, CONSTS.ADMIN_EMAIL, to_email,
                #     subject=EMAIL_CONSTS.FINISHED_TITLE, content= EMAIL_CONSTS.FINISHED_CONTENT.format(results_url=CONSTS.WEBSERVER_RESULTS_URL_EXT, process_id=job_id))
                if to_email: 
                    send_email(CONSTS.SMTP_SERVER, CONSTS.ADMIN_EMAIL,
                            to_email, subject=EMAIL_CONSTS.INIT_TITLE,
                            content= EMAIL_CONSTS.INIT_CONTENT.format(results_url=CONSTS.WEBSERVER_PROCESS_STATE_URL_EXT, process_id=job_id))
            return 0 #child process
        else:
            return 1 #parent process
            
        
            
        #return job_run_output.stdout.decode('utf-8').split('.')[0], ''
	
    @staticmethod
    def _write_cmds_file(cmds_file, parameters, working_dir, run_number, to_email):
    
        # the queue does not like very long commands so I use a dummy delimiter (!@#) to break the commands for q_submitter
        new_line_delimiter = '!@#'

        required_modules_as_str = ' '.join(CONSTS.REQUIRED_MODULES)
        with open(cmds_file, 'w') as f:
            f.write(f'module load mamba/mamba-2.1.1;')
            f.write('bash activate /guidance/guidance_server_python/guidance_env;')
            f.write(f'cd {working_dir};')
            f.write(new_line_delimiter)
            f.write(f'python3 {CONSTS.MAIN_SCRIPT} {parameters} > {os.path.join(working_dir, "std.out")}\t{CONSTS.GUIDANCE_JOB_PREFIX}_{run_number}')
            if to_email:
                f.write(new_line_delimiter)
                f.write(f"python3 {CONSTS.EMAIL_SCRIPT} {CONSTS.SMTP_SERVER} {CONSTS.ADMIN_EMAIL} {to_email} --subject '{EMAIL_CONSTS.FINISHED_TITLE}' --content '{EMAIL_CONSTS.FINISHED_CONTENT.format(results_url=CONSTS.WEBSERVER_RESULTS_URL_EXT, process_id=run_number)}'")
            f.write('\n')
        f.close()
