import pathlib
import subprocess
import os
import sys
from SharedConsts import EMAIL_CONSTS
import SharedConsts as CONSTS
from utils import logger, get_job_logger, load_running_parameters_text

# script/ lives next to this file both locally and on the server (MAIN_SCRIPT and
# EMAIL_SCRIPT already resolve it the same way via CONSTS.SCRIPTS_DIR), so reuse that
# constant here instead of hardcoding the server's absolute path.
if CONSTS.SCRIPTS_DIR not in sys.path:
    sys.path.append(CONSTS.SCRIPTS_DIR)
from email_sender import send_email

# The TAU SLURM REST API (script/submit_slurm.py) is only reachable from the TAU
# network and requires SLURM_API_KEY, so there's no queue to submit to from a local
# dev machine. Set USE_SLURM_QUEUE=false in the local .env to run jobs as a plain
# subprocess instead; it defaults to true so the server keeps submitting to SLURM.
USE_SLURM_QUEUE = os.getenv('USE_SLURM_QUEUE', 'true').strip().lower() not in ('0', 'false', 'no')

if USE_SLURM_QUEUE:
    from submit_slurm import submit_job_to_Q


class GuidanceJobSubmitter:

    @staticmethod
    def submit_job(working_dir, run_parameters, is_daily_test=False):

        job_id = str(pathlib.Path(working_dir).stem)
        parameters = f"{os.path.join(working_dir, 'VARS.json')} {os.path.join(working_dir, 'FORM.json')}"
        to_email = run_parameters
        cmds_file = os.path.join(working_dir, 'qsub.cmds')
        GuidanceJobSubmitter._write_cmds_file(cmds_file, parameters, working_dir, job_id)

        # a simple command when using shebang header (#!) in q_submitter_power.py
        # replace by JS
        submission_cmd = f'cd {working_dir}\n python3 {CONSTS.MAIN_SCRIPT} {parameters} > {os.path.join(working_dir, "std.out")}'
        if is_daily_test:
            write_daily_test_cmd = f"python {CONSTS.WRITE_DAILY_TEST_SCRIPT} {CONSTS.DAILY_TEST_DIR} {job_id}"
            submission_cmd = f'{submission_cmd}\n {write_daily_test_cmd}'

        pid = os.fork()
        if pid == 0:
            job_logger = get_job_logger(job_id)
            if job_logger:
                job_logger.info(f'submission_cmd={submission_cmd}')

            if USE_SLURM_QUEUE:
                try:
                    job_run_output = submit_job_to_Q(working_dir, submission_cmd)
                except Exception as e:
                    job_run_output = ''
                    if job_logger:
                        job_logger.info(f'submit_job_to_Q failed: {e}')
            else:
                # No queue locally: run the pipeline synchronously (blocks until it's
                # done) and use the exit code as the success/failure signal, matching
                # the truthy job-id / falsy '' contract submit_job_to_Q returns.
                result = subprocess.run(submission_cmd, shell=True)
                job_run_output = job_id if result.returncode == 0 else ''

            if job_logger:
                job_logger.info(f'process id ={job_run_output}')

            if not job_run_output:
                error = f'submit job {job_id} failed'
                logger.error(error)
                if job_logger:
                    job_logger.error(error)
                # create failure file
                fail_page = os.path.join(CONSTS.WEBSERVER_RESULTS_DIR, job_id, f'GUIDANCE_{job_id}.END_FAIL')
                open(fail_page, 'a').close()
                if to_email:
                    send_email(CONSTS.SMTP_SERVER, CONSTS.ADMIN_EMAIL,
                           to_email, subject=EMAIL_CONSTS.CRASHED_TITLE,
                           content= EMAIL_CONSTS.CRASHED_CONTENT.format(results_url=CONSTS.WEBSERVER_RESULTS_URL_EXT, process_id=job_id))
                send_email(CONSTS.SMTP_SERVER, CONSTS.ADMIN_EMAIL,
                       CONSTS.ADMIN_EMAIL, subject=EMAIL_CONSTS.CRASHED_TITLE,
                       content= EMAIL_CONSTS.CRASHED_CONTENT.format(results_url=CONSTS.WEBSERVER_RESULTS_URL_EXT, process_id=job_id))
            elif USE_SLURM_QUEUE and to_email:
                # Queued asynchronously: let the user know it's been submitted, with the
                # full running parameters. The "results are ready" email is sent
                # separately by the pipeline itself once the job actually finishes on
                # the compute node (guidance3.utils.common.send_finish_email_to_user).
                content = EMAIL_CONSTS.INIT_CONTENT.format(results_url=CONSTS.WEBSERVER_PROCESS_STATE_URL_EXT, process_id=job_id)
                running_params = load_running_parameters_text(working_dir)
                if running_params:
                    content += f"\n\nRunning Parameters:\n{running_params}"
                send_email(CONSTS.SMTP_SERVER, CONSTS.ADMIN_EMAIL,
                        to_email, subject=EMAIL_CONSTS.INIT_TITLE,
                        content=content)
            return 0 #child process
        else:
            return 1 #parent process

    @staticmethod
    def _write_cmds_file(cmds_file, parameters, working_dir, run_number):

        # the queue does not like very long commands so I use a dummy delimiter (!@#) to break the commands for q_submitter
        new_line_delimiter = '!@#'

        required_modules_as_str = ' '.join(CONSTS.REQUIRED_MODULES)
        with open(cmds_file, 'w') as f:
            f.write(f'module load {required_modules_as_str};')
            f.write('bash activate /guidance/guidance_server_python/guidance_env;')
            f.write(f'cd {working_dir};')
            f.write(new_line_delimiter)
            f.write(f'python3 {CONSTS.MAIN_SCRIPT} {parameters} > {os.path.join(working_dir, "std.out")}\t{CONSTS.GUIDANCE_JOB_PREFIX}_{run_number}')
            f.write('\n')
        f.close()
