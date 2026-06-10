import os.path
import random
import subprocess
from datetime import datetime
import os
import re
import sys
import zipfile
import time
import smtplib
from email.message import EmailMessage
from email.utils import formataddr
from subprocess import run, PIPE
import logging
from guidance3 import constants as CONST
from guidance3.utils.timing import timeit

#@timeit
def sample_from_empirical_distribution(distribution_file_name, out_sample_file_name, sample_size):
    op_vals = []
    op_density = []
    op_prob = []
    op_CDF = []
    try:
        with open(distribution_file_name, 'r') as file:
            lines = file.readlines()
            for line in lines[1:]:  # Skip the first line
                values = line.split()
                op_vals.append(float(values[1]))
                op_density.append(float(values[2]))
    except Exception as e:
        raise Exception(f"Guidance::SampleFromEmpiricDistribution:cannot open IN: {distribution_file_name} {e}\n")

    for i in range(len(op_vals) - 1):
        op_prob.append(op_density[i] * (op_vals[i + 1] - op_vals[i]))
        op_CDF.append(sum(op_prob[:i + 1]))

    try:
        with open(out_sample_file_name, 'a') as out_file:
            j = 0
            while j < sample_size:
                rand_num = random.random()
                i = 0
                while rand_num > op_CDF[i]:
                    i += 1
                rand_num = random.random()
                dist_rand_num = op_vals[i] + rand_num * (op_vals[i + 1] - op_vals[i])
                out_file.write(f"{dist_rand_num}\n")
                j += 1
    except Exception as e:
        raise Exception(f"Guidance::SampleFromEmpiricDistribution:cannot open OUT: {out_sample_file_name} {e}\n")

    return op_vals

#@timeit
def sample_from_uniform_dist(start, end, out_sample_file_name, sample_size):
    sample = []
    try:
        with open(out_sample_file_name, 'a') as out_file:
            j = 0
            while j < sample_size:
                op_rand = random.uniform(start, end)
                sample.append(op_rand)
                out_file.write(f"{op_rand}\n")
                j += 1
    except Exception as e:
        raise Exception(f"Guidance::sample_from_uniform_dist:cannot open OUT: {out_sample_file_name} {e}\n")
    return sample


def flag_that_finished_ok(config):
    if config.isServer == 1:
        ends_ok_path = os.path.join(config.WorkingDir, f"GUIDANCE_{config.run_number}.END_OK")
    else:
        ends_ok_path = os.path.join(config.WorkingDir, "ENDS_OK")

    with open(ends_ok_path, "w"):
        pass

    if (config.PROGRAM in {"GUIDANCE", "GUIDANCE3"}) and (config.Seq_Type == "Codons"):
        prot_scr_path = os.path.join(config.WorkingDir, f"{config.Output_Prefix}_res_pair.PROT.scr")
        prot_res_scr_path = os.path.join(config.WorkingDir, f"{config.Output_Prefix}_res_pair_res.PROT.scr")

        if (os.path.getsize(prot_scr_path) / 1048576) > 100:
            with zipfile.ZipFile(f"{prot_scr_path}.zip", 'w', zipfile.ZIP_DEFLATED) as zip_file:
                zip_file.write(prot_scr_path)
            os.remove(prot_scr_path)

        if (os.path.getsize(prot_res_scr_path) / 1048576) > 100:
            with zipfile.ZipFile(f"{prot_res_scr_path}.zip", 'w', zipfile.ZIP_DEFLATED) as zip_file:
                zip_file.write(prot_res_scr_path)
            os.remove(prot_res_scr_path)


def send_administrator_mail_on_error(message, config):
    email_subject = f"SYSTEM ERROR has occurred on GUIDANCE: {config.run_url}"
    email_message = f"Hello,\\n\\nUnfortunately a system SYSTEM ERROR has occurred on GUIDANCE: {config.run_url}.\\nERROR: {message}."
    admin_email = CONST.ADMIN_EMAIL
    # Activate in case the cluster node fails to communicate with the net
    msg = "{}/sendEmail.pl -f 'bioSequence@tauex.tau.ac.il' -t '{}' -u '{}' -xu '{}' -xp '{}' -s '{}' -m '{}'".format(
        config.send_email_dir,
        "bioSequence@tauex.tau.ac.il",
        email_subject,
        config.userName,
        config.userPass,
        config.smtp_server,
        email_message
    )
    os.chdir(config.send_email_dir)
    email_system_return = subprocess.getoutput(msg)
    return email_system_return


def print_time():
    now = datetime.now()
    formatted_time = now.strftime("%H:%M:%S %d-%m-%Y")
    return formatted_time

def exit_on_error(which_error, error_msg, config):
    error_definition = "<font size=+1 color='red'>ERROR! GUIDANCE session has been terminated:</font><br />\n"
    sys_error = "<font size=+1 color='red'>A SYSTEM ERROR OCCURRED!</font><br />Please try to run GUIDANCE again in a few minutes.<br />We apologize for the inconvenience.<br />\n"

    if config.isServer != 1:
        sys_error = "Guidance error\n"
        error_definition = "Guidance error: "

    if config.isServer == 1:
        with open(config.server_output, 'a') as output_file, open(
                f'{config.OutLogFile}', "a") as log_file:
            if which_error == 'user_error':
                log_file.write(f"\nEXIT on error:\n{error_msg}\n")
                output_file.write(f"{error_definition} {error_msg}")
            elif which_error == 'sys_error':
                send_administrator_mail_on_error(error_msg, config)
                log_file.write(f"\n{error_msg}\n")
                output_file.write(f"{sys_error}")

        # Finish the output page
        time.sleep(10)
        with open(os.path.join(config.WorkingDir, config.server_output), 'r') as output_file:
            output = output_file.readlines()

        # Remove the refresh commands from the output page
        with open(os.path.join(config.WorkingDir, config.server_output), 'w') as output_file:
            for line in output:
                if "TTP-EQUIV=\"REFRESH\"" in line or "CONTENT=\"NO-CACHE\"" in line:
                    continue
                elif re.match(r'(.*)RUNNING(.*)', line):
                    output_file.write(
                        re.match(r'(.*)RUNNING(.*)', line).group(1) + "FAILED" + re.match(r'(.*)RUNNING(.*)',
                                                                                          line).group(2))
                else:
                    output_file.write(line)

        try:
            with open(f'{config.WorkingDir}/errors.txt', "w") as f:
                f.write(error_msg.replace("<br>", "\n"))
        except Exception as e:
            with open(f'{config.OutLogFile}', "a") as log_file:
                log_file.write(f"Could not write errors.txt: {e}\n")

        if config.user_mail != "": #TODO change on the server
            send_mail_on_error(config)

        with open(f'{config.OutLogFile}', "a") as log_file:
            log_file.write(f"\nExit Time: {print_time()}\n")
        os.chmod(config.WorkingDir, 0o755)


    else:
        if which_error == 'user_error':
            with open(f'{config.OutLogFile}', "a") as log_file:
                log_file.write(f"\nEXIT on error:\n{error_msg}\n")
            print(f"ERROR: {error_msg}")
        elif which_error == 'sys_error':
            with open(f'{config.OutLogFile}', "a") as log_file:
                log_file.write(f"\n{error_msg}\n")
            print(f"ERROR: {error_msg}")
            print(sys_error + "\n")

    if config.PROGRAM == "GUIDANCE" and config.isServer == 1:  # Zip BP dir on server
        # Tar and remove the BP dir
        cmd = f"cd {config.WorkingDir};tar -czf {config.Output_Prefix}_BP_Dir.tar.gz ./BP"
        # print(f"{cmd}\n")
        os.system(cmd)
        if os.path.exists(f"{config.Output_Prefix}_BP_Dir.tar.gz"):
            os.system(f"rm -r -f {config.BootStrap_Dir}")

    sys.exit(1)


def send_mail_on_error(config):
    email_subject = "Your GUIDANCE run for {} FAILED".format(config.usrSeq_File)
    HttpPath = "{}{}".format(config.run_url, config.output_page)
    email_message = "Hello,\n\nUnfortunately your GUIDANCE run (number {}) has failed.\nPlease have a look at {} for further details\n\nSorry for the inconvenience\nGUIDANCE Team".format(
        config.run_number, HttpPath)

    send_email_script = f'{CONST.BIN_DIR}/perl/sendEmail.pl'

    cmd = [
        send_email_script,
        '-f', CONST.ADMIN_EMAIL,
        '-t', config.user_mail, #TODO - config.user_mail
        '-u', email_subject,
        '-xu', config.userName,
        '-xp', config.userPass,
        '-s', config.smtp_server,
        '-m', email_message
    ]

    cmd = ' '.join(cmd)

    with open(f'{config.OutLogFile}', "a") as log_file:
        log_file.write(f"MESSAGE:{email_message}\nCOMMAND:{cmd}\n")
        os.chdir(f'{config.send_email_dir}')
        # Execute the command
        try:
            email_system_return = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            if 'successfully' not in email_system_return.stdout:
                log_file.write(
                    f"send_mail: The message was not sent successfully. system returned: {email_system_return.stdout}\n")
        except subprocess.TimeoutExpired:
            log_file.write("send_mail: timed out after 30s — SMTP server may be unreachable\n")
        except Exception as e:
            log_file.write(f"send_mail: failed — {e}\n")

def subtract_time_from_now(begin_time_str, time_str):
    """Return elapsed time string since begin_time_str + time_str.

    begin_time_str is expected in the format 'HH:MM:SS' and time_str in 'DD-MM-YYYY',
    producing a combined string 'HH:MM:SS DD-MM-YYYY' that is parsed and compared
    against the current time.
    """
    combined = f"{begin_time_str} {time_str}"
    match = re.match(r'(\d+):(\d+):(\d+) (\d+)-(\d+)-(\d+)', combined)
    if not match:
        return f'error: could not parse time string: {combined}\n'
    hour, minute, second, day, month, year = match.groups()
    try:
        begin = datetime(int(year), int(month), int(day), int(hour), int(minute), int(second))
    except ValueError as e:
        return f'error: invalid date/time values: {e}\n'
    delta = datetime.now() - begin
    total_seconds = int(delta.total_seconds())
    if total_seconds < 0:
        return 'error: begin time is in the future\n'
    hours, remainder = divmod(total_seconds, 3600)
    minutes, secs = divmod(remainder, 60)
    return f"{str(hours).zfill(2)}:{str(minutes).zfill(2)}:{str(secs).zfill(2)}"


def compare_time(time1_dict, time2_dict):
    """Compare two time dictionaries and return elapsed time string.

    Each dict has keys: Year, Month, Day, Hour, Minute, Second (string or int values).
    Returns ("yes", elapsed_str) on success or an error string on failure.
    """
    try:
        t1 = datetime(
            int(time1_dict['Year']), int(time1_dict['Month']), int(time1_dict['Day']),
            int(time1_dict['Hour']), int(time1_dict['Minute']), int(time1_dict['Second'])
        )
        t2 = datetime(
            int(time2_dict['Year']), int(time2_dict['Month']), int(time2_dict['Day']),
            int(time2_dict['Hour']), int(time2_dict['Minute']), int(time2_dict['Second'])
        )
    except (ValueError, KeyError) as e:
        return f'error: invalid time dict values: {e}\n'
    delta = t2 - t1
    total_seconds = int(delta.total_seconds())
    if total_seconds < 0:
        return f'error: t1 is after t2\n'
    hours, remainder = divmod(total_seconds, 3600)
    minutes, secs = divmod(remainder, 60)
    elapsed = f"{str(hours).zfill(2)}:{str(minutes).zfill(2)}:{str(secs).zfill(2)}"
    return ("yes", elapsed)


def convert_current_time(date_dictionary):
    """Populate date_dictionary with the current date/time components."""
    current_time = datetime.now()
    date_dictionary['Year'] = current_time.year
    date_dictionary['Month'] = current_time.month
    date_dictionary['Day'] = current_time.day
    date_dictionary['Hour'] = current_time.hour
    date_dictionary['Minute'] = current_time.minute
    date_dictionary['Second'] = current_time.second


def calculate_time_difference(begin_time_str):
    """Return human-readable elapsed time since begin_time_str (format: 'DD/MM/YYYY HH:MM:SS')."""
    begin = datetime.strptime(begin_time_str, '%d/%m/%Y %H:%M:%S')
    delta = datetime.now() - begin
    hours, remainder = divmod(int(delta.total_seconds()), 3600)
    minutes, seconds = divmod(remainder, 60)
    days = delta.days
    return f"{days}d {hours}h {minutes}m {seconds}s" if days > 0 else f"{hours}h {minutes}m {seconds}s"


def convert_num(input_num):
    if input_num < 10:
        return f"0{input_num}"
    else:
        return str(input_num)


def print_message_to_output(msg, config):
    try:
        # with open(f"{config.WorkingDir}{config.server_output}", "a") as output_file:
        with open(f"{config.server_output}", "a") as output_file:
            output_file.write(f"\n<ul><li>{msg}</li></ul>\n")
    except Exception as e:
        print(f"Failed to open output file: {e}\n")
        sys.exit()


def print_initial_running_progress(config):
    config.progress_report = "ProgressReport.html"
    config.alt_msa_status = "MSA_STATUS.txt"

    with open(config.WorkingDir + config.alt_msa_status, "a") as ALT_STATUS:
        ALT_STATUS.write("<ul class=\"in_progress\"><li>Generating alternative alignments</li></ul>\n")

    with open(config.WorkingDir + config.progress_report, "a") as PROGRESS:
        PROGRESS.write("<p><font face=Verdana size=2>\n")

        if config.Redirect_From_MAFFT != "1":
            PROGRESS.write("<ul class=\"in_progress\"><li>Generating the base alignment</li></ul>\n")

        if config.PROGRAM == "HoT":
            PROGRESS.write("<ul class=\"in_progress\"><li>Constructing guide tree</li></ul>\n")
        else:
            PROGRESS.write("<ul class=\"in_progress\"><li>Constructing bootstrap guide-trees</li></ul>\n")
        PROGRESS.write("<ul class=\"in_progress\"><li>Generating alternative alignments</li></ul>\n")

        if config.PROGRAM == "GUIDANCE":
            PROGRESS.write("<ul class=\"in_progress\"><li>Calculating GUIDANCE scores</li></ul>\n")

        if config.PROGRAM == "HoT":
            PROGRESS.write("<ul class=\"in_progress\"><li>Calculating HoT scores</li></ul>\n")

        if config.PROGRAM == "GUIDANCE3":
            PROGRESS.write("<ul class=\"in_progress\"><li>Calculating GUIDANCE3 scores</li></ul>\n")
            PROGRESS.write("<ul class=\"in_progress\"><li>Running the model and selecting the best MSA</li></ul>\n")

        PROGRESS.write("</font>\n")

#@timeit
def update_progress(progress_file, message):
    with open(progress_file, "r") as progress:
        data = progress.readlines()
    with open(progress_file, "w") as progress:
        for line in data:
            if message in line:
                line = line.replace("in_progress", "finished")
                if "(estimated time" in line:
                    line = line.replace("(estimated time", "").split(")")[1]
                progress.write(line)
            elif "Started generating alternative alignments" in message and "Generating alternative alignments" in line:
                line = line.replace("Generating alternative alignments", message)
                progress.write(line)
            elif "Finished generating" in message and "Started generating alternative alignments" in line:
                line = line.replace("in_progress", "finished")
                line = line.replace("Started generating alternative alignments", message)
                progress.write(line)
            elif "Finished Calculating" in message and "Calculating" in line:
                line = line.replace("in_progress", "finished")
                line = line.replace("Calculating", "Finished Calculating")
                progress.write(line)
            elif "Finished running the model" in message and "Running the model" in line:
                line = line.replace("in_progress", "finished")
                line = line.replace("Running the model", "Finished running the model")
                progress.write(line)
            else:
                progress.write(line)

def send_finish_email_to_user(config):
    # Set up logging
    logging.basicConfig(filename=f"{config.WorkingDir}/log.txt", level=logging.INFO)

    email_subject = ""
    # http_path = f"http://guidance-dev.tau.ac.il/results/{vars['run_number']}"
    base_http_path = "http://guidance-dev.tau.ac.il/"
    http_path = base_http_path + f"/guidance/results/{config.run_number}"

    if config.JOB_TITLE:
        email_subject = f"Your Guidance results for {config.JOB_TITLE} are ready"
    elif config.usrSeq_File:
        email_subject = f"Your Guidance results for {config.usrSeq_File} are ready"
    else:
        email_subject = f"Your Guidance results for run number {config.run_number} are ready"

    email_message = f"""Hello,

The results for your Guidance run are ready at:
{http_path}

Running Parameters:
Job Title: {config.JOB_TITLE}
Sequences File: {config.usrSeq_File}
MSA Algorithm: {config.MSA_Program}
Number of Bootstraps: {config.Bootstraps}
Scoring Method: {config.PROGRAM}

Please note: the results will be kept on the server for three months.

Thanks,
GUIDANCE Team"""

    # Set up email
    msg = EmailMessage()
    msg.set_content(email_message)
    msg['Subject'] = email_subject
    msg['From'] = formataddr(('GUIDANCE Team', 'admin@example.com'))  # Replace with ADMIN_EMAIL
    msg['To'] = config.user_mail

    # Send email
    try:
        with smtplib.SMTP(config.smtp_server) as server:
            server.login(config.userName, config.userPass)
            server.send_message(msg)
            logging.info("Email sent successfully.")
    except Exception as e:
        logging.error(f"Failed to send email: {e}")

    # Log the email message and command
    log_msg = f"MESSAGE: {email_message}\n"
    logging.info(log_msg)

    # Run external command (if needed)
    email_command = [
        'perl', f'{CONST.BIN_DIR}/perl/sendEmail.pl',
        '-f', 'admin@example.com',  # Replace with GENERAL_CONSTANTS::ADMIN_EMAIL
        '-t', config.user_mail,
        '-u', email_subject,
        '-xu', config.userName,
        '-xp', config.userPass,
        '-s', config.smtp_server,
        '-m', email_message
    ]

    result = run(email_command, stdout=PIPE, stderr=PIPE, text=True)
    if "successfully" not in result.stdout:
        logging.error(f"send_mail: The message was not sent successfully. System returned: {result.stdout}")

