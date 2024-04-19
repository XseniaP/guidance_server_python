import os.path
import random
import subprocess
from datetime import datetime
import os
import re
import sys
import zipfile
import time
import guidance_CONSTANTS
from time_decorator import timeit
@timeit
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

@timeit
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


def flag_that_finished_ok(args_library):
    if args_library.isServer == 1:
        ends_ok_path = os.path.join(args_library.WorkingDir, f"GUIDANCE_{args_library.run_number}.END_OK")
    else:
        ends_ok_path = os.path.join(args_library.WorkingDir, "ENDS_OK")

    with open(ends_ok_path, "w"):
        pass

    if (args_library.PROGRAM in {"GUIDANCE", "GUIDANCE2"}) and (args_library.Seq_Type == "Codons"):
        prot_scr_path = os.path.join(args_library.WorkingDir, f"{args_library.Output_Prefix}_res_pair.PROT.scr")
        prot_res_scr_path = os.path.join(args_library.WorkingDir, f"{args_library.Output_Prefix}_res_pair_res.PROT.scr")

        if (os.path.getsize(prot_scr_path) / 1048576) > 100:
            with zipfile.ZipFile(f"{prot_scr_path}.zip", 'w', zipfile.ZIP_DEFLATED) as zip_file:
                zip_file.write(prot_scr_path)
            os.remove(prot_scr_path)

        if (os.path.getsize(prot_res_scr_path) / 1048576) > 100:
            with zipfile.ZipFile(f"{prot_res_scr_path}.zip", 'w', zipfile.ZIP_DEFLATED) as zip_file:
                zip_file.write(prot_res_scr_path)
            os.remove(prot_res_scr_path)


def send_administrator_mail_on_error(message, args_library):
    email_subject = f"SYSTEM ERROR has occurred on GUIDANCE: {args_library.run_url}"
    email_message = f"Hello,\\n\\nUnfortunately a system SYSTEM ERROR has occurred on GUIDANCE: {args_library.run_url}.\\nERROR: {message}."
    admin_email = guidance_CONSTANTS.ADMIN_EMAIL
    # Activate in case the cluster node fails to communicate with the net
    # msg = "ssh bioseq@lecs \" cd {}; perl sendEmail.pl -f 'bioSequence@tauex.tau.ac.il' -t '{}' -u '{}' -xu '{}' -xp '{}' -s '{}' -m '{}'\"".format(VARS['send_email_dir'], admin_email, email_subject, VARS['userName'], VARS['userPass'], VARS['smtp_server'], email_message)
    msg = "{}/sendEmail.pl -f 'bioSequence@tauex.tau.ac.il' -t '{}' -u '{}' -xu '{}' -xp '{}' -s '{}' -m '{}'".format(
        args_library.send_email_dir,
        "bioSequence@tauex.tau.ac.il",
        email_subject,
        args_library.userName,
        args_library.userPass,
        args_library.smtp_server,
        email_message
    )
    # print("MESSAGE:{}\nCOMMAND:{}".format(email_message, msg))
    os.chdir(args_library.send_email_dir)
    # email_system_return = os.popen(msg).read()
    email_system_return = subprocess.getoutput(msg)
    return email_system_return


def print_time():
    now = datetime.now()
    formatted_time = now.strftime("%H:%M:%S %d-%m-%Y")
    return formatted_time

def exit_on_error(which_error, error_msg, args_library):
    error_definition = "<font size=+1 color='red'>ERROR! GUIDANCE session has been terminated:</font><br />\n"
    sys_error = "<font size=+1 color='red'>A SYSTEM ERROR OCCURRED!</font><br />Please try to run GUIDANCE again in a few minutes.<br />We apologize for the inconvenience.<br />\n"

    if args_library.isServer == 0:
        sys_error = "Guidance error\n"
        error_definition = "Guidance error: "

    if args_library.isServer == 1:
        # with open(os.path.join(args_library.WorkingDir, args_library.server_output), 'a') as output_file, open(
        #         f'{args_library.OutLogFile}', "a") as log_file:
        with open(args_library.server_output, 'a') as output_file, open(
                f'{args_library.OutLogFile}', "a") as log_file:
            if which_error == 'user_error':
                log_file.write(f"\nEXIT on error:\n{error_msg}\n")
                output_file.write(f"{error_definition} {error_msg}")
            elif which_error == 'sys_error':
                send_administrator_mail_on_error(error_msg, args_library)
                log_file.write(f"\n{error_msg}\n")
                output_file.write(f"{sys_error}")

        # Finish the output page
        time.sleep(10)
        with open(os.path.join(args_library.WorkingDir, args_library.server_output), 'r') as output_file:
            output = output_file.readlines()

        # Remove the refresh commands from the output page
        with open(os.path.join(args_library.WorkingDir, args_library.server_output), 'w') as output_file:
            for line in output:
                if "TTP-EQUIV=\"REFRESH\"" in line or "CONTENT=\"NO-CACHE\"" in line:
                    continue
                elif re.match(r'(.*)RUNNING(.*)', line):
                    output_file.write(
                        re.match(r'(.*)RUNNING(.*)', line).group(1) + "FAILED" + re.match(r'(.*)RUNNING(.*)',
                                                                                          line).group(2))
                else:
                    output_file.write(line)

            output_file.write(
                    f"<hr> <h4 class=footer><p align='center'>\nQuestions and comments are welcome! Please <span class=\"admin_link\"><a href=\"mailto:bioSequence@tauex.tau.ac.il?subject=GUIDANCE%20Run%20Number%20{args_library.run_number}\">contact us</a></span></p></h4>\n<div id=\"bottom_links\"> <!-- links before the footer --><span class=\"bottom_link\"> <a href=\"{args_library.home_URL}\" target=\"_blank\">Home</a> &nbsp;|&nbsp;<a href=\"{args_library.overview_URL}\" target=\"_blank\">Overview</a> &nbsp;|&nbsp;<a href=\"{args_library.gallery_URL}\" target=\"_blank\">Gallery</a> &nbsp;|&nbsp;<a href=\"/credits.html\" target=\"_blank\">Credits</a> </span> <br /> </div>")

            output_file.write("</body>\n")
            output_file.write("</html>\n")

        if args_library.user_email != "":
            send_mail_on_error(args_library)

        with open(f'{args_library.OutLogFile}', "a") as log_file:
            log_file.write(f"\nExit Time: {print_time()}\n")
        os.chmod(args_library.WorkingDir, 0o755)

    else:
        if which_error == 'user_error':
            with open(f'{args_library.OutLogFile}', "a") as log_file:
                log_file.write(f"\nEXIT on error:\n{error_msg}\n")
            print(f"ERROR: {error_msg}")
        elif which_error == 'sys_error':
            with open(f'{args_library.OutLogFile}', "a") as log_file:
                log_file.write(f"\n{error_msg}\n")
            print(f"ERROR: {error_msg}")
            print(sys_error + "\n")

    if args_library.PROGRAM == "GUIDANCE" and args_library.isServer == 1:  # Zip BP dir on server
        # Tar and remove the BP dir
        cmd = f"cd {args_library.WorkingDir};tar -czf {args_library.Output_Prefix}_BP_Dir.tar.gz ./BP"
        # print(f"{cmd}\n")
        os.system(cmd)
        if os.path.exists(f"{args_library.Output_Prefix}_BP_Dir.tar.gz"):
            os.system(f"rm -r -f {args_library.BootStrap_Dir}")

    sys.exit()


def send_mail_on_error(args_library):
    email_subject = "Your GUIDANCE run for {} FAILED".format(args_library.usrSeq_File)
    HttpPath = "{}{}".format(args_library.run_url, args_library.output_page)
    email_message = "Hello,\n\nUnfortunately your GUIDANCE run (number {}) has failed.\nPlease have a look at {} for further details\n\nSorry for the inconvenience\nGUIDANCE Team".format(
        args_library.run_number, HttpPath)

    send_email_script = './sendEmail.pl'

    cmd = [
        send_email_script,
        '-f', guidance_CONSTANTS.ADMIN_EMAIL,
        '-t', args_library.user_email,
        '-u', email_subject,
        '-xu', args_library.userName,
        '-xp', args_library.userPass,
        '-s', args_library.smtp_server,
        '-m', email_message
    ]

    cmd = ' '.join(cmd)

    with open(f'{args_library.OutLogFile}', "a") as log_file:
        log_file.write(f"MESSAGE:{email_message}\nCOMMAND:{cmd}\n")
        os.chdir(f'{args_library.send_email_dir}')
        # Execute the command
        email_system_return = subprocess.run(cmd, capture_output=True, text=True)

        # Check if the email was sent successfully
        if 'successfully' not in email_system_return.stdout:
            log_file.write(
                f"send_mail: The message was not sent successfully. system returned: {email_system_return.stdout}\n")

def subtract_time_from_now(begin_time_str, time_str):
    begin_time_str += " " + time_str

    match = re.match(r'(\d+):(\d+):(\d+) (\d+)-(\d+)-(\d+)', begin_time_str)
    if match:
        hour, minute, second, day, month, year = match.groups()

    date1_dict = {'Year': year, 'Month': month, 'Day': day, 'Hour': hour, 'Minute': minute, 'Second': second}
    date2_dict = {'Year': '', 'Month': '', 'Day': '', 'Hour': '', 'Minute': '', 'Second': ''}
    convert_current_time(date2_dict)

    time_difference = compare_time(date1_dict, date2_dict)

    if "error" in time_difference[0]:
        return time_difference[0]
    else:
        return time_difference[1]


def compare_time(time1_dict, time2_dict):
    days_each_month = {'01': 31, '02': 28, '03': 31, '04': 30, '05': 31, '06': 30,
                       '07': 31, '08': 31, '09': 30, '10': 31, '11': 30, '12': 31}
    # time_difference = 0
    # no_of_days_passed = 0

    if time1_dict['Month'] == time2_dict['Month']:  # same month
        if time1_dict['Day'] == time2_dict['Day']:  # same day
            if time2_dict['Hour'] >= time1_dict['Hour']:  # compare hour: h2 > h1
                time_difference = calculate_time_difference(time1_dict['Hour'], time2_dict['Hour'],
                                                            time1_dict['Minute'], time2_dict['Minute'],
                                                            time1_dict['Second'], time2_dict['Second'], 0)
            else:
                return f'error: H1 is: {time1_dict["Hour"]}, H2 is: {time2_dict["Hour"]}. It is the same day, so it is impossible that H1 > H2 \n'
        else:  # different day
            if time2_dict['Day'] >= time1_dict['Day']:
                no_of_days_passed = time2_dict['Day'] - time1_dict['Day']
                time_difference = calculate_time_difference(time1_dict['Hour'], time2_dict['Hour'],
                                                            time1_dict['Minute'], time2_dict['Minute'],
                                                            time1_dict['Second'], time2_dict['Second'],
                                                            no_of_days_passed)
            else:
                return f'error: D1 is: {time1_dict["Day"]}, D2 is: {time2_dict["Day"]}. it is impossible in the same month that D1>D2 \n'
    else:  # different month
        if time2_dict['Month'] - time1_dict['Month'] > 1 or time2_dict['Month'] - time1_dict['Month'] < 0:
            return f'error: M1 is: {time1_dict["Month"]}, M2 is: {time2_dict["Month"]}. The program doesnt allow a difference bigger than 1 month.\n'
        else:  # 1 month difference
            no_of_days_passed = time2_dict['Day'] + days_each_month[time1_dict['Month']] - time1_dict['Day']
            time_difference = calculate_time_difference(time1_dict['Hour'], time2_dict['Hour'], time1_dict['Minute'],
                                                        time2_dict['Minute'], time1_dict['Second'],
                                                        time2_dict['Second'], no_of_days_passed)

    return ("yes", time_difference)


def convert_current_time(date_dictionary):
    current_time = datetime.now()
    # 2023 - 12 - 25 16: 19:38.993586
    # 2023 12 25 16

    date_dictionary['Year'] = current_time.year
    date_dictionary['Month'] = convert_num(current_time.month)
    date_dictionary['Day'] = convert_num(current_time.day)
    date_dictionary['Hour'] = convert_num(current_time.hour)
    date_dictionary['Minute'] = convert_num(current_time.minute)
    date_dictionary['Second'] = convert_num(current_time.second)


def calculate_time_difference(hour1, hour2, minute1, minute2, second1, second2, days_passed):
    reduce_minute = "no"
    reduce_hour = "no"
    reduce_day = "no"

    # Seconds
    if second2 >= second1:
        seconds_passed = second2 - second1
    else:
        seconds_passed = 60 + second2 - second1
        reduce_minute = "yes"

    # Minutes
    if minute2 >= minute1:
        minutes_passed = minute2 - minute1
    else:
        minutes_passed = 60 + minute2 - minute1
        reduce_hour = "yes"

    if reduce_minute == "yes":
        if minutes_passed == 0:
            minutes_passed = 59
        else:
            minutes_passed -= 1

    # Hours
    if hour2 >= hour1:
        hours_passed = hour2 - hour1
    else:
        hours_passed = 24 + hour2 - hour1
        reduce_day = "yes"

    if reduce_hour == "yes":
        if hours_passed == 0:
            hours_passed = 23
        else:
            hours_passed -= 1

    # Days
    if days_passed > 0:
        if reduce_day == "yes":
            days_passed -= 1
        hours_passed += 24 * days_passed

    hours_passed = str(hours_passed).zfill(2)
    minutes_passed = str(minutes_passed).zfill(2)
    seconds_passed = str(seconds_passed).zfill(2)

    return f"{hours_passed}:{minutes_passed}:{seconds_passed}"

def convert_num(input_num):
    if input_num < 10:
        return f"0{input_num}"
    else:
        return str(input_num)


def print_message_to_output(msg, args_library):
    try:
        # with open(f"{args_library.WorkingDir}{args_library.server_output}", "a") as output_file:
        with open(f"{args_library.server_output}", "a") as output_file:
            output_file.write(f"\n<ul><li>{msg}</li></ul>\n")
    except Exception as e:
        print(f"Failed to open output file: {e}\n")
        sys.exit()


def print_initial_running_progress(args_library):
    args_library.progress_report = "ProgressReport.html"
    args_library.alt_msa_status = "MSA_STATUS.txt"

    with open(args_library.WorkingDir + args_library.alt_msa_status, "a") as ALT_STATUS:
        ALT_STATUS.write("<ul class=\"in_progress\"><li>Generating alternative alignments</li></ul>\n")

    with open(args_library.WorkingDir + args_library.progress_report, "a") as PROGRESS:
        PROGRESS.write("<p><font face=Verdana size=2>\n")

        if args_library.Redirect_From_MAFFT != "1":
            PROGRESS.write("<ul class=\"in_progress\"><li>Generating the base alignment</li></ul>\n")

        PROGRESS.write("<ul class=\"in_progress\"><li>Constructing bootstrap guide-trees</li></ul>\n")
        # PROGRESS.write("REPLACE")
        PROGRESS.write("<ul class=\"in_progress\"><li>Generating alternative alignments</li></ul>\n")

        if args_library.PROGRAM == "GUIDANCE":
            PROGRESS.write("<ul class=\"in_progress\"><li>Calculating GUIDANCE scores</li></ul>\n")

        if args_library.PROGRAM == "HoT":
            PROGRESS.write("<ul class=\"in_progress\"><li>Calculating HoT scores</li></ul>\n")

        if args_library.PROGRAM == "GUIDANCE2":
            PROGRESS.write("<ul class=\"in_progress\"><li>Calculating GUIDANCE2 scores</li></ul>\n")

        if args_library.PROGRAM == "GUIDANCE3":
            PROGRESS.write("<ul class=\"in_progress\"><li>Calculating GUIDANCE3 scores</li></ul>\n")

        PROGRESS.write("</font>\n")

@timeit
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
                # line = re.sub(r"<Calculating [a-zA-Z0-10]+ scores>", message, line)

            # elif "Calculating" in message and "Calculating" in line:
            #     line = re.sub(r"<Calculating \s+ scores for tree # \d+>", message, line)
                # line = line.replace("Calculating GUIDANCE2 scores", message)
                # progress.write(line)
            # elif "Finished calculating" in message and "Calculating" in line:
            #     line = line.replace("in_progress", "finished")
                # line = line.replace("Started", "Finished")
                # line = re.sub(r"<Calculating \s+ scores for tree # \d+>", message, line)
                # progress.write(line)
            else:
                progress.write(line)

