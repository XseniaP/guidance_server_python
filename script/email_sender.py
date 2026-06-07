import logging
import os
logger = logging.getLogger('main') # use logger instead of printing

def send_email(smtp_server, sender, receiver, subject='', content=''):
    from email.mime.text import MIMEText
    from smtplib import SMTP
    try:
        msg = MIMEText(content)
        msg['Subject'] = subject
        msg['From'] = sender
        msg['To'] = receiver
        s = SMTP(smtp_server)
        ans = s.send_message(msg)
        s.quit()
        return 'ok'
    except:
        return 'error'
    

if __name__ == '__main__':
    # This block will be executed only when you run it as your main program.
    # If this module is being imported from another script, this block won't be executed, however the function will be available...
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('smtp_server', help='SMTP server')
    parser.add_argument('sender', help='From whom the email will be sent')
    parser.add_argument('receiver', help='To whom the email will be sent')
    parser.add_argument('--subject', help='The subject of the email')
    parser.add_argument('--content', help='The content of the email', default='')
    parser.add_argument('--run_number', help='Job run number', default='')
    args = parser.parse_args()

    error = False
    if args.run_number != '':
        error_path = os.path.join('/lsweb/josef_results/effectidor', args.run_number, 'error.txt')
        if os.path.exists(error_path): 
            with open (error_path, "r") as f:
                for line in f:
                    if line.find('returned non-zero exit status 1') != -1:
                        error = True
                    elif line.find('error:') != -1:
                        error = True
            f.close()
    if not error:
        ans = send_email(args.smtp_server, args.sender, args.receiver, args.subject, args.content)
    else:
        ans = 'error in run'
    print (ans, end="")
    
# python email_sender.py mxout.tau.ac.il evolseq@tauex.tau.ac.il josefspr@gmail.com --subject test --content testing  