
import SharedConsts as CONSTS
from utils import send_email


run_number = '1111'
email = 'josefspr@gmail.com'
send_email(smtp_server=CONSTS.SMTP_SERVER,
    sender=CONSTS.ADMIN_EMAIL,
    receiver=f'{email}',
    subject=f'{CONSTS.WEBSERVER_NAME.upper()} - your job has been submitted! (Run number: {run_number})',
    content='this is content')
