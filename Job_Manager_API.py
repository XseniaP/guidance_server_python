import os
import shutil
import uuid
import json
from time import time
from random import randint
import pandas as pd
from InputValidator import InputValidator
from Job_Manager_Thread_Safe_Guidance import Job_Manager_Thread_Safe_Guidance
from utils import send_email, State, logger, LOGGER_LEVEL_JOB_MANAGE_API, get_job_logger
from SharedConsts import EMAIL_CONSTS
import SharedConsts as CONSTS
logger.setLevel(LOGGER_LEVEL_JOB_MANAGE_API)


class Job_Manager_API:
    def __init__(self, max_number_of_process: int, upload_root_path: str, input_file_name: str, func2update_html):
        self.__input_file_name = input_file_name
        self.__upload_root_path = upload_root_path
        self.__j_manager = Job_Manager_Thread_Safe_Guidance(max_number_of_process, upload_root_path, input_file_name,
                                                             self.__process_state_changed, self.__process_state_changed)
        self.input_validator = InputValidator()
        self.__func2update_html = func2update_html

    def __build_and_send_mail(self, process_id, subject, content, email_address):
        job_logger = get_job_logger(process_id)
        try:
            '''
            send_email('mxout.tau.ac.il', 'TAU Evolseq <evolseq@tauex.tau.ac.il>',
                       email_address, subject=subject,
                       content= content)
            '''
            send_email(CONSTS.SMTP_SERVER, CONSTS.ADMIN_EMAIL,
                       email_address, subject=subject,
                       content= content)
            log_msg = f'sent email to {email_address}'
            logger.info(log_msg)
            if job_logger: 
                job_logger.info(log_msg)
        except:
            logger.exception(f'failed to sent email to {email_address}')
            if job_logger: 
                job_logger.exception(f'failed to sent email to {email_address}')

    def __process_state_changed(self, process_id, state, email_address, job_name, job_prefix):
        log_msg = f'called __process_state_changed state={state.str()}'
        logger.info(log_msg)
        if job_logger: 
            job_logger.info(log_msg)
        if email_address != None and job_name != 'daily_test':
            if state == State.Finished:
               self.__build_and_send_mail(process_id, EMAIL_CONSTS.FINISHED_TITLE, EMAIL_CONSTS.FINISHED_CONTENT.format(results_url=CONSTS.WEBSERVER_RESULTS_URL_EXT, process_id=process_id), email_address)
               
            elif state == State.Crashed:
                self.__build_and_send_mail(process_id, EMAIL_CONSTS.CRASHED_TITLE, EMAIL_CONSTS.CRASHED_CONTENT.format(results_url=CONSTS.WEBSERVER_RESULTS_URL_EXT, process_id=process_id), email_address)
        else:
            log_msg = f'process_id = {process_id} email_address is None, state = {state}, job_name = {job_name}'
            logger.warning(log_msg)
            job_logger = get_job_logger(process_id)
            if job_logger:
                job_logger.info(log_msg)
                
        self.__func2update_html(process_id, state)

    def __validate_email_address(self, email_address):
        if len(email_address) > 100:
            return False
        if '@' in email_address and '.' in email_address:
            return True
        return False

    def get_new_process_id(self):
        #return str(uuid.uuid4()) commented out by JS
        time_str = str(round(time()))
        rand_str = str(randint(1000,9999))
        return f'{time_str}{rand_str}'

    def add_guidance_process(self, process_id: str, email_address: str, job_name: str):
        job_logger = get_job_logger(process_id)
        log_msg = f'add_guidance_process: process_id = {process_id} email_address = {email_address}'
        logger.info(log_msg)
        if job_logger:
                job_logger.info(log_msg)
        if email_address:
            is_valid_email = self.__validate_email_address(email_address)
        else:
            is_valid_email = True
        if is_valid_email:
            logger.info(f'validated email address')

            self.__j_manager.add_guidance_process(process_id, email_address, job_name)

            if email_address and job_name != 'daily_test':
                self.__build_and_send_mail(process_id, EMAIL_CONSTS.INIT_TITLE, EMAIL_CONSTS.INIT_CONTENT.format(results_url=CONSTS.WEBSERVER_PROCESS_STATE_URL_EXT, process_id=process_id), email_address)
        
            return True
        logger.warning(f'process_id = {process_id}, can\'t add process: is_valid_email = {is_valid_email}')
        if job_logger:
            job_logger.warning(log_msg)
        return False
        
    def add_postprocess(self, process_id: str, species_list: list, k_threshold: float):
        parent_folder = os.path.join(self.__upload_root_path, process_id)
        if os.path.isdir(parent_folder):
            self.__j_manager.add_postprocess(process_id, k_threshold, species_list)
            return True
        logger.warning(f'process_id = {process_id} don\'t have a folder')
        return None
    
    def get_guidance_job_state(self, process_id):
        #return self.__j_manager.get_guidance_job_state(process_id)
        resultsDir = os.path.join(CONSTS.WEBSERVER_RESULTS_DIR, process_id)
        if not os.path.exists(resultsDir):
           return None
        elif not os.path.exists(os.path.join(resultsDir, 'qsub.cmds')):
           return State.Init
        elif os.path.exists(os.path.join( CONSTS.WEBSERVER_RESULTS_DIR, process_id, f'GUIDANCE_{process_id}.END_FAIL')):
           return State.Crashed
        elif not os.path.exists(os.path.join( CONSTS.WEBSERVER_RESULTS_DIR, process_id, f'GUIDANCE_{process_id}.END_OK')):
           return State.Running
        else:
           return State.Finished
        
    def get_processes_state_dict(self):
        return self.__j_manager.get_processes_state_dict()
        
    #def get_waiting_list(self):
    #   return self.__j_manager.get_waiting_list()
    
    def clean_internal_state(self):
        self.__j_manager.clean_internal_state()
        
    def get_waiting_list(self):
        return self.__j_manager.get_waiting_list()
        
    def add_process_from_waiting_list(self):
        self.__j_manager.add_process_from_waiting_list()
       
   
