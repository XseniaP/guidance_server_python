import os
import SharedConsts as sc
from GuidanceJobSubmitter import GuidanceJobSubmitter
from Job_Manager_Thread_Safe import Job_Manager_Thread_Safe
from utils import logger


class Job_Manager_Thread_Safe_Guidance:
    def __init__(self, max_number_of_process: int, upload_root_path: str, input_file_name: str, func2update_html_kraken, func2update_html_postprocess):
        self.__input_file_name = input_file_name
        self.__func2update_html_kraken = func2update_html_kraken
        self.__func2update_html_postprocess = func2update_html_postprocess
        self.__search_engine = GuidanceJobSubmitter()
        function2call_processes_changes_state = {
            sc.GUIDANCE_JOB_PREFIX: self.__func2update_html_kraken,
            sc.POSTPROCESS_JOB_PREFIX: self.__func2update_html_postprocess
        }
        function2append_process = {
            sc.GUIDANCE_JOB_PREFIX: self.__guidance_process,
            sc.POSTPROCESS_JOB_PREFIX: self.__postprocess_process
        }
        paths2verify_process_ends = {
            #when the job crashes/ finished this file path will be checked to set the change to finished if file exists of crashed if file doesn't.
            #for a string of: '' it won't set the state
            sc.GUIDANCE_JOB_PREFIX: lambda process_id: os.path.join( upload_root_path, process_id, f'GUIDANCE_{process_id}.END_OK'),
            sc.POSTPROCESS_JOB_PREFIX: lambda process_id: os.path.join(os.path.join(upload_root_path, process_id), sc.FINAL_OUTPUT_FILE_NAME)
        }
        self.__job_manager = Job_Manager_Thread_Safe(max_number_of_process, upload_root_path, function2call_processes_changes_state, function2append_process, paths2verify_process_ends)

    def __guidance_process(self, process_folder_path: str, email_address, job_name):
        logger.info(f'process_folder_path = {process_folder_path}')
        #file2fltr = os.path.join(process_folder_path, self.__input_file_name)
        pbs_id, _ = self.__search_engine.submit_job(process_folder_path, None)
        return pbs_id
    
    def __postprocess_process(self, process_folder_path: str, k_threshold, species_list):
        logger.info(f'process_folder_path = {process_folder_path}')
        pbs_id = '' #run_post_process(process_folder_path, k_threshold, species_list)
        return pbs_id
        
    def __get_state(self, process_id, job_prefix):
        state = self.__job_manager.get_job_state(process_id, job_prefix)
        if state:
            return state
        logger.warning(f'process_id = {process_id}, job_prefix = {job_prefix} not in __job_manager')
        return None
        
    def get_guidance_job_state(self, process_id):
        return self.__get_state(process_id, sc.GUIDANCE_JOB_PREFIX)
        
    def get_postprocess_job_state(self, process_id):
        return self.__get_state(process_id, sc.POSTPROCESS_JOB_PREFIX)
    
    def add_guidance_process(self, process_id: str, email_address, job_name):
        logger.info(f'process_id = {process_id}')
        self.__job_manager.add_process(process_id, sc.GUIDANCE_JOB_PREFIX, email_address, job_name)
    
    def add_postprocess(self, process_id: str, k_threshold, species_list):
        logger.info(f'process_id = {process_id}')
        self.__job_manager.add_process(process_id, sc.POSTPROCESS_JOB_PREFIX, k_threshold, species_list)
    
    def get_job_state(self, process_id: str, job_prefix: str):
        logger.info(f'process_id = {process_id} job_prefix = {job_prefix}')
        return self.__job_manager.get_job_state(process_id, job_prefix)

    def get_processes_state_dict(self):
        return self.__job_manager.get_processes_state_dict()
    
    def get_job_name(self, process_id):
        return self.__job_manager.get_job_name(process_id)
        
    def clean_internal_state(self):
        self.__job_manager.clean_internal_state()
        
    def get_waiting_list(self):
        return self.__job_manager.get_waiting_list()
        
    def add_process_from_waiting_list(self):
        self.__job_manager.add_process_from_waiting_list()