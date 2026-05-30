# python/python-3.6.7

import os
import logging
import shutil
import datetime

if __name__ == '__main__':
        from sys import argv

        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('path',
                            help='A path to a folder in which the sweeps analysis will be written.',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
        parser.add_argument('jobID', type=str,
                            help='The size of a window to which a score will be computed')
                            
        args = parser.parse_args()
        
        results_dir = os.path.join('/bioseq/data/results/Guidance', args.jobID)
        endFile = os.path.join( results_dir, f'GUIDANCE_{args.jobID}.END_OK')
        if os.path.exists(endFile): 
            status = 'PASS'
        else: 
            status = 'FAIL'
        
        results_url = f'http://guidance.tau.ac.il/results/{args.jobID}'
        date = datetime.datetime.today().strftime('%d%m%Y')
        with open(os.path.join(args.path, f'guidance_{date}.txt'), "w") as f:
            f.write(f'{status},{results_url}')
        f.close()
        
        

        
       

