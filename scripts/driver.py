import sys
import os
sys.path.insert(0, '../build/')
from ase.io import read
import numpy as np
import argparse
from subprocess import Popen, PIPE, run
from pathlib import Path

ROOT = os.path.dirname(__file__)
back2root = ['&&','cd',ROOT]

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--rm', action='store_true',
                        help='remove entries in all workspaces')
    parser.add_argument('--init', action='store_true',
                        help='init all signac project subfolders')
    parser.add_argument('--run', action='store_true',
                        help='run all signac project subfolders')
    parser.add_argument('-np','--parallel', type=int, default=-1,
                        help='set the number of process to use when running a project')
    args = parser.parse_args()

    if args.parallel == -1:
        parallel = ''
    else:
        parallel = '--parallel {}'.format(args.parallel)

    if args.rm:
        print('Remove files from workspaces')
        run('rm -rf ./*/workspace/*',shell=True)

    if args.init:
        print('Initialize all projects')
        for path in Path('./').rglob('init.py'):
            print('    ',os.path.dirname(path))
            move = ['cd',os.path.dirname(path),'&&']
            command = ' '.join(move+['python', os.path.abspath(path)]+back2root)
            run(command, shell=True)

    if args.run:
        print('Run all project')
        for path in Path('./').rglob('project.py'):
            print('    ',os.path.dirname(path))
            move = ['cd',os.path.dirname(path),'&&']
            command = ' '.join(move+['python', os.path.abspath(path), 'run', parallel, '--progress']+back2root)
            run(command ,shell=True)