#!/bin/env python

import glob, os, sys, fnmatch, getpass, subprocess
from shutil import copyfile

folder = sys.argv[1]

path = os.path.dirname(os.path.realpath(sys.argv[0]))

if len(sys.argv) > 2:
	HOST = sys.argv[2]
else:
	HOST = "bluefin"

RunMatches = ''
#for root, dirnames, filenames in os.walk(path + '/' + folder):
#	for filename in fnmatch.filter(filenames, '4fqi.pdb'):
#		os.system('python ClustaloGrishin.py ' +  root + '/alignment.aln ' + root[34:])


for root, dirnames, filenames in os.walk(path + '/' + folder):
	for filename in fnmatch.filter(filenames, 'run_relax.slurm'):
		
		#copyfile('./relaxscore.txt', root + '/relaxscore.txt')
		copyfile('./sc_parser.bash', root + '/sc_parser.bash')
		RunMatches += ('cd /dors/meilerlab' + root + '\n')
		RunMatches += ('sbatch run_relax.slurm' + '\n')


#HOST="bluefin"
#HOST="login.accre.vanderbilt.edu"
# Ports are handled in ~/.ssh/config since we use OpenSSH
COMMAND=RunMatches

ssh = subprocess.Popen(["ssh", "%s" % HOST, COMMAND],
                       shell=False,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE)
result = ssh.stdout.readlines()
if result == []:
    error = ssh.stderr.readlines()
    print("ERROR: %s" % error, file=sys.stderr)
else:
    print(result)


