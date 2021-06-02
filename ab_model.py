import glob, os, sys, fnmatch, getpass, subprocess

pdb = sys.argv[1][0:4].upper()
folderh = sys.argv[1].replace(' ', '_')[0:6].upper()
folderl = sys.argv[1].replace(' ', '_')[0:5].upper() + sys.argv[1][6].upper()
folderall = './RosettaABTest/' + sys.argv[1].replace(' ', '_').upper()
if os.path.isdir(folderall):
	sys.exit()
os.system('mkdir ' + folderall)
os.system('( cd ' + folderall + ' && python ../../clean_pdb.py ' + folderh.replace('_', ' ') + ' )')
os.system('( cd ' + folderall + ' && python ../../clean_pdb.py ' + folderl.replace('_', ' ') + ' )')
os.system('( cd ' + folderall + ' && cat *.fasta > ' + pdb + '.fasta )')
with open(folderall + '/' + pdb + '.fasta', 'r') as file:
	filedata = file.read()

filedata = filedata.replace(folderh, 'heavy').replace(folderl, 'light')

with open(folderall + '/' + pdb + '.fasta', 'w') as file:
	file.write(filedata)

os.system('mkdir ' + folderall + '/grafting')

with open('./abH3.flags', 'r') as file:
	flags = file.read()
with open(folderall + '/abH3.flags', 'w') as file:
	file.write(flags)

os.system('( cd ' + folderall + ' && /home/kodalip/Pranav/RosettaTutorials/Rosetta_AB/main/source/bin/antibody.linuxgccrelease -fasta ' + pdb + '.fasta -antibody::grafting_database /home/kodalip/Pranav/RosettaTutorials/Rosetta_AB/tools/antibody/ -antibody::blastp /dors/meilerlab/apps/Linux2/x86_64/blast/2.7.1/bin/blastp -antibody:n_multi_templates 1 -exclude_pdb ' + pdb.lower() + ' -allow_omega_mismatches_for_north_clusters )')

os.system('mkdir ' + folderall + '/H3_modeling')

with open('./run_H3.slurm', 'r') as file:
	RunH3 = file.read()

RunH3 = RunH3.replace('OPTIONSP', os.path.dirname(os.path.abspath(folderall + '/abH3.flags')))

RunH3 = RunH3.replace('PATHL', pdb)

with open(folderall + '/run_H3.slurm', 'w') as file:
	file.write(RunH3)

HOST = "login.accre.vanderbilt.edu"
RunMatches = ''
root = os.path.dirname(os.path.abspath(folderall + '/abH3.flags'))
RunMatches += ('cd /dors/meilerlab' + root + '\n')
RunMatches += ('sbatch run_H3.slurm' + '\n')
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

