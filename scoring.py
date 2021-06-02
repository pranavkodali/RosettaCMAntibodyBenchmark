import os, fnmatch, pickle, operator, math, subprocess

q=[]

for dirname in os.listdir(os.getcwd()):
	if len(dirname) == 7 and os.path.isdir(os.path.join(os.getcwd(), dirname)):
		if dirname[4] == '_':
			if os.path.isdir(os.path.join(os.path.join(os.getcwd(), dirname), 'ScorePlots')) is not True and (os.path.isfile(os.path.join(os.path.join(os.getcwd(), dirname), 'all/relax99_S_0010_0001.pdb.gz')) is True or os.path.isfile(os.path.join(os.path.join(os.getcwd(), dirname), 'all/relax99_S_0007_0001.pdb.gz')) is True):
				print(dirname)
				p = subprocess.Popen('python runrmsd.py ' + dirname + ' 1', shell=True)
				#p.wait()
				q.append(p)
	if len(q) > 10:
		exitcodes = [p.wait() for p in q]
		q = []

exitcodes = [p.wait() for p in q]
