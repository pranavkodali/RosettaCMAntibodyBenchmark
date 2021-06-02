#!/bin/env python

import glob, os, sys, urllib, urllib.request, urllib.error, urllib.parse, webbrowser, requests, gzip, shutil
from shutil import copyfile

folder = sys.argv[1]
folderall = sys.argv[2]

if not os.listdir('./' + folder):
	pass
else:

		if os.path.exists('./' + folder +'/all.fasta'):
			pass
		else:
			os.system('( cd ' + folder + ' && cat *.fasta > all.fasta )')

		
		url = 'http://www.genome.jp/tools-bin/clustalw/html'
		with open('./' + folder +'/all.fasta', 'rb') as allFile: 
			uploadText = allFile.read()
		data = {'pwalignment': 'fast','type': 'protein','sequence': uploadText,'gapopen': '10','ktuple': '1','window': '5','pairgap': '3','topdiags': '5','score': 'percent','gapext': '0.05','pwgapopen': '10.0','pwgapext': '0.1','pwmatrix': 'blosum','matrix': 'blosum'}
		r = requests.post('https://www.genome.jp/tools-bin/clustalw', data = data)
		with open('./' + folder +'/clustalwresults.html', 'w') as f:
                        f.write(r.text)



		f = open(glob.glob('./' + folder +'/clustalwresults.html')[0], 'r')
		text = f.read()

		start = 0
		start = text.find('Aligning', start)
		end = text.find('Guide', start)
		while text.find('Sequences (', start) < end and text.find('Sequences (', start) > -1:  
			idx = text.find('Sequences (', start)
			clsidx = text.find(':', idx)
			one = int(text[idx+11:clsidx])
			idx = clsidx
			clsidx = text.find(')', idx)
			two = int(text[idx+1:clsidx])
			idx = text.find(':', clsidx)
			clsidx = text.find('S', idx)
			if clsidx - idx > 10:
				clsidx = text.find('G', idx)
			score = float(text[idx+2:clsidx-1])
			start = clsidx - 1

		struct = {(a, b, c): 0 for a in range(one + 1) for b in range(two) for c in range(0)}
		maxentries = one
		if maxentries > 5:
			start = 0
			start = text.find('Aligning', start)
			end = text.find('Guide', start)
			while text.find('Sequences (', start) < end and text.find('Sequences (', start) > -1: 
				idx = text.find('Sequences (', start)
				clsidx = text.find(':', idx)
				one = int(text[idx+11:clsidx])
				idx = clsidx
				clsidx = text.find(')', idx)
				two = int(text[idx+1:clsidx])
				idx = text.find(':', clsidx)
				clsidx = text.find('S', idx)
				if clsidx - idx > 10:
					clsidx = text.find('G', idx)
				score = float(text[idx+2:clsidx-1])
				struct[one, two, 0] = score
				struct[two, one, 0] = score
				start = clsidx - 1

			minimum = 100
			for x in range(1,maxentries):
				for y in range(x + 1, maxentries + 1):
					if struct[x, y, 0] < minimum:
						minimum = struct[x, y, 0]
						minx = x
						miny = y

			minimumcombo = 100000000000000000000;

			for p in range(1,maxentries - 1):
				if p == minx or p == miny:
					pass
				else:
					for q in range(p+1,maxentries):
						if q == minx or q == miny or q == p:
							pass
						else:
							for r in range(q+1, maxentries):
								if r == minx or r==miny or r==q or r==p:
									pass
								else:
									test = struct[one, two, 0] * struct[one, p, 0] * struct[one, q, 0] * struct[one, r, 0] * struct[two, p, 0] * struct[two, q, 0] * struct[two, r, 0] * struct[p, q, 0] * struct[p, r, 0] * struct[q, r, 0]
									if test < minimumcombo:
										minimumcombo = test
										minp = p
										minq = q
										minr = r

			output = 'ideal: ' + repr(minx) + ',' + repr(miny) + ',' + repr(minp) + ',' + repr(minq) + ',' + repr(minr) + ','

			start = 0
			start = text.find('Pearson', start)
			end = text.find('Start of', start)
			while text.find('Sequence', start) < end and text.find('Sequence', start) > -1:  
				idx = text.find('Sequence', start)
				clsidx = text.find(':', idx )
				numberfound = int(text[idx+9:clsidx])
				if numberfound == minx:
					minxname = text[clsidx+2:clsidx+8]
				elif numberfound == miny:
					minyname = text[clsidx+2:clsidx+8]
				elif numberfound == minp:
					minpname = text[clsidx+2:clsidx+8]
				elif numberfound == minq:
					minqname = text[clsidx+2:clsidx+8]
				elif numberfound == minr:
					minrname = text[clsidx+2:clsidx+8]
				start = clsidx + 1
			fileList = os.listdir('./' + folder)
			for i in [minxname, minyname, minpname, minqname, minrname]:
				with gzip.open('./templates/' + i[0:4].lower() + '.pdb.gz', 'rb') as f_in:
					with open('./' + folderall +'/all/' + i[0:4].lower() + '.pdb', 'wb') as f_out:
						shutil.copyfileobj(f_in, f_out)
			with gzip.open('./templates/' + minxname[0:4].lower() + '.pdb.gz', 'rb') as f_in:
				with open('./' + folderall +'/1/' + minxname[0:4].lower() + '.pdb', 'wb') as f_out:
					shutil.copyfileobj(f_in, f_out)
			with gzip.open('./templates/' + minyname[0:4].lower() + '.pdb.gz', 'rb') as f_in:
				with open('./' + folderall +'/2/' + minyname[0:4].lower() + '.pdb', 'wb') as f_out:
					shutil.copyfileobj(f_in, f_out)
			with gzip.open('./templates/' + minpname[0:4].lower() + '.pdb.gz', 'rb') as f_in:
				with open('./' + folderall +'/3/' + minpname[0:4].lower() + '.pdb', 'wb') as f_out:
					shutil.copyfileobj(f_in, f_out)
			with gzip.open('./templates/' + minqname[0:4].lower() + '.pdb.gz', 'rb') as f_in:
				with open('./' + folderall +'/4/' + minqname[0:4].lower() + '.pdb', 'wb') as f_out:
					shutil.copyfileobj(f_in, f_out)
			with gzip.open('./templates/' + minrname[0:4].lower() + '.pdb.gz', 'rb') as f_in:
				with open('./' + folderall +'/5/' + minrname[0:4].lower() + '.pdb', 'wb') as f_out:
					shutil.copyfileobj(f_in, f_out)
		else:
			os.remove('./' + folder + '/' + 'all.fasta')
			os.remove('./' + folder + '/' + 'results.html')


