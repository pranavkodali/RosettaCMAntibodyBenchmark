import os, sys, gzip, re, requests, csv, tsne
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import scipy.cluster.hierarchy as sch
import plotly as py
import plotly.graph_objs as go
from Bio.PDB import *
from Bio.Alphabet import *
from Bio.Seq import *


directory = os.path.dirname(os.path.abspath(__file__)) + '/templates'

pdb_list = ['1cz8','4ttd','2qsc','4ydl','2qqk','4ydi','4ydk','5jz7','6b0g','5xhv','3grw','4tsa','4xnq','3go1','1vge','5ig7','2a9m','3idg','5i19','5cgy','4kq3','5awn','5y9k','3x3f','3n9g','4jdv','5whj','4olz','5nb5','4olx','3auv','5kn5','1n0x','4olv','4olu','5igx','5vl7','2xzc','5veb','3bn9','5f96','2fl5','6bhz','5f3h','4j6r','5ukp','3l5x','5uko','3fn0','4jn2','3r1g','4hpo','3k2u','4hf5','3eot','5e8e','2fjh','5uby','5ubz','6ayz','2r56','5kw9','4jy5','4mxv','3h42','5tfw','5ogi','3hae','4n0y','4hfw','4hfu','3eo0','3eo9','5v7r','5j13','4nwt','5gru','5bo1','4xtr','3juy','5tf1','4ywg','4hs8','5bk1','4nrx','4nrz','5wuv','1igm','5vob','5jo4','4jzn','4jzo','4jzj','5xku','1ce1','4r96','5ik3','4dke','5ur8','4ers','3kdm','4qci','5tfs','4py7','3eyf','6ehy','6ehx','1w72','3skj','1za3','5itb','3n85','3c09','2v7n','5wcd','5jof','5wca','4hpy','2fb4','3qrg','5gks','3l95','2yk1','4dkf','3ghe','4s1r','1b2w','2d7t','5ea0','4nki','5u0u','5u0r','5tgb','4jam','3kr3','3na9','4zd3','5uea','5f6h','4eow','5ifh','5vvf','4hk0','3tnm','5dd5','4lsu','4lss','5ucb','5t4z','5bk2','4k94','5bk5','1jps','5fcu','3nab','3nac','3ma9','4g80','4hk3','5uek','4cni','4od2','3u6r','5iwl','3piq','4uu9','6be3','5n4j','5xxy','4lex','3u7y','4nyl','8fab','3mly','4zs6','5jrp','4leo','6erx','3se9','4ioi','4k9e','1bvk','5k9j','5sy8','4llu','4lly','4odx','3pnw','2xra','3uls','4qf1','5f6i','5mvz','5te7','4qhl','4jfx','4imk','5l6y','5jr1','3oaz','2uzi','4qhu','4buh','3bdy','1gc1','5lbs','3ux9','3hr5','3p0y','5gzo','4ky1','4gsd','4yjz','1iqd','4hg4','4oqt','4jlr','2vxv','4edw','3kym','4gxv','6b3k','5bjz','6b3d','5k59','5x8m','4lst','3hi6','4liq','3dvg','5fuz','6bp2','3dvn','2eiz','2g75','5ukq','4lkx','5n2k','5bk3','4fze','4g7v','2qr0','2r0k','5d72','5vic','4o5l','4n9g','4y5x','5b71','3naa','3b2u','4hha','5u4r','4y5v','4ygv','5ukn','4hcr','4v1d','3pp4','5xaj','2hwz','2aj3','6b0w','4rrp','5feh','5f9o','6b0e','5ggs','6b0a','5drz','5drw','3ztn','3g6a','5cex','5gmq','5f9w','4ps4','4jpk','4jpi','1nl0','4f57','4jha','3u30','3hc3','6b08','4f58','2jb5','4nzu','3u4b','4irz','4z0x','4ojf','4ud3','4lmq','3nfs','4ogy','4uta','5w6g','4zfg','5w6c','5bvj','3efd','4xc1','5u15','5ibu','5tru','5udc','3p30','5trp','1l7i','4xgz','4wuk','4s1s','3hmx','3giz','5dd0','2fjf','1rhh','4mwf','4ut9','1dql','5ggq','5ggt','4ris','5ggv','4ptu','5u68','4fqq','4yk4','4nik','4fqi','4jfz','5d7s','3lmj','1wt5','4fqc','3qeh','2vxq','3sob','4hie','5n4g','5kg9','5c6w','4tsb','5vig','5czx','5czv','3hc4','3hc0','2wuc','4h8w','5ibt','2b2x','5tdn','4k3j','3gkw','3mlw','4zyk','3sqo','4dn3','5eu7','6axl','5alb','2qqn','2b1h','6axk','4xml','3u0t','5bv7','3mme','5waw','4qxg','5chn','4u6v','2xwt','2f5a','2nyy','3gjf','5hi4','4p59','5t29','4kmt','1u6a','5bzd','4i18','4g6f','5ush','2ny1','3idx','2r8s','4xvj','5tdo','3aaz','5g64','5mes','4evn','5grz','5aam','5grx','5gry','3wd5','5grv','5grw','4xvt','5lsp','2a9n','4npy','1nfd','4jqi','6c6z','5usi','1dfb','4nm4','5d1z','3zl4','1tzi','1tzh','4np4','4m6o','4m6n','5vf6','6az2','3nh7','5sx5','5usl','4om0','4om1','3uji','5d6c','4j4p','1mhp','4ypg','5tzt','5ob5','3tcl','3m8o','2yss','5d9q','3s34','5yax','6b14','3fzu','4rx4','2agj','5wym','3sdy','5fha','5fhb','4olw','4ot1','4rav','5vsi','5y11','3so3','4tsc','5tz2','5t6l','4dgy','4qhm','1fve','1fvd','4m62','5dum','6b0h','3i9g','5dur','3ncj','5nhw','4m1d','3inu','5i1l','5i1i','5i1h','5i1e','1aqk','4oaw','5i1a','5i1c','4hwe','1rz7','4z5r','4d9q','5i8c','5it2','2hfg','5c2b','5cin','4iml','5cil','2cmr','2ghw','3pgf','5i15','5jw5','5i17','5i16','3lzf','4lf3','5wk2','3oau','3qhf','5v7u','1s3k','4nnp','4jg1','5esv','5w08','4x7s','5w05','5w06','5u7o','4g5z','5u3p','2xtj','1dn0','4r7d','4xxd','2j6e','4yhy','4yhz','3g04','4jkp','5t5b','2r0l','4yhl','4yho','1hez','1g9m','4xx1','4nug','4i77','3gbn','3lh2','3qwo','1bj1','1dee','5gs1','5vag','4ut7','5cck','5iq9','3nps','5j75','5iq7','2nz9','5uem','4m5y','5uoe','5bzw','5anm','4lkc','5i1d','6ehw','5ty6','5i1g','4fql','2eks','5f72','3ujj','5cd3','4wuu','4x4x','4x4y','3w9e','5tpl','3uc0','4q2z','4xak','5u3j','5c8j','5u3o','5u3l','3mac','5dr5','4uv7','4uv4','4d9l','4r8w','3qcu','6azx','3qot','6azz','2h9g','6azm','5cus','5eii','2yc1','4hs6']

errors = 0
used = 0

full_seq = {}
heavy_seq = {}
light_seq = {}
L1_seq = {}
L2_seq = {}
L3_seq = {}
H1_seq = {}
H2_seq = {}
H3_seq = {}
FR_L1_seq = {}
FR_L2_seq = {}
FR_L3_seq = {}
FR_L4_seq = {}
FR_H1_seq = {}
FR_H2_seq = {}
FR_H3_seq = {}
FR_H4_seq = {}

start_tag = 'Aligning...'
end_tag = 'Guide tree file created:'

def plot_corr(df, name, size=10):
	corr = df.corr
	print(corr)
	fig, ax = plt.subplots(figsize = (size, size))
	#print(ax.shape)
	cax = ax.matshow(df, cmap = 'RdYlGn')
	plt.xticks(range(len(df.columns)), df.columns, rotation = 90)
	plt.yticks(range(len(df.columns)), df.columns)
	cbar = fig.colorbar(cax, ticks=[-1, 0, 1], aspect = 40, shrink = 0.8)
	plt.savefig('./cluster/' + name + '_graph.pdf')


def read_fasta_file(file_name):
	""" return array of sequences from FASTA formatted file

	file_name -- relative or absolute path to FASTA file
	"""
	seqArray=[]
	seqArrayNames=[]
	seq=""
	for l in open(file_name):
		if l.startswith('>'):
			if not "" == seq:
				seqArray.append(seq)
				seq=""
			seqArrayNames.append(l.rstrip())
		else:
			seq=seq+l.rstrip()
	seqArray.append(seq)
	return seqArray

def find_between(s, first, last):
	start = s.index(first) + len(first)
	end = s.index(last, start)
	return int(s[start:end])

def find_after(s, first):
	if 'Not Aligned' in s:
		return 0
	else:
		return int(round(float(s.split(first, 1)[1].split('\n', 1)[0])))

def write_fasta_file(name, dicts):
	""" writes fasta file

	file_name -- name of file, will also be header / sequence name
	data -- the single sequence to be written
	prefix -- path prepended to the filename, explicitly have a directory separator if needed
	"""
	with open('./cluster/' + name + '.fasta', 'w') as f:
		for keys, data in dicts.items():
			f.write('>%s\n' % keys)
			f.write(data)
			f.write('\n\n')
	with open('./cluster/' + name + '.fasta', 'r') as f:
		uploadText = f.read()
	'''url = 'http://www.genome.jp/tools-bin/clustalw/html'
	data = {'pwalignment': 'fast','type': 'protein','sequence': uploadText,'gapopen': '10','ktuple': '1','window': '5','pairgap': '3','topdiags': '5','score': 'percent','gapext': '0.05','pwgapopen': '10.0','pwgapext': '0.1','pwmatrix': 'blosum','matrix': 'blosum'}
	r = requests.post('https://www.genome.jp/tools-bin/clustalw', data = data)
	with open('./cluster/' + name + '.html', 'w') as f:
                f.write(r.text)'''
	start_tag_found = False
	end_tag_found = False
	with open('./cluster/' + name + '.html', 'r') as in_file:
		with open('./cluster/' + name + '.txt', 'w') as out_file:
			for line in in_file:
				if not start_tag_found:
					if line.strip() == start_tag:
						start_tag_found = True
				else:
					if end_tag in line:
						end_tag_found = True
					if end_tag_found == False:
						if len(line) > 1:
							out_file.write(line)
	data = np.zeros((used, used))
	#data = [[0 for x in range(used)] for y in range(used)]
	with open('./cluster/' + name + '.txt', 'r') as f:
		for line in f:
			one = find_between(line, 'Sequences (', ':') - 1
			two = find_between(line, ':', ')') - 1
			score = find_after(line, 'Score:')
			#print('one ' + str(one))
			#print(two)
			data[one][two] = score
			data[two][one] = score
			data[one][one] = 100
			data[two][two] = 100
	with open('./cluster/' + name + '.csv', 'wb') as f:
		writer = csv.writer(f)
		writer.writerows(data)
	#Y = tsne.tsne(data, 2, 10)
	#trace = go.Scatter(x = Y[:, 0], y = Y[:, 1], mode = 'markers')
	#pylab.scatter(Y[:, 0], Y[:, 1], 20, labels)
	#pylab.show()
	#data.sort()
	#for array in data:
		#array.sort()
	#print(data)
	#data = data[np.mean(data,axis=1).argsort()]
	#trace = go.Heatmap(z = data[:, np.argmax(data, axis=1)].tolist())
	#trace = go.Heatmap(z = np.array(sorted(data, cmp=lambda x, y: list(x).index(100) - list(y).index(100))).tolist())
	#trace = go.Heatmap(z = data.tolist())
	#graphdata = [trace]
	#py.offline.plot(graphdata, filename= os.path.dirname(os.path.abspath(__file__)) + '/cluster/' + name + '_graph.html', auto_open = False)
	df = pd.DataFrame(np.array(data).transpose())
	X = df.values
	d = sch.distance.pdist(X)
	L = sch.linkage(d, method='complete')
	ind = sch.fcluster(L, 0.5 * d.max(), 'distance')
	columns = [df.columns.tolist()[i] for i in list((np.argsort(ind)))]
	df = df.reindex(columns, axis = 1)
	'''X = df.values
	d = sch.distance.pdist(X)
	L = sch.linkage(d, method='complete')
	ind = sch.fcluster(L, 0.5 * d.max(), 'distance')
	columns = [df.columns.tolist()[i] for i in list((np.argsort(ind)))]
	df = df.reindex(columns, axis = 0)'''
	#print(df)
	plot_corr(df, name, size=100)
					
			


def IdentifyCDRs(light_chain, heavy_chain):
	""" Determination of loops in both light and heavy chains

	light_chain -- string representation of light chain
	heavy_chain -- string representation of heavy chain
	"""

	L1_pattern=r'C[A-Z]{1,17}(WYL|WLQ|WFQ|WYQ|WYH|WVQ|WVR|WWQ|WVK|WYR|WLL|WFL|WVF|WIQ|WYR|WNQ|WHL|WHQ|WYM|WYY)'
	L3_pattern=r'C[A-Z]{1,15}(L|F|V|S)G[A-Z](G|Y)'
	H1_pattern=r'C[A-Z]{1,16}(W)(I|V|F|Y|A|M|L|N|G)(R|K|Q|V|N|C|G)(Q|K|H|E|L|R)' # Jeff's mod for ATHM set
   #H1_pattern=r'C[A-Z]{1,16}(W)(I|V|F|Y|A|M|L|N|G)(R|K|Q|V|N|C)(Q|K|H|E|L|R)'
	H3_pattern=r'C[A-Z]{1,33}(W)(G|A|C)[A-Z]{1,2}(Q|S|G|R)'

	''' Identift CDR region and return them as dict with keys: 'FR_H1', 'FR_H2', 'FR_H3', 'FR_H4', 'FR_L1', 'FR_L2', 'FR_L3', 'FR_L4', 'H1', 'H2', 'H3', 'L1', 'L2', 'L3'
	'''
	light_first = light_chain[:65] if len(light_chain) > 130 else light_chain[:60]
	heavy_first = heavy_chain[:70] if len(heavy_chain) > 140 else heavy_chain[:60]

	# LL II GG HH TT

	## L1

	len_FR_L1=0
	res = re.search(L1_pattern,light_first)
	FR_L1 = L1 = False
	if res:
		L1 = res.group()[1:-3]
		L1_start = light_chain.index(L1)
		L1_end = L1_start + len(L1) - 1
		print(("L1 detected: %s (%d residues at positions %d to %d)" % (L1, len(L1), L1_start, L1_end)))
		FR_L1 = light_chain[:L1_start]
		if len(FR_L1) >  24:
			len_FR_L1 = len(FR_L1) - 24
			FR_L1 = light_chain[len_FR_L1:L1_start]
	else:
		print("L1 detected: False")
		raise ValueError()

	light_second = light_chain[L1_end+16+7:L1_end+16+7+80] if len(light_chain) > 130 else light_chain[L1_end+16+7:]

	## L3

	L3 = False
	res = re.search(L3_pattern,light_second)
	if res:
		L3 = res.group()[1:-4]
		L3_start = light_chain.index(L3)
		L3_end = L3_start + len(L3) - 1
		print(("L3 detected: %s ( %d residues at positions %d to %d)" % (L3, len(L3), L3_start, L3_end)))
	else:
		print("L3 detected: False")
		raise ValueError()

	if L1 and L3:
		#L1_start = light_chain.index(L1)
		#L1_end = L1_start + len(L1) - 1

		L2_start = L1_end + 16
		L2_end = L2_start + 7 - 1

		L2 = light_chain[L2_start:L2_start+7]  # L2 is identified here. Current implementation can deal with only 7-resiue L2
		print(("L2 detected: %s (%d residues at positions %d to %d)" % (L2, len(L2), L2_start, L2_end)))

		#FR_L1 = light_chain[:L1_start]
		FR_L2 = light_chain[ L1_end + 1  :  L1_end + 1 + 15					]
		FR_L3 = light_chain[ L2_end + 1  :  L2_end + 1 + L3_start - L2_end - 1 ]
		FR_L4 = light_chain[ L3_end + 1  :  L3_end + 1 + 12					]

		print(("FR_L1: ", FR_L1))
		print(("FR_L2: ", FR_L2))
		print(("FR_L3: ", FR_L3))
		print(("FR_L4: ", FR_L4))
		print(("L segments: ",FR_L1,L1,FR_L2,L2,FR_L3,L3,FR_L4))

		# Light chain sub-type classification. This is useful in the future. But currently this is not used.
		# ... skipped, see Google doc for details
		# FR classification by AHo. This might be useful in the future. But currently this is not used.
		# ... skipped, see Google doc for details


	# HH EE AA VV YY

	## H1
	res = re.search(H1_pattern, heavy_first)
	H1 = False
	len_FR_H1 = 0
	if res:
		H1 = res.group()[4:-4]
		H1_start = heavy_chain.index(H1)
		H1_end = H1_start + len(H1) - 1
		print(("H1 detected: %s (%d residues at positions %d to %d)" % (H1, len(H1), H1_start, H1_end)))
		FR_H1 = heavy_chain[:H1_start]
		if len(FR_H1) >  25:
			len_FR_H1 = len(FR_H1) - 25
			FR_H1 = heavy_chain[len_FR_H1:H1_start]
	else:
		print("H1 detected: False")
		raise ValueError()


	heavy_second = heavy_chain[H1_end+33+15:H1_end+33+15+95+len_FR_H1] if len(heavy_chain) > 140 else heavy_chain[H1_end+33+15:]

	## H3
	H3 = False #H3_and_stem=False
	res = re.search(H3_pattern,heavy_second)
	if res:
		H3 = res.group()[3:-4] #H3_and_stem = res.group()[0:-4]
		H3_start = heavy_chain.index(H3)
		H3_end = H3_start + len(H3) - 1
		print(("H3 detected: %s (%d residues at positions %d to %d)" % (H3, len(H3), H3_start, H3_end)))
	else:
		print("H3 detected: False")
		raise ValueError()


	if H1 and H3:
		#H1_start = heavy_chain.index(H1)
		#H1_end = H1_start + len(H1) - 1
		H2_start = H1_end + 15
		H2_end = H3_start - 33
		H2 = heavy_chain[H2_start:H2_start + H2_end-H2_start+1]
		print(("H2 detected: %s (%d residues at positions %d to %d)" % (H2, len(H2), H2_start, H2_end)))

		#FR_H1 = heavy_chain[:H1_start]

		#if len(FR_H1) >  26:
		#	FR_H1 = light_chain[20:H1_start]

		FR_H2 = heavy_chain[H1_end + 1: H1_end + 1 + H2_start - H1_end - 1]
		FR_H3 = heavy_chain[H2_end + 1: H2_end + 1 + H3_start - H2_end - 1]
		FR_H4 = heavy_chain[H3_end + 1: H3_end + 1 + 12]

		print(("FR_H1: ", FR_H1))
		print(("FR_H2: ", FR_H2))
		print(("FR_H3: ", FR_H3))
		print(("FR_H4: ", FR_H4))
		print(("H segments: ",FR_H1,H1,FR_H2,H2,FR_H3,H3,FR_H4))

	if not (L1 and L3 and H1 and H3):
		if not L1: print(('ERROR: CDR L1 cannot be recognized !!!  L1 pattern: %s' % L1_pattern)) # C[A-Z]{1,17}(WYL|WLQ|WFQ|WYQ|WYH|WVQ|WVR|WWQ|WVK|WYR|WLL|WFL|WVF|WIQ|WYR|WNQ|WHL|WHQ|WYM|WYY)'
		if not L3: print(('ERROR: CDR L3 cannot be recognized !!!  L3 pattern: %s' % L3_pattern)) # C[A-Z]{1,15}(L|F|V|S)G[A-Z](G|Y)
		if not H1: print(('ERROR: CDR H1 cannot be recognized !!!  H1 pattern: %s' % H1_pattern)) # C[A-Z]{1,16}(W)(I|V|F|Y|A|M|L|N|G)(R|K|Q|V|N|C)(Q|K|H|E|L|R)
		if not H3: print(('ERROR: CDR H3 cannot be recognized !!!  H3 pattern: %s' % H3_pattern)) # C[A-Z]{1,33}(W)(G|A|C)[A-Z](Q|S|G|R)
		raise ValueError()
		#sys.exit(1)

	res = dict(L1=L1, L2=L2, L3=L3, H1=H1, H2=H2, H3=H3,  FR_L1=FR_L1, FR_L2=FR_L2, FR_L3=FR_L3, FR_L4=FR_L4,  FR_H1=FR_H1, FR_H2=FR_H2, FR_H3=FR_H3, FR_H4=FR_H4)
	#if Options.verbose: print 'L1: %(L1)s\nL2: %(L2)s\nL3: %(L3)s\nH1: %(H1)s\nH2: %(H2)s\nH3: %(H3)s' % res

	return res


def int_(s):
	""" Removes all capital letters from its argument, returns int value of remainder
	"""
	return int( re.sub('[A-Z]', '', s) )  #v = int( re.sub('[A-Z]', '', new_number_FR_L1) )  # $new_number_FR_L1[$i] =~ s/[A-Z]//	 #new_number_FR_L1[i] = string.translate(new_number_FR_L1[i], None, string.ascii_letters)

def Extract_FR_CDR_Sequences(L1='', L2='', L3='', H1='', H2='', H3='', FR_L1='', FR_L2='', FR_L3='', FR_L4='', FR_H1='', FR_H2='', FR_H3='', FR_H4=''):
	""" Find Cothia numbering for loops and conserved regions
	"""

	# LIGHT CHAIN
	# FR_L1	How can we handle missing residue in C/N-terminals?
	if re.search( r'[A-Z][QE][A-Z]{9}[A-Z][A-Z]{4}[LVIMF][A-Z]C', FR_L1): # Change G to [A-Z] (3G04)
		if   len(FR_L1) == 19: new_number_FR_L1="5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
		elif len(FR_L1) == 20: new_number_FR_L1="4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
		elif len(FR_L1) == 21: new_number_FR_L1="3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
		elif len(FR_L1) == 22: new_number_FR_L1="2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
		elif len(FR_L1) == 23: new_number_FR_L1="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
		elif len(FR_L1) == 24:
			#new_number_FR_L1="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
			FR_L1 = FR_L1[1:]  # Remove 0th residue 10/24/2012
			new_number_FR_L1="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
		else: print(("ERROR: FR_L1 matches [A-Z][QE][A-Z]{9}[A-Z][A-Z]{4}[LVIMF][A-Z]C but length",len(FR_L1),"is not between 19 and 24"))

	elif re.search( r'[A-Z][QE][A-Z]{8}[A-Z][A-Z]{4}[LVIMF][A-Z]C', FR_L1):  # Change G to [A-Z] (3G04)
		if   len(FR_L1) == 19: new_number_FR_L1="4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
		elif len(FR_L1) == 20: new_number_FR_L1="3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
		elif len(FR_L1) == 21: new_number_FR_L1="2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
		elif len(FR_L1) == 22: new_number_FR_L1="1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
		elif len(FR_L1) == 23:
			FR_L1 = FR_L1[1:]  # Remove 0th residue 10/24/2012
			new_number_FR_L1="1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
		elif len(FR_L1) == 24:
			FR_L1 = FR_L1[2:]  # Remove -1st and 0th residue 10/24/2012
			new_number_FR_L1="1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
		else: print(("ERROR: FR_L1 matches [A-Z][QE][A-Z]{8}[A-Z][A-Z]{4}[LVIMF][A-Z]C but length",len(FR_L1),"is not between 19 and 24"))
	else:
		print('ERROR: Current code could not assign Chothia numbering of FR_L1 in the query sequence!!! Exiting...')

	# L1
	if   len(L1) ==  8: new_number_L1="24,25,26,27,28,29,30,34"
	elif len(L1) ==  9: new_number_L1="24,25,26,27,28,29,30,33,34"
	elif len(L1) == 10: new_number_L1="24,25,26,27,28,29,30,32,33,34"
	elif len(L1) == 11: new_number_L1="24,25,26,27,28,29,30,31,32,33,34"
	elif len(L1) == 12: new_number_L1="24,25,26,27,28,29,30,30A,31,32,33,34"
	elif len(L1) == 13: new_number_L1="24,25,26,27,28,29,30,30A,30B,31,32,33,34"
	elif len(L1) == 14: new_number_L1="24,25,26,27,28,29,30,30A,30B,30C,31,32,33,34"
	elif len(L1) == 15: new_number_L1="24,25,26,27,28,29,30,30A,30B,30C,30D,31,32,33,34"
	elif len(L1) == 16: new_number_L1="24,25,26,27,28,29,30,30A,30B,30C,30D,30E,31,32,33,34"
	elif len(L1) == 17: new_number_L1="24,25,26,27,28,29,30,30A,30B,30C,30D,30E,30F,31,32,33,34"
	else: print(("ERROR: L1 length",len(L1),"is not between 8 and 17"))

	# FR_L2
	if len(FR_L2) == 15: new_number_FR_L2="35,36,37,38,39,40,41,42,43,44,45,46,47,48,49"
	else: print(("ERROR: FR_L2 length",len(FR_L2),"is not 15"))

	# L2
	if   len(L2) ==  7: new_number_L2="50,51,52,53,54,55,56"
	elif len(L2) ==  8: new_number_L2="50,51,52,53,54,54A,55,56"
	elif len(L2) ==  9: new_number_L2="50,51,52,53,54,54A,54B,55,56"
	elif len(L2) == 10: new_number_L2="50,51,52,53,54,54A,54B,54C,55,56"
	elif len(L2) == 11: new_number_L2="50,51,52,53,54,54A,54B,54C,54D,55,56"
	else: print(("ERROR: L2 length",len(L2),"is not between 7 and 11"))

	# FR_L3
	if   len(FR_L3) == 32: new_number_FR_L3="57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88"
	elif len(FR_L3) == 33: new_number_FR_L3="57,58,59,60,61,62,63,64,65,66,66A,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88"
	elif len(FR_L3) == 34: new_number_FR_L3="57,58,59,60,61,62,63,64,65,66,66A,66B,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88"
	else: print(("ERROR: FR_L3 length",len(FR_L3),"is not between 32 and 34"))

	# L3
	if   len(L3) ==  5: new_number_L3="89,90,91,92,97"
	elif len(L3) ==  6: new_number_L3="89,90,91,92,93,97"
	elif len(L3) ==  7: new_number_L3="89,90,91,92,93,94,97"
	elif len(L3) ==  8: new_number_L3="89,90,91,92,93,94,95,97"
	elif len(L3) ==  9: new_number_L3="89,90,91,92,93,94,95,96,97"
	elif len(L3) == 10: new_number_L3="89,90,91,92,93,94,95,95A,96,97"
	elif len(L3) == 11: new_number_L3="89,90,91,92,93,94,95,95A,95B,96,97"
	elif len(L3) == 12: new_number_L3="89,90,91,92,93,94,95,95A,95B,95C,96,97"
	elif len(L3) == 13: new_number_L3="89,90,91,92,93,94,95,95A,95B,95C,95D,96,97"
	elif len(L3) == 14: new_number_L3="89,90,91,92,93,94,95,95A,95B,95C,95D,95E,96,97"
	elif len(L3) == 15: new_number_L3="89,90,91,92,93,94,95,95A,95B,95C,95D,95E,95F,96,97"
	else: print(("ERROR: L3 length",len(L3),"is not between 5 and 15"))

	# FR_L4
	new_number_FR_L4="98,99,100,101,102,103,104,105,106,107,108,109"

	# HEAVY CHAIN
	# FR_H1	How can we handle missing residue in C/N-terminals?
	if   len(FR_H1) == 16: new_number_FR_H1="10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
	elif len(FR_H1) == 17: new_number_FR_H1="9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
	elif len(FR_H1) == 18: new_number_FR_H1="8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
	elif len(FR_H1) == 19: new_number_FR_H1="7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
	elif len(FR_H1) == 20: new_number_FR_H1="6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
	elif len(FR_H1) == 21: new_number_FR_H1="5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
	elif len(FR_H1) == 22: new_number_FR_H1="4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
	elif len(FR_H1) == 23: new_number_FR_H1="3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
	elif len(FR_H1) == 24: new_number_FR_H1="2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
	elif len(FR_H1) == 25: new_number_FR_H1="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
	elif len(FR_H1) == 26:
		#new_number_FR_H1="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
		new_number_FR_H1="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
		FR_H1 = FR_H1[1:]  # Remove 0th residue 10/24/2012
	else: print(("ERROR: FR_H1 length",len(FR_H1),"is not between 16 and 26"))

	# H1
	if   len(H1) ==  6: new_number_H1="26,27,32,33,34,35"
	elif len(H1) ==  7: new_number_H1="26,27,28,32,33,34,35"
	elif len(H1) ==  8: new_number_H1="26,27,28,29,32,33,34,35"
	elif len(H1) ==  9: new_number_H1="26,27,28,29,30,32,33,34,35"
	elif len(H1) == 10: new_number_H1="26,27,28,29,30,31,32,33,34,35"
	elif len(H1) == 11: new_number_H1="26,27,28,29,30,31,31A,32,33,34,35"
	elif len(H1) == 12: new_number_H1="26,27,28,29,30,31,31A,31B,32,33,34,35"
	elif len(H1) == 13: new_number_H1="26,27,28,29,30,31,31A,31B,31C,32,33,34,35"
	else: print(("ERROR: H1 length",len(H1),"is not between 6 and 13"))

	# FR_H2
	if len(FR_H2) == 14: new_number_FR_H2="36,37,38,39,40,41,42,43,44,45,46,47,48,49"
	else: print(("ERROR: FR_H2 length",len(FR_H2),"is not 14"))

	# H2
	if   len(H2) == 12: new_number_H2="50,51,52,57,58,59,60,61,62,63,64,65"
	elif len(H2) == 13: new_number_H2="50,51,52,56,57,58,59,60,61,62,63,64,65"
	elif len(H2) == 14: new_number_H2="50,51,52,55,56,57,58,59,60,61,62,63,64,65"
	elif len(H2) == 15: new_number_H2="50,51,52,54,55,56,57,58,59,60,61,62,63,64,65"
	elif len(H2) == 16: new_number_H2="50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65"
	elif len(H2) == 17: new_number_H2="50,51,52,52A,53,54,55,56,57,58,59,60,61,62,63,64,65"
	elif len(H2) == 18: new_number_H2="50,51,52,52A,52B,53,54,55,56,57,58,59,60,61,62,63,64,65"
	elif len(H2) == 19: new_number_H2="50,51,52,52A,52B,52C,53,54,55,56,57,58,59,60,61,62,63,64,65"
	elif len(H2) == 20: new_number_H2="50,51,52,52A,52B,52C,52D,53,54,55,56,57,58,59,60,61,62,63,64,65"
	elif len(H2) == 21: new_number_H2="50,51,52,52A,52B,52C,52D,52E,53,54,55,56,57,58,59,60,61,62,63,64,65"
	elif len(H2) == 22: new_number_H2="50,51,52,52A,52B,52C,52D,52E,52F,53,54,55,56,57,58,59,60,61,62,63,64,65"
	else: print(("ERROR: H2 length",len(H2),"is not between 12 and 22"))

	# FR_H3
	if   len(FR_H3) == 30: new_number_FR_H3="66,67,68,69,70,71,72,73,76,77,78,79,80,81,82,82A,82B,82C,83,84,85,86,87,88,89,90,91,92,93,94"
	elif len(FR_H3) == 31: new_number_FR_H3="66,67,68,69,70,71,72,73,74,76,77,78,79,80,81,82,82A,82B,82C,83,84,85,86,87,88,89,90,91,92,93,94"
	elif len(FR_H3) == 32: new_number_FR_H3="66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,82A,82B,82C,83,84,85,86,87,88,89,90,91,92,93,94"
	else: print(("ERROR: FR_H3 length",len(FR_H3),"is not between 30 and 32"))

	# H3
	if   len(H3) ==  3: new_number_H3="95,101,102"
	elif len(H3) ==  4: new_number_H3="95,96,101,102"
	elif len(H3) ==  5: new_number_H3="95,96,97,101,102"
	elif len(H3) ==  6: new_number_H3="95,96,97,98,101,102"
	elif len(H3) ==  7: new_number_H3="95,96,97,98,99,101,102"
	elif len(H3) ==  8: new_number_H3="95,96,97,98,99,100,101,102"
	elif len(H3) ==  9: new_number_H3="95,96,97,98,99,100,100A,101,102"
	elif len(H3) == 10: new_number_H3="95,96,97,98,99,100,100A,100B,101,102"
	elif len(H3) == 11: new_number_H3="95,96,97,98,99,100,100A,100B,100C,101,102"
	elif len(H3) == 12: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,101,102"
	elif len(H3) == 13: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,101,102"
	elif len(H3) == 14: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,101,102"
	elif len(H3) == 15: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,101,102"
	elif len(H3) == 16: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,101,102"
	elif len(H3) == 17: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,101,102"
	elif len(H3) == 18: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,101,102"
	elif len(H3) == 19: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,101,102"
	elif len(H3) == 20: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,101,102"
	elif len(H3) == 21: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,101,102"
	elif len(H3) == 22: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,101,102"
	elif len(H3) == 23: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,101,102"
	elif len(H3) == 24: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,101,102"
	elif len(H3) == 25: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,101,102"
	elif len(H3) == 26: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,101,102"
	elif len(H3) == 27: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,101,102"
	elif len(H3) == 28: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,101,102"
	elif len(H3) == 29: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,100U,101,102"
	elif len(H3) == 30: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,100U,100V,101,102"
	elif len(H3) == 31: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,100U,100V,100W,101,102"
	elif len(H3) == 32: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,100U,100V,100W,100X,101,102"
	elif len(H3) == 33: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,100U,100V,100W,100X,100Y,101,102"
	elif len(H3) == 34: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,100U,100V,100W,100X,100Y,100Z,101,102"
	else: print(("ERROR: H3 length",len(FR_H3),"is not between 3 and 34"))

	# FR_H4
	new_number_FR_H4="103,104,105,106,107,108,109,110,111,112,113,114"

	try:
		(new_number_L1 and new_number_L2 and new_number_L3 and new_number_H1 and new_number_H2 and new_number_H3
		 and new_number_FR_H1 and new_number_FR_H2 and new_number_FR_H3 and new_number_FR_H4
		 and new_number_FR_L1 and new_number_FR_L2 and new_number_FR_L3 and new_number_FR_L4 )
	except:
		print("Numbering failed.  Exiting.")
		sys.exit(1)

	# Converting all new_number_* vars in to a lists
	new_number_L1=new_number_L1.split(',');  new_number_L2=new_number_L2.split(',');  new_number_L3=new_number_L3.split(',');
	new_number_H1=new_number_H1.split(',');  new_number_H2=new_number_H2.split(',');  new_number_H3=new_number_H3.split(',');
	new_number_FR_L1=new_number_FR_L1.split(',');  new_number_FR_L2=new_number_FR_L2.split(',');  new_number_FR_L3=new_number_FR_L3.split(',');  new_number_FR_L4=new_number_FR_L4.split(',');
	new_number_FR_H1=new_number_FR_H1.split(',');  new_number_FR_H2=new_number_FR_H2.split(',');  new_number_FR_H3=new_number_FR_H3.split(',');  new_number_FR_H4=new_number_FR_H4.split(',');

	print_seq_FR_L1, print_seq_FR_L2, print_seq_FR_L3, print_seq_FR_L4, print_seq_FR_L4_extra = '', '', '', '', ''
	print_seq_FR_H1, print_seq_FR_H2, print_seq_FR_H3, print_seq_FR_H4, print_seq_FR_H4_extra = '', '', '', '', ''
	numbering_L, numbering_H = {}, {}

	# OUTPUT FOR LIGHT CHAIN. This should be save in the file 'numbering_L.txt'.
	for i, s in enumerate(FR_L1): #FR_L1
		numbering_L[ new_number_FR_L1[i] ] = s  #+='%s %s\n' % (s, new_number_FR_L1[i])
		v = int_(new_number_FR_L1[i])
		if v >= 10 and v <= 23: print_seq_FR_L1 += s

	for i, s in enumerate(L1): numbering_L[ new_number_L1[i] ] = s  # +='%s %s\n' % (s, new_number_L1[i])  # L1

	for i, s in enumerate(FR_L2):  #FR_L2
		numbering_L [ new_number_FR_L2[i] ] = s  #+='%s %s\n' % (s, new_number_FR_L2[i])
		v = int_(new_number_FR_L2[i])
		if (v >= 35 and v <= 38) or (v >= 45 and v <= 49): print_seq_FR_L2 += s

	for i, s in enumerate(L2): numbering_L[ new_number_L2[i] ] = s  #+='%s %s\n' % (s, new_number_L2[i])  # L2

	for i, s in enumerate(FR_L3):  #FR_L3
		numbering_L[ new_number_FR_L3[i] ] = s  # +='%s %s\n' % (s, new_number_FR_L3[i])
		v = int_(new_number_FR_L3[i])
		if (v >= 57 and v <= 66) or (v >= 71 and v <= 88): print_seq_FR_L3 += s

	for i, s in enumerate(L3): numbering_L[ new_number_L3[i] ] = s  #+='%s %s\n' % (s, new_number_L3[i])  # L3

	for i, s in enumerate(FR_L4):  #FR_L4
		numbering_L[ new_number_FR_L4[i] ] = s  #+='%s %s\n' % (s, new_number_FR_L4[i])
		v = int_(new_number_FR_L4[i])
		if v >= 98 and v <= 107: print_seq_FR_L4 += s
		if v >= 98 and v <= 101: print_seq_FR_L4_extra += s


	# OUTPUT FOR HEAVY CHAIN. This should be save in the file 'numbering_H.txt'.
	for i, s in enumerate(FR_H1):  #FR_H1
		numbering_H[ new_number_FR_H1[i] ] = s  #+='%s %s\n' % (s, new_number_FR_H1[i])
		v = int_(new_number_FR_H1[i])
		if v >= 10 and v <= 25: print_seq_FR_H1 += s

	for i, s in enumerate(H1): numbering_H[ new_number_H1[i] ] = s  #+='%s %s\n' % (s, new_number_H1[i])  # H1

	for i, s in enumerate(FR_H2):  #FR_H2
		numbering_H[ new_number_FR_H2[i] ] = s  #+='%s %s\n' % (s, new_number_FR_H2[i])
		v = int_(new_number_FR_H2[i])
		if (v >= 36 and v <= 39) or (v >= 46 and v <= 49): print_seq_FR_H2 += s

	for i, s in enumerate(H2): numbering_H[ new_number_H2[i] ] = s  #+='%s %s\n' % (s, new_number_H2[i])  # H2

	for i, s in enumerate(FR_H3):  #FR_H3
		numbering_H[ new_number_FR_H3[i] ] = s  #+='%s %s\n' % (s, new_number_FR_H3[i])
		v = int_(new_number_FR_H3[i])
		if v >= 66 and v <= 94: print_seq_FR_H3 += s

	for i, s in enumerate(H3): numbering_H[ new_number_H3[i] ] = s  #+='%s %s\n' % (s, new_number_H3[i])  # H3

	for i, s in enumerate(FR_H4):  #FR_H4
		numbering_H[ new_number_FR_H4[i] ] = s  #+='%s %s\n' % (s, new_number_FR_H4[i])
		v = int_(new_number_FR_H4[i])
		if v >= 103 and v <= 112: print_seq_FR_H4 += s
		if v >= 103 and v <= 106: print_seq_FR_H4_extra += s

	#if Options.verbose: print 'numbering_L:\n', numbering_L
	#if Options.verbose: print 'numbering_H:\n', numbering_H

	FRL = print_seq_FR_L1 + print_seq_FR_L2 + print_seq_FR_L3 + print_seq_FR_L4
	FRH = print_seq_FR_H1 + print_seq_FR_H2 + print_seq_FR_H3 + print_seq_FR_H4

	if len(FRL) != 58 and len(FRL) != 60 and len(FRL) != 61 and len(FRL) != 63 :
		print("ERROR: Current DB does not cover the length of FRL of your query.")
		print(("ERROR: FRL length of your query:", len(FRL)))
		print("ERROR: DB: 58 or 60")
		sys.exit(1)

	if len(FRH) != 63 and len(FRH) != 65 and len(FRH) != 66 and len(FRH) != 68:
		print("ERORR: Current DB does not cover the length of FRL of your query.")
		print(("ERORR: FRH length of your query:", len(FRH)))
		print("ERORR: DB: 63 or 65")
		sys.exit(1)

	return dict(FRL = print_seq_FR_L1 + print_seq_FR_L2 + print_seq_FR_L3 + print_seq_FR_L4,
				FRH = print_seq_FR_H1 + print_seq_FR_H2 + print_seq_FR_H3 + print_seq_FR_H4,
				light = FR_L1 + L1 + FR_L2 + L2 + FR_L3 + L3 + print_seq_FR_L4,
				heavy = FR_H1 + H1 + FR_H2 + H2 + FR_H3 + H3 + print_seq_FR_H4,
				numbering_L=numbering_L, numbering_H=numbering_H)


for root, dirs, files in os.walk(directory):
	for file in files:
		print(os.path.splitext(os.path.splitext(file)[0])[0])
		if os.path.splitext(os.path.splitext(file)[0])[0] in pdb_list:
			code = os.path.splitext(file)[0][:4]	
			print(code)
			PDBfile = gzip.open(os.path.join(directory, file), 'rb')
			parser = PDBParser()
			ppb = PPBuilder()
			structure = parser.get_structure(code, PDBfile)
			model = structure[0]
			heavy = ''
			light = ''
			for pp in ppb.build_peptides(model['H']):
				heavy += str(pp.get_sequence())
				print(pp.get_sequence())
			for pp in ppb.build_peptides(model['L']):
				light += str(pp.get_sequence())
				print(pp.get_sequence())
			try:
				CDRs = IdentifyCDRs(light, heavy)
				full_seq[code] = heavy + '/' + light
				heavy_seq[code] = heavy
				light_seq[code] = light
				L1_seq[code] = CDRs['L1']
				L2_seq[code] = CDRs['L2']
				L3_seq[code] = CDRs['L3']
				H1_seq[code] = CDRs['H1']
				H2_seq[code] = CDRs['H2']
				H3_seq[code] = CDRs['H3']
				FR_L1_seq[code] = CDRs['FR_L1']
				FR_L2_seq[code] = CDRs['FR_L2']
				FR_L3_seq[code] = CDRs['FR_L3']
				FR_L4_seq[code] = CDRs['FR_L4']
				FR_H1_seq[code] = CDRs['FR_H1']
				FR_H2_seq[code] = CDRs['FR_H2']
				FR_H3_seq[code] = CDRs['FR_H3']
				FR_H4_seq[code] = CDRs['FR_H4']
				used += 1
			except ValueError as e:
				errors += 1


write_fasta_file('full', full_seq)
write_fasta_file('heavy', heavy_seq)
write_fasta_file('light', light_seq)
write_fasta_file('L1', L1_seq)
write_fasta_file('L2', L2_seq)
write_fasta_file('L3', L3_seq)
write_fasta_file('H1', H1_seq)
write_fasta_file('H2', H2_seq)
write_fasta_file('H3', H3_seq)
write_fasta_file('FR_L1', FR_L1_seq)
write_fasta_file('FR_L2', FR_L2_seq)
write_fasta_file('FR_L3', FR_L3_seq)
write_fasta_file('FR_L4', FR_L4_seq)
write_fasta_file('FR_H1', FR_H1_seq)
write_fasta_file('FR_H2', FR_H2_seq)
write_fasta_file('FR_H3', FR_H3_seq)
write_fasta_file('FR_H4', FR_H4_seq)
print(used)
