#!/bin/env python

import glob, os, sys, time, fnmatch, getpass, subprocess, urllib, urllib.request, urllib.error, urllib.parse, webbrowser, re, shutil, gzip, operator
from shutil import copyfile

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

def write_fasta_file(file_name, data, prefix):
	""" writes fasta file

	file_name -- name of file, will also be header / sequence name
	data -- the single sequence to be written
	prefix -- path prepended to the filename, explicitly have a directory separator if needed
	"""
	with open(prefix+file_name+'.fasta', 'w') as f: f.write('>%s\n' % file_name);  f.write(data); f.write('\n')

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
		sys.exit(1)

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

def Blast(folder, folderall, minscore, maxscore, threshold):
	os.system('mkdir ' + folder)
	if os.path.exists('./' + folder +'/results.html'): #Skip Blast if results already exist
		pass
	else: #Submit fragment to Blast and download results
		Blasturl = urllib.request.Request('https://blast.ncbi.nlm.nih.gov/Blast.cgi?')
		BlastFile = open('./' + folderall +'/' + folder + '.fasta', 'rb')
		BlastuploadText = BlastFile.read()
		BlastFile.close()
		Blastdata = urllib.parse.urlencode({'QUERY': BlastuploadText,'PROGRAM': 'blastp','DATABASE': 'pdb','MAX_NUM_SEQ': '5000','HITLIST_SIZE': '100','DESCRIPTIONS': '5000','ALIGNMENTS': '100', 'SHORT_QUERY_ADJUST': 'on', 'CMD': 'Put'}).encode('utf-8')
		BlastIDresults = urllib.request.urlopen(Blasturl, Blastdata)
		with open('./' + folder +'/IDresults.html', 'wb') as f:
			f.write(BlastIDresults.read())
			f.close
		BlastResultsFile = open('./' + folder +'/IDresults.html', 'r')
		HTMLresults = BlastResultsFile.read()
		start = 0
		if HTMLresults.find('<input name="RID" value="', start) > 0:
			idx = HTMLresults.find('<input name="RID" value="', start)
			closing_anchor = HTMLresults.find('"', idx+25)
			scoreID = HTMLresults[idx+25:closing_anchor]
		time.sleep(30)
		Blastdata = urllib.parse.urlencode({'CMD': 'Get','RID': scoreID,'DESCRIPTIONS': '5000',}).encode('utf-8')
		Blastresults = urllib.request.urlopen(Blasturl, Blastdata)
		with open('./' + folder +'/results.html', 'wb') as f:
			f.write(Blastresults.read())
			f.close()

	with open('./' + folder +'/results.html', 'rb') as f:
		text = f.read().decode('utf-8')

	start = 0
	while text.find('Status=WAITING', start) > 0: #Check if Blast is finished
		print('Not Ready')
		time.sleep(10)
		f.close()
		Blastresults = urllib.request.urlopen(Blasturl, Blastdata)
		with open('./' + folder +'/results.html', 'wb') as f:
			f.write(Blastresults.read())
			f.close
		with open('./' + folder +'/results.html', 'rb') as f:
			text = f.read().decode('utf-8')
		start = 0

	start = 0
	q = []
	added = 0
	while text.find('<td class="c1 l lim">', start) > 0: #Look through results file for list of PDBs and score
		idx = text.find('<td class="c1 l lim">', start)
		closing_anchor = text.find('</a>', idx)
		title = text[closing_anchor-6:closing_anchor]
		title = title.replace('_', ' ')
		altidx = text.find('c7', idx - 65)
		altclose = text.find('%', altidx)
		score = int(float(text[altidx+4:altclose]))
		dr = ''
		if title[0:4].lower() in antibody_list and title[0:4].lower() != folderall[0:4].lower() and minimum < score < maximum: # Filter by database
			print(title[0:4].lower() )
			if folder == 'heavy' and len(score_list[title[0:4].lower()]) == 0:
				score_list[title[0:4].lower()].append(score)
			if folder == 'light' and len(score_list[title[0:4].lower()]) == 1:
				score_list[title[0:4].lower()].append(score)
		start = idx + 40

folderall = sys.argv[1].replace(' ', '_').upper()
folderh = sys.argv[1].replace(' ', '_')[0:6].upper()
folderl = sys.argv[1].replace(' ', '_')[0:5].upper() + sys.argv[1][6].upper()
for the_file in os.listdir('./' + folderall):
	file_path = os.path.join('./' + folderall, the_file)
	try:
		if os.path.isfile(file_path):
			os.unlink(file_path)
	except Exception as e:
		print(e)
os.system('( cd ' + folderall + ' && ~schmits/alt/bin/python3.6 ../clean_pdb.py ' + folderh.replace('_', ' ') + ' )')
os.system('( cd ' + folderall + ' && ~schmits/alt/bin/python3.6 ../clean_pdb.py ' + folderl.replace('_', ' ') + ' )')
score_list = {}

pdb = sys.argv[1][0:4].upper()

minimum = 0
maximum = 98





os.system('mkdir ' + folderall + '/Single_Control')

antibody_list = ['4ydl', '4ydi', '4ydk', '4ydj', '1vge', '5ig7', '4kq3', '5y9k', '3x3f', '1n0x', '1i9r', '5veb', '2fl5', '6bhz', '3k2u', '2r56', '4x0k', '5xhv', '3eo9', '5wuv', '4py7', '3eyf', '5cd3', '5wcd', '5wca', '3qrg', '2d7t', '5u0r', '4lsu', '4lss', '1jps', '4x7s', '4zs6', '5jrp', '4evn', '1bvk', '4odx', '5x8m', '5mvz', '4zso', '3sob', '5jr1', '3pgf', '3oaz', '3oau', '5igx', '1iqd', '1mim', '4edw', '4gxv', '5d7s', '3grw', '1cly', '3u30', '4fze', '5d72', '3b2u', '4hha', '6b0w', '6b0e', '6b0g', '6b0a', '6b0h', '5cex', '4rrp', '4f57', '1w72', '6b08', '4f58', '4nzr', '4nzu', '4uta', '6ani', '5u15', '3p30', '4xgz', '1za6', '4ut9', '1dql', '5ggq', '5ggs', '5ggt', '4ut7', '5ggv', '5vic', '5vig', '2wuc', '3lh2', '4ers', '3fzu', '2qqk', '4ogy', '5hi4', '4kmt', '4hie', '4nm4', '1h0d', '6c6z', '3zl4', '3n9g', '5d6c', '5tqa', '1jv5', '4tsa', '4tsc', '4tsb', '5t6l', '1fve', '1fvd', '5nhw', '4m1d', '4n0y', '5v7r', '5v7u', '4g5z', '5cck', '3skj', '3kdm', '4qci', '5tfs', '5tpp', '5d9q', '4dtg', '3w9e', '5c8j', '5dr5', '4r8w', '3qcu', '3gbn', '5f7e', '2yc1', '4ttd', '4lmq', '5ewi', '5cgy', '6ele', '5awn', '5i9q', '5whj', '1ikf', '5jw5', '4jdv', '5bo1', '5vob', '5jo4', '5ik3', '5w08', '5jof', '2qsc', '5gks', '4nki', '5f6h', '4eow', '3tnm', '5n4g', '5n4j', '2fb4', '4zd3', '5k9j', '2xra', '5te7', '5l6y', '2qqn', '3p0y', '5gzo', '3kym', '2b1h', '6b3d', '2uzi', '3g04', '4lst', '3hi6', '2eiz', '4k7p', '5bk2', '4g7v', '2qr0', '6fax', '4o5l', '3nfs', '4v1d', '3d85', '4yhl', '5drz', '5drw', '4jha', '1rz7', '4irz', '5f6i', '5tru', '5udc', '5trp', '1dee', '1rhh', '5xxy', '4mxv', '1wt5', '3qeh', '4uu9', '5ur8', '4zfg', '5tdn', '5eu7', '4g6a', '4g6f', '4g6m', '4xvj', '5mes', '3wd5', '4xvt', '3u0t', '5lsp', '2a9n', '4m6o', '4m6n', '3nh7', '4om0', '4om1', '5gmq', '4ypg', '3m8o', '5vsi', '5tpl', '5dum', '5dur', '5ush', '5usl', '5i1l', '5tgb', '5i1h', '5i1e', '5i1d', '5i1g', '5i1a', '5i1c', '4hwe', '1axs', '2cmr', '5i19', '5i15', '5i17', '5i16', '5wk2', '5uea', '2xtj', '4dke', '4dkf', '4jkp', '1gc1', '2h9g', '5bzd', '2yss', '4m5y', '4hpo', '5bzw', '2eks', '4xak', '1dn0', '4h8w', '2vh5', '5jz7', '5itb', '4olz', '4olx', '5ob5', '4olw', '4olv', '4olu', '5vl7', '5f9w', '3l5x', '5f9o', '4hf5', '5e8e', '5uby', '5ubz', '4yk4', '3h42', '5f96', '2zkh', '4u6v', '4hfw', '4hfu', '5tf1', '4nrx', '4nrz', '4jzn', '4jzo', '5esv', '5esa', '5tfw', '3fn0', '5e08', '2a9m', '1igm', '4hpy', '3idg', '4m62', '5ea0', '3kr3', '3bn9', '2g75', '5bk1', '5bk3', '5bk5', '4g80', '4ywg', '4ioi', '3uls', '4ptu', '4z5r', '4ky1', '4gsd', '3piq', '4hg4', '5fgc', '4oqt', '5xku', '5k59', '6bp2', '1opg', '5ucb', '4y5x', '3u4b', '4z0x', '5lbs', '5w6g', '5w6c', '4xc1', '3tcl', '1l7i', '4s1r', '4s1s', '3hmx', '4r26', '4ris', '5bjz', '3lmj', '2vxv', '2vxq', '5sx5', '3hc4', '3hc0', '5bvj', '5ibt', '5ibu', '4y5v', '4xml', '5bv7', '4xmp', '2xwt', '3gjf', '5vsh', '2ny1', '5tdo', '5y11', '5g64', '5j13', '4npy', '5d1x', '5d1z', '4np4', '4jam', '1g9m', '3idx', '4dn3', '4rx4', '4fqj', '4qxg', '5kmv', '4fnl', '4oaw', '4xnq', '4imk', '4xny', '4iml', '3lzf', '5b71', '5u7o', '4ps4', '8fab', '4yhy', '4yhz', '4nug', '4p59', '4r7d', '4yho', '5feh', '5vag', '5i1i', '3ma9', '4lkx', '5uoe', '3sdy', '4lkc', '5ty6', '4wuu', '4wuk', '4x4x', '4x4y', '3gkw', '3mac', '5i8o', '1ad0', '1ad9', '4hs8', '5tzu', '4hs6', '4jy5', '3u6r', '5kn5', '2xzc', '4j6r', '5u4r', '6ayz', '4zyk', '3qhf', '1t3f', '3sqo', '4dgy', '3hae', '2yk1', '3ghe', '3na9', '5vvf', '5uem', '5kw9', '3naa', '3nab', '3nac', '5uek', '4cni', '6be3', '3mlw', '2hfg', '6erx', '4ut6', '4llu', '4lly', '4qf1', '4qhl', '5c2b', '4qhu', '3bdy', '4al8', '5t33', '4hk0', '4hk3', '5sy8', '4ocr', '3dvg', '5fuz', '3dvn', '5it2', '6axl', '6axk', '4ygv', '4hcr', '1hez', '2aj3', '4xx1', '3ztn', '3g6a', '4i77', '3g6j', '3go1', '3hc3', '2jb5', '1uj3', '5l0q', '3n85', '5umn', '4fq2', '4rav', '2fjh', '3giz', '2fjf', '4mwf', '4fqq', '4fqi', '4fqk', '4fql', '4fqc', '5n2k', '5ezi', '5czx', '5czv', '1nl0', '3r1g', '3dif', '3u7y', '4nyl', '5alb', '3mme', '5waw', '5chn', '2f5a', '3mly', '4nwt', '5szf', '5t29', '1u6a', '4jpw', '3dgg', '3aaz', '5ifh', '5usi', '5anm', '4j4p', '5dhv', '5tzt', '2agj', '5fha', '4ye4', '5fhb', '4ot1', '3so3', '5tz2', '1gaf', '4qhm', '3ncj', '3inu', '4wv1', '1aqk', '5i8k', '5i8c', '5cin', '5cil', '4nnp', '3qot', '5vzy', '3se9', '4xxd', '2j6e', '4llw', '5t5b', '3mxw', '1dfb', '5iq9', '3nps', '5iq7', '3ujj', '3uji', '6az2', '5cus', '4d9q', '5u3j', '5u3o', '5u3n', '5u3m', '5u3l', '5u3p', '4uv7', '4uv4', '4d9l', '6azx', '6azz', '6azm', '5eii']
for antibody in antibody_list:
	score_list[antibody] = []

heavy_chain = read_fasta_file('./' + folderall +'/' + folderh + '.fasta')[0]
light_chain = read_fasta_file('./' + folderall +'/' + folderl + '.fasta')[0]

CDRs = IdentifyCDRs(light_chain, heavy_chain)
CDRs.update( Extract_FR_CDR_Sequences(**CDRs) )

write_fasta_file('heavy', CDRs['heavy'], './' + folderall +'/')
write_fasta_file('light', CDRs['light'], './' + folderall +'/')

Blast('heavy', folderall, minimum, maximum, 10)
Blast('light', folderall, minimum, maximum, 10)

for key, value in score_list.items():
	score_list[key] = sum(value)
print(score_list)
sorted_score_list = sorted(score_list, key=score_list.__getitem__)
sorted_score_list.reverse()
print(sorted_score_list)

#Get top 15 best combined file to run clustering
for key in sorted_score_list[:1]:
	best_template = key
	with gzip.open('./templates/' + key + '.pdb.gz', 'rb') as f_in:
		with open('./' + folderall + '/Single_Control/' + key + '.pdb', 'wb') as f_out:
			shutil.copyfileobj(f_in, f_out)
	os.system('~schmits/alt/bin/python3.6 pdb2fasta.py ' + './' + folderall + '/Single_Control/' + key + '.pdb > '+ './' + folderall + '/Single_Control/' + key + '.fasta')
	time.sleep(0.05)

q = []

name = 'Single_Control'

fname = os.path.join(sys.path[0] + '/' + folderall, name)
for root, dirs, files in os.walk(fname):
	for file in files:
		if file.endswith('.pdb'):
			os.system('~schmits/alt/bin/python3.6 pdb2fasta.py ' + fname +'/' + file + ' > '+ fname +'/' + file[0:4] + '.fasta')
heavy_light = CDRs['heavy'] + '/' + CDRs['light']
write_fasta_file(pdb, heavy_light.replace('\n', ''), fname + '/')
os.system('( cd ' + fname + ' && cat *.fasta > alignment.fasta )')
os.system('clustalo -i ' + fname + '/alignment.fasta -o ' + fname + '/alignment.aln')
os.remove(os.path.join('./' + folderall +'/', folderh + '.fasta'))
os.remove(os.path.join('./' + folderall +'/', folderl + '.fasta'))
os.remove(os.path.join('./' + folderall +'/', folderh + '.pdb'))
os.remove(os.path.join('./' + folderall +'/', folderl + '.pdb'))
referenceContents = ' '
pdbCode = []
pdbContents = []
isReference = 0
with open(fname + '/alignment.aln') as file:
	for line in file.readlines():
		if line.find('>') == -1 and isReference == 0:
			pdbContents[len(pdbContents) - 1] =  pdbContents[len(pdbContents) - 1] + line
		elif line.find('>') == -1 and isReference == 1:
			referenceContents = referenceContents + line
		else:
			if line.find(pdb.upper()) > -1 or line.find(pdb.lower()) > -1:
				isReference = 1
			else:
				isReference = 0
				pdbCode.append(line[1:].replace('\n', ''))
				pdbContents.append(' ')
pdbContents[0] = pdbContents[0].replace('\n','')
pdbContents[0] = pdbContents[0] + '-'

referenceContents = referenceContents.replace('\n','')
referenceContents = referenceContents.replace('/','-/')

newFileName = (pdb.lower() + "_" + pdbCode[0][:4].lower() + ".grishin")
f = open(os.path.join(fname, newFileName), 'w')
f.write("## {} {}.pdb\n".format(pdb.lower(), best_template))
f.write("#\n")
f.write("scores_from_program: 0\n")
f.write("0{}\n".format(referenceContents))
f.write("0{}\n".format(pdbContents[0]))
copyfile('./relax.options', fname + '/relax.options')
f.close()
os.system('( cd ' + fname + ' && /dors/meilerlab/apps/rosetta/rosetta-3.7/main/source/bin/relax.default.linuxgccrelease -s ' + best_template + '.pdb @relax.options -nstruct 1 )')
while not os.path.exists(os.path.join(fname, (best_template + "_0001.pdb.gz"))):
	#print(os.path.join(fname, (pdb + "_0001.pdb.gz")))
	time.sleep(1)
os.rename(os.path.join(fname, (best_template + "_0001.pdb.gz")), os.path.join(fname, (best_template + ".pdb.gz")))
os.system('( cd ' + fname + ' && /dors/meilerlab/apps/rosetta/rosetta-3.7/main/source/bin/partial_thread.default.linuxgccrelease -in:file:fasta {} -in:file:alignment {} -in:file:template_pdb {} )'.format(pdb.upper() + '.fasta',newFileName,best_template + ".pdb.gz"))
while not os.path.exists(os.path.join(fname, best_template + ".pdb.pdb")):
	time.sleep(1)
with open('./rosetta_cm.options', 'r') as file:
	Options = file.read()
	newOptions = Options.replace('InsertFastaHere', pdb + '.fasta')

with open('./rosetta_cm.xml', 'r') as file:
	XML = file.read()
	XML = XML.replace('3FRAGSFILE', pdb + '_3.frags')
	XML = XML.replace('9FRAGSFILE', pdb + '_9.frags')

with open('./run_array.slurm', 'r') as file:
	RunSlurm = file.read()

with open('./run_relax.slurm', 'r') as file:
	RunRelax = file.read()


with open(fname + '/rosetta_cm.options', 'w') as file:
	file.write(newOptions)
RunAllSlurm = RunSlurm.replace('PATHL', 'SingleControl')
RunAllRelax = RunRelax.replace('PATHL', 'SingleControlRelax')
RunAllSlurm = RunAllSlurm.replace('OPTIONSP', fname + '/rosetta_cm.options')
RunAllRelax = RunAllRelax.replace('OPTIONSP', fname + '/relax.options')
RunAllSlurm = RunAllSlurm.replace('NUMBEROFTRIALS', '2')
RunAllSlurm = RunAllSlurm.replace('THREADS', '10')
RunAllRelax = RunAllRelax.replace('THREADS', '10')

with open(fname + '/run_array.slurm', 'w') as file:
	file.write(RunAllSlurm)
with open(fname + '/run_relax.slurm', 'w') as file:
	file.write(RunAllRelax)
copyfile('./relax.options', fname + '/relax.options')
copyfile('./sc_parser.bash', fname + '/sc_parser.bash')
XMLstring = ' <Template pdb="' + best_template + '.pdb.pdb" cst_file="AUTO" weight=   1.000 />'
allXML = XMLstring + '\n'
with open(fname  + '/rosetta_cm.xml', 'w') as file:
	file.write(XML.replace('THREADING', allXML))

shutil.rmtree('heavy')
shutil.rmtree('light')


