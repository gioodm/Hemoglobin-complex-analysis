#!/usr/bin/python3

__author__ = "Giorgia Del Missier"

import sys


MASA = {'A':115., 'L':170., 'R':225., 'K': 200., 'N':160., 'M':185., 'D':150., \
		'F':210., 'C':135., 'P':145., 'Q': 180., 'S':115., 'E':190., 'T':140., \
		'G':75., 'W': 255., 'H':195., 'Y': 230., 'I':175., 'V':155.}

OneToThreeAAs = {'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP', 'C':'CYS', 'Q':'GLN', 'E':'GLU', \
				 'G':'GLY', 'H':'HIS', 'I':'ILE', 'L':'LEU', 'K':'LYS', 'M':'MET', 'F':'PHE', \
				 'P':'PRO', 'S':'SER', 'T':'THR', 'W':'TRP', 'Y':'TYR', 'V':'VAL'}


def parseDSSP(dsspfile):
	fdssp = open(dsspfile)
	dssp, count = {}, 0

	for line in fdssp:
		if line.find('  #  RESIDUE') == 0:
			count = 1
			continue

		if count == 1:
			num = line[5:10].strip()
			ch = line[11]
			res = line[13]
			ss = line[16]
			asa = float(line[35:38])
			phi, psi = float(line[103:109]), float(line[109:115])
			if ss == ' ': ss = 'C'

			try:
				dssp[(int(num),ch)]=[res, ss, asa, phi, psi]
			except:
				pass
		
	return dssp

	
def printRASA(dssp):
	keys = dssp.keys()
	list(keys).sort()

	for key in keys:
		print(key, dssp[key][2] / MASA[dssp[key][0]])


def RASAChain(dssp, ch):
	RASA_list = []
	for key, value in dssp.items():
		if key[1] == ch:
			rasa = min(1, dssp[key][2]/ MASA[dssp[key][0]])
			RASA_list.append((key[0], rasa, value[0]))
			RASA_list.sort()

	return RASA_list


def compareRASA(l1, l2, threshold = 0.10):
	significative_list = []

	for i in range(len(l1)):
		diff = round(l2[i][1] - l1[i][1], 4) 
		if diff > threshold:
			res = OneToThreeAAs[l1[i][2]]
			significative_list.append((l1[i][0], diff, res))

	return significative_list
	

def extract(l, res):
	for e in l:
		if e[0] == res:
			return e[1]


def print_significative(significative, tetralist, monolist, chain):
	print("CHAIN\tRESn\tRSA(M)\tRSA(C)\tRSA(M)-RSA(C)")

	for res, diff, restype in significative:
		cRASA = round(extract(tetralist, res), 4)
		mRASA = round(extract(monolist, res), 4)
		print(chain + "\t" + restype + str(res) + "\t" + str(mRASA) + "\t" + str(cRASA) + "\t" + str(diff))

				

if __name__ == '__main__':

	if len(sys.argv) == 2:
		dsspfile = sys.argv[1]
		dssp = parseDSSP(dsspfile)
		print(printRASA(dssp))

	elif len(sys.argv) == 4:
		dsspfile1 = sys.argv[1]
		dsspfile2 = sys.argv[2]
		dssp1 = parseDSSP(dsspfile1)
		dssp2 = parseDSSP(dsspfile2)
		chain = sys.argv[3]
		tetralist = RASAChain(dssp1, chain)
		monolist = RASAChain(dssp2, chain)
		list_significatives = compareRASA(tetralist, monolist)

		print_significative(list_significatives, tetralist, monolist, chain)
	

	else:
		print('parsedssp2.py dsspfile -o dsspchain chain')

