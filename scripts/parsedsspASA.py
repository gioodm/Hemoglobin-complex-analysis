#!/usr/bin/python3

__author__ = "Giorgia Del Missier"

import sys


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

			dssp[(num,ch)]=[res, ss, asa, phi, psi]

	return dssp

		
def getTotalASA(dssp):
	total_ASA = 0.0
	for asa in dssp.values():
		total_ASA += asa[2]
	return total_ASA


def getTotalASAChain(dssp, ch):
	total_ASA = 0.0
	for key in dssp:
		if key[1] == ch:
			total_ASA += dssp[key][2]
	return total_ASA

				

if __name__ == '__main__':

	if len(sys.argv) == 2:
		dsspfile=sys.argv[1]
		dssp = parseDSSP(dsspfile)
		print('Total ASA: ' + str(getTotalASA(dssp)))

	elif len(sys.argv) == 4:
		tetramer = sys.argv[1]
		tetramer_dssp = parseDSSP(tetramer)
		trimer = sys.argv[2]
		trimer_dssp = parseDSSP(trimer)
		missing_chain = sys.argv[3]

		for chain in "ABCD":
			if chain != missing_chain:
				trimer_ASA = getTotalASAChain(trimer_dssp, chain)
				tetramer_ASA = getTotalASAChain(tetramer_dssp, chain)
				final_ASA = trimer_ASA - tetramer_ASA
				print("Surface of interaction between " + str(chain) + " and " + str(missing_chain) + " equal to: " + str(final_ASA))

	else:
		print('parsedssp.py dsspfile -o dssptrimer missing_chain')

