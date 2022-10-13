#!/usr/bin/env python

import sys
from Bio import Entrez
from Bio import AlignIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Align.Applications import ClustalwCommandline
from Bio import motifs
from Bio.Align import AlignInfo

if len(sys.argv) < 2:
	print("Missing query string")
	exit(-1)
stage = 100
if len(sys.argv) > 2:
	stage = int(sys.argv[2])
	
query = sys.argv[1]

#fetch sequences
seqfile = 'searchresults.fasta'
Entrez.email = "anupam06@pitt.edu"
records = Entrez.read(Entrez.esearch(db='protein',term=query))
fasta = open(seqfile,'w')
fasta.write(Entrez.efetch(db='protein',id=records['IdList'],rettype='fasta',retmode='text').read())
fasta.close()
#stage
if stage <= 70:
	sys.exit()
	
#align sequences
alignfile = 'align.aln'
cline = ClustalwCommandline('clustalw', infile=seqfile,outfile=alignfile)
cline()

#stage
if stage <= 80:
	sys.exit()
	
align = AlignIO.read(alignfile,'clustal')
summary = AlignInfo.SummaryInfo(align)

cons = summary.gap_consensus(threshold=0)
print(cons)

#stage
if stage <= 90:
	sys.exit()
	

#figure out what sequence match consensus in non gap regions
cnts = []
for seq in align:
	cnt = 0
	for i in range(len(cons)):
		if cons[i] != '-' and cons[i] == seq[i]:
			cnt += 1
	cnts.append((seq,cnt))

cnts = sorted(cnts, key=lambda nc: (-nc[1],nc[0].description))

for (s,n) in cnts[:10]:
	print(s.description,n)
	print(s.seq)


