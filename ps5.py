#! /usr/bin/env python

from collections import Counter
import pybedtools

#problem 1.1 On what chromosome is the region with the largest start position 
#(2nd column) in lamina.bed?

filename = '/vol3/home/johnsonr/data-sets/bed/lamina.bed'

lamina = pybedtools.BedTool(filename)

bigstart = 0
bigchrom = ''
for record in lamina:
	start = int(record.start)
	chrom = record.chrom
	if start < bigstart: continue
	elif start > bigstart:
		bigstart = start
		bigchrom = chrom
print 'answer-1:', bigchrom 

#problem 1.2 What is the region with the largest end position on chrY in lamina.bed?

bigend = 0
for record in lamina:
	start = record.start
	chrom = record.chrom
	end = record.end
	region_length = record.end - record.start
	if chrom != 'chrY': continue
	elif chrom == 'chrY':
		if end < bigend: continue
		elif end > bigend:
			bigend = end
print 'answer-2:', chrom,':', start,'-',end

#problem 2.1 Which of the first 10 sequence records has the largest number of 'C' residues in the sequence? Report its record name 

filename = '/vol3/home/johnsonr/data-sets/fastq/SP1.fq'

def parse_fastq(filename):
	line_num = 0
	num_records = 0
	for line in open(filename):
		line_type = line_num % 4
		if line_type == 0:
			name = line.strip()
		elif line_type == 1:
			seq = line.strip()
		elif line_type == 3:
			quals = line.strip()
			yield name, seq, quals
		line_num += 1		 
def sum_quals(quals):
	sum = 0 
	for char in quals:
		sum += ord(char)
	return sum

def reverse_complement(seq):
	comps = []
	for char in reversed(seq):
		if char == 'A':
			comps.append('T')
		elif char == 'T':
			comps.append('A')
		elif char == 'G':
			comps.append('C')
		elif char == 'C':
			comps.append('G')
	return comps

maxC = 0
namemaxC = ''
num_records = 0
sum = 0
maxsum = 0
revcomp10 = []
for name, seq, quals in parse_fastq(filename):
	countC = Counter(seq)
	num_records += 1
	revcomp = reverse_complement(seq)
	if num_records <= 10:
		if countC['C'] > maxC:
			maxC = countC['C']
			namemaxC = name
		revcomp10.append(revcomp)
	qual_sum = sum_quals(quals)
	if qual_sum > maxsum:
		maxsum = qual_sum

print 'answer-3:', namemaxC	
print 'answer-4:', maxsum
print 'answer-5:', revcomp10

