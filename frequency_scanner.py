# Frequency_Scanner takes fastq files and calculates the frequency of every type of single nucletodie change 
# Copyright (C) 2019 Ryan C Shean 

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>. 


import glob
import subprocess
import argparse
import sys


def up_count(s1, s2, c1, c2):
	global a_count
	global c_count
	global t_count
	global g_count

	global a_to_c
	global a_to_t
	global a_to_g 

	global c_to_a 
	global c_to_t
	global c_to_g

	global t_to_a
	global t_to_c 
	global t_to_g

	global g_to_a
	global g_to_c
	global g_to_t

	if s1 == 'A':
		a_count += c1
		if s2 == 'C':
			a_to_c += c2
		elif s2 == 'T':
			a_to_t += c2
		elif s2 == 'G':
			a_to_g += c2
	elif s1 == 'C':
		c_count += c1
		if s2 == 'A':
			c_to_a += c2
		elif s2 == 'T':
			c_to_t += c2
		elif s2 == 'G':
			c_to_g += c2
	elif s1 == 'T':
		t_count += c1
		if s2 == 'A':
			t_to_a += c2
		elif s2 == 'C':
			t_to_c += c2
		elif s2 == 'G':
			t_to_g += c2
	elif s1 == 'G':
		g_count += c1
		if s2 == 'A':
			g_to_a += c2
		elif s2 == 'C':
			g_to_c += c2
		elif s2 == 'T':
			g_to_t += c2


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Frequency Scanner is used to calculate relative frequencies nucleotide transitions and transversions of low frequency random mutations. Takes a reference fasta and aligns all fasta files in the working directory to the reference and calculates mutational frequeinces based off this. If you are using paired end reads you need to have them fully merged at input.')
	parser.add_argument('reference_fasta', help='You need a reference fasta to align all your reads to. For best results this should be the majority consensus. Although you will get decent output so long as your input reads can align to this file')
	parser.add_argument('-qual', help='Specify a qualtiy filter, any read that has any base with a quality score lower than this will be completely removed. Default: 34')
	parser.add_argument('-af', help='Mutations with a frequency greater than this will not be counted')

	args = parser.parse_args()

	if args.qual:
		qual = args.qual
	else:
		qual = 34
	if args.af:
		af = args.af
	else:
		af = 0.05

	ref_fasta = args.reference_fasta
	subprocess.call('bwa index ' + ref_fasta, shell=True)

	for name in glob.glob('*.fastq'):
		base = name.split('.')[0]
		out = '-'.join(base.split('-')[0:3]) + '-' + str(qual)
		subprocess.call('fastq_quality_filter -q ' + str(qual) + ' -p 100 -i ' + name + ' -o ' + out + '.fastq', shell=True)
		
		subprocess.call('bwa mem ' + ref_fasta + ' ' + out + '.fastq > ' + base + '.bam',shell=True)
		subprocess.call('samtools sort ' + base + '.bam -o ' + base + '_sorted.bam', shell=True)
		subprocess.call('./bin/bam-readcount -f ' + ref_fasta + ' ' + base + '_sorted.bam -w 0 > ' + base + '.vcf', shell=True)


	for filepath in glob.glob('*.vcf'):

		a_count = 0	
		c_count = 0
		t_count = 0
		g_count = 0

		a_to_c = 0
		a_to_t = 0
		a_to_g = 0

		c_to_a = 0
		c_to_t = 0
		c_to_g = 0

		t_to_a = 0
		t_to_c = 0
		t_to_g = 0

		g_to_a = 0
		g_to_c = 0
		g_to_t = 0

		for line in open(filepath):
			# not a comment at the header
			line_split = line.split('\t')
			temp_dp = float(line_split[3])
			temp_ref = line_split[2]
			temp_a_count = int(line_split[5].split(':')[1])
			temp_c_count = int(line_split[6].split(':')[1])
			temp_g_count = int(line_split[7].split(':')[1])
			temp_t_count = int(line_split[8].split(':')[1])
			
			if temp_ref == 'A':
				a_count += temp_a_count
				if temp_t_count / temp_dp <= 0.05:
					a_to_t += temp_t_count
				if temp_c_count / temp_dp <= 0.05:
					a_to_c += temp_c_count
				if temp_g_count / temp_dp <= 0.05:
					a_to_g += temp_g_count

			elif temp_ref == 'C':
				c_count += temp_c_count
				if temp_a_count / temp_dp <= 0.05:
					c_to_a += temp_a_count
				if temp_t_count / temp_dp <= 0.05:
					c_to_t += temp_t_count
				if temp_g_count / temp_dp <= 0.05:
					c_to_g += temp_g_count

			elif temp_ref == 'G':
				g_count += temp_g_count
				if temp_a_count / temp_dp <= 0.05:
					g_to_a += temp_a_count
				if temp_c_count / temp_dp <= 0.05:
					g_to_c += temp_c_count
				if temp_t_count / temp_dp <= 0.05:
					g_to_t += temp_t_count

			elif temp_ref == 'T':
				t_count += temp_t_count
				if temp_a_count / temp_dp <= 0.05:
					t_to_a += temp_a_count
				if temp_c_count / temp_dp <= 0.05:
					t_to_c += temp_c_count
				if temp_g_count / temp_dp <= 0.05:
					t_to_g += temp_g_count


		total = a_count + c_count + t_count + g_count# + a_to_c + a_to_t + a_to_g + c_to_a + c_to_t + c_to_g + t_to_a + t_to_c + t_to_g + g_to_a + g_to_c + g_to_t
		print('Total Bases counted: ' + str(total))
		if total != 0:
			total = total / 10000.0
		else:
			total = 1

		print()
		print('A count: ' + str(a_count/total))
		print('\tA -> C: ' + str(a_to_c/total))
		print('\tA -> T: ' + str(a_to_t/total))
		print('\tA -> G: ' + str(a_to_g/total))
		print('')
		print('C count: ' + str(c_count/total))
		print('\tC -> A: ' + str(c_to_a/total))
		print('\tC -> T: ' + str(c_to_t/total))
		print('\tC -> G: ' + str(c_to_g/total))
		print('')
		print('T count: ' + str(t_count/total))
		print('\tT -> A: ' + str(t_to_a/total))
		print('\tT -> C: ' + str(t_to_c/total))
		print('\tT -> G: ' + str(t_to_g/total))
		print('')
		print('G count: ' + str(g_count/total))
		print('\tG -> A: ' + str(g_to_a/total))
		print('\tG -> C: ' + str(g_to_c/total))
		print('\tG -> T: ' + str(g_to_t/total))
		info = filepath.split('.')[0]
		print(info)
		samp = info
		con =''# info.split('_')[0]
		prep =''# info.split('_')[1]
		passage =''# info.split('_#')[1]
		g = open(base + '.csv', 'w')
		g.write(samp + ',' + prep + ',' + con + ',total,' + str(total) +','+ str(af)+ ',' + passage +'\n')

		g.write(samp + ',' + prep + ',' + con + ',A,' + str(a_count) +','+ str(af)+ ',' + passage +'\n')
		g.write(samp + ',' + prep + ',' + con + ',AC,' + str(a_to_c/total) +','+ str(af) + ',' + passage +'\n')
		g.write(samp + ',' + prep + ',' + con + ',AT,' + str(a_to_t/total) + ','+ str(af) + ',' + passage +'\n')
		g.write(samp + ',' + prep + ',' + con + ',AG,' + str(a_to_g/total) + ','+ str(af) + ',' + passage +'\n')

		g.write(samp + ',' + prep + ',' + con + ',C,' + str(c_count) +','+ str(af)+ ',' + passage +'\n')
		g.write(samp + ',' + prep + ',' + con + ',CA,' + str(c_to_a/total) + ','+ str(af) + ',' + passage +'\n')
		g.write(samp + ',' + prep + ',' + con + ',CT,' + str(c_to_t/total) + ','+ str(af) + ',' + passage +'\n')
		g.write(samp + ',' + prep + ',' + con + ',CG,' + str(c_to_g/total) + ','+ str(af) + ',' + passage +'\n')

		g.write(samp + ',' + prep + ',' + con + ',T,' + str(t_count) +','+ str(af)+ ',' + passage +'\n')
		g.write(samp + ',' + prep + ',' + con + ',TA,' + str(t_to_a/total) + ','+ str(af) + ',' + passage +'\n')
		g.write(samp + ',' + prep + ',' + con + ',TC,' + str(t_to_c/total) +','+ str(af) +  ',' + passage +'\n')
		g.write(samp + ',' + prep + ',' + con + ',TG,' + str(t_to_g/total) + ','+ str(af) + ',' + passage +'\n')

		g.write(samp + ',' + prep + ',' + con + ',G,' + str(g_count) +','+ str(af)+ ',' + passage +'\n')
		g.write(samp + ',' + prep + ',' + con + ',GA,' + str(g_to_a/total) + ','+ str(af) + ',' + passage +'\n')
		g.write(samp + ',' + prep + ',' + con + ',GC,' + str(g_to_c/total) + ','+ str(af) + ',' + passage +'\n')
		g.write(samp + ',' + prep + ',' + con + ',GT,' + str(g_to_t/total) + ','+ str(af) + ',' + passage +'\n')
		g.close()

	subprocess.call('cat *.vcf > final_output.csv', shell = True )



