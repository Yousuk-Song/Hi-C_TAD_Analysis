#!/home/juyoung/anaconda3/envs/FANC/bin/python

"""
The purpose of this python3 script is to process hi-c data using fan-c 
Author: Yousuk Song
Last updated date: 2023.06.15
"""

import argparse
import datetime
import os
import fanc
import fanc.plotting as fancplot
from fanc.architecture.comparisons import hic_pca
import math
import genomic_regions as gr

def parse_arguments():
	parser = argparse.ArgumentParser(description="Input mcool and resoltuon for fanc analysis")
	required = parser.add_argument_group('required arguments')
	required.add_argument('-m', '--MCOOL',
			  required=True,
			  type=str,
			  help="Input Hi-C mcool file")
	required.add_argument('-r', '--RESOLUTION',
                          required=True,
                          type=str,
                          help="Resolution required for analysis (1kb, 2kb, 5kb, 10kb, 25kb, 100kb, 250kb, 500kb, 1mb)")
	args = parser.parse_args()
	return args


def get_current_datetime():
	return str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))


def print_log(msg):
	print("[" + get_current_datetime() + "] " + str(msg))


def mcool2cool(mcool, res):
	hic = f"{mcool}@{res}"
	return hic


def boundries(mcool, chrom, pos1, pos2, res, ins_bin, insulation_file):
	print_log(f'Calculate boundaries for {mcool} at {chrom}:{pos1}-{pos2}')
	cool =  mcool2cool(mcool, res)
	boundaries_file = f'{insulation_file}_boundaries'
	os.system(f'fanc boundaries {insulation_file} {boundaries_file} -w {ins_bin} -x')


def insulation(mcool, chrom, pos1, pos2, res, ins_bin):
	print_log(f'Calculate insulation score for {mcool} at {chrom}:{pos1}-{pos2}')
	cool =  mcool2cool(mcool, res)
	insulation_file = f'{mcool}.insulation_1mb_{chrom}:{pos1}_{pos2}'
	os.system(f'fanc insulation {cool} {insulation_file} -w {ins_bin} -r {chrom}:{pos1}-{pos2} --impute -o bed')

def plot(mcool, chrom, pos1, pos2, res, ins_bin, insulation_file, boundaries_file):
	print_log(f'Plot tad for {mcool} at {chrom}:{pos1}-{pos2}')
	final_plot = f'{mcool}_{res}_tads_insulatuon_{ins_bin}_{chrom}:{pos1}-{pos2}.png'
	cool = mcool2cool(mcool, res)
	os.system(f'fancplot -o {final_plot} \
			{chrom}:{pos1}-{pos2} \
			-p triangular {cool} -vmin 0 -vmax 0.03 \
			-p line {insulation_file} -l "{ins_bin}" \
			-p bar {boundaries_file}')

def loop(mcool, chrom, pos1, pos2, res):
	cool = mcool2cool(mcool, res)
	loop_file = f'{mcool}_{chrom}_{pos1}_{pos2}.loops'
	filtered_loop_file = loop_file.replace('.loops', '.filtered.loops')
	merged_loop_file = loop_file.replace('.loops', '.merged.loops')
	bedpe_loop_file = f'{merged_loop_file}.bedpe'
	# loop calling
	os.system(f'fanc loops {cool} {loop_file} -t 2 -f')
	# loop filter
	os.system(f'fanc loops {loop_file} {filtered_loop_file} --rh-filter -d 5 -o 5 -f')
	os.system(f'rm {loop_file}')
	# loop merge
	os.system(f'fanc loops {filtered_loop_file} {merged_loop_file} -j --remove-singlets -f')
	# export loop to bedpe
	os.system(f'fanc loops {merged_loop_file} -b {bedpe_loop_file} -f') 
	os.system(f'rm {merged_loop_file}')

def compartment(mcool, ref):
	os.system(f'fanc compartments {mcool}@1mb {mcool}_1mb.ab -f')
	os.system(f'fanc compartments -g {ref} -d {mcool}.domains.bed {mcool}_1mb.ab -f')

def len2mb(length):
	return f'{math.ceil(length/1000000)}mb'

chromosome_lengths = {
    'chr1': 249250621,
    'chr2': 243199373,
    'chr3': 198022430,
    'chr4': 191154276,
    'chr5': 180915260,
    'chr6': 171115067,
    'chr7': 159138663,
    'chr8': 146364022,
    'chr9': 141213431,
    'chr10': 135534747,
    'chr11': 135006516,
    'chr12': 133851895,
    'chr13': 115169878,
    'chr14': 107349540,
    'chr15': 102531392,
    'chr16': 90354753,
    'chr17': 81195210,
    'chr18': 78077248,
    'chr19': 59128983,
    'chr20': 63025520,
    'chr21': 48129895,
    'chr22': 51304566,
    'chrX': 155270560,
    'chrY': 59373566
}

if __name__ == '__main__':
	args = parse_arguments()
	mcool = args.MCOOL
	res = args.RESOLUTION
	ref = '/mnt/lustre/tmp_juyoung/juyoung/ref/hg19_HPV/hg19_HPV.fa'
	gtf = '/home/juyoung/Integration/KNCC_Hi-C/Hi-C_Script/TAD_Calling/gencode.v43lift37.basic.annotation.protein_coding.gff3'
	ins_bin = '1mb'
	compartment(mcool, ref)
	for chrom in chromosome_lengths:
		length = chromosome_lengths[chrom]
		pos1 = '0mb'
		pos2 = len2mb(length)
		insulation(mcool, chrom, pos1, pos2, res, ins_bin) # insulation calling
		insulation_file = f'{mcool}.insulation_1mb_{chrom}:{pos1}_{pos2}_{ins_bin}.bed'
		boundries(mcool, chrom, pos1, pos2, res, ins_bin, insulation_file) # boundaries calling
		boundaries_file = f'{insulation_file}_boundaries'
		plot(mcool, chrom, pos1, pos2, res, ins_bin, insulation_file, boundaries_file)
#		loop(mcool, chrom, pos1, pos2, res)








