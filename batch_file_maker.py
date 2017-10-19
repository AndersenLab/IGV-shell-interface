import pandas as pd
import os

#Have a QTL, as part of CEGWAS output they do fine mapping, look at most correlated variants within QTL, have confidence interval. Has a ref and alt group. So that would be regions across multiple strains. 
#Other use case validate VCF outputs from various variant callers. Pull out indel calls of various sizes, pick one strain that is ref and one that is alt, refs are 0/0, alts are 1/1 and hets are 0/1. Shouldn't be any het calls. Dan did heterozygous polarization, push their genotype into the closest homozygous genome. 
#Input: VCF, # of variants to take pictures, flag => needs alt, ref, and heterozygous (might be none) call.
#Add metadata into the filenames
#VCF file formats are supposed to be standardized.


def import_testfile(testfile = 'for_garth.vcf'): # Reads in a VCF file and extracts the positions to snapshot
	df = pd.read_table(testfile, skiprows=57)
	df = df.rename(columns={'#CHROM':'CHROM'})
	df['END'] = df['INFO'].map(lambda x: x.split(';')[4][4:])
	df = df.astype(str)
	return(df)

def strain_module(strain = False): # Finds and loads a .bam and .bai file from the elegansvariation database
	if not strain:
		strain = input('Strain (e.g. N2): ')
	df = pd.read_table('https://www.elegansvariation.org/data/download/bam.sh', skiprows=7, header=None, sep=' ')
	df = df.ix[:,list(~df.isnull().all())]
	bam = df[df[2] == strain + '.bam'][11].values[0]
	bai = df[df[2] == strain + '.bam.bai'][7].values[0]
	result = ['load ' + bam + ' index=' + bai + ' name=' +strain]
	return(result)

def goto_module(chromosome=False, start=False, end=False): # Makes a list of the positions to visit
	if not chromosome:
		chromosome = input('Chromosome (e.g. IV): ')
	if not start:
		start = input('start bp (e.g. 1): ')
	if not end:
		end = input('end bp (e.g. 100): ')
	goto = ['goto chr' + chromosome + ':' + start + '-' + end, 'snapshot']
	return(goto)

def genome_module(genome, snapshot_directory, new=True): # Loads the genome and sets up the batch file, assuming one genome per batch file for now.
	batch = []
	if new:
		batch = batch + ['new']
	if snapshot_directory:
		batch = batch + ['snapshotDirectory '+snapshot_directory]
	if genome:
		batch += ['genome '+genome]
	return(batch)

def construct(testfile='for_garth.vcf', genome='WS245', snapshot_directory=False, stack_tracks=True): # Brings all the modules together into a batch file
	if not snapshot_directory:
		snapshot_directory = os.getcwd()
	df = import_testfile(testfile)
	strains = [s for s in df.columns if s not in ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'END']]
	gotos = df.apply(lambda x: goto_module(chromosome=x['CHROM'], start=x['POS'], end=x['END']), axis = 1).sum()
	batch = genome_module(genome, snapshot_directory)

	if not stack_tracks:
		for strain in strains:
			batch += strain_module(strain)
			batch += gotos
	elif stack_tracks:
		for strain in strains:
			batch += strain_module(strain)
		batch += gotos

	batch = '\n'.join(batch)
	file = open('testfile.txt','w') 
	file.write(batch) 
	file.close()

construct()

'''
Example:
new
snapshotDirectory Desktop
genome WS245
load https://northwestern.box.com/shared/static/5wtb0enhny1ypq8budpi0gczq2ww2p5z.bam index=https://northwestern.box.com/shared/static/inofdwap2coam2a44qg9t0zp4nam5jso.bai
goto chrI:11042686-11043034
snapshot
goto chrI:11051345-11051565
snapshot
goto chrI:11061918-11075187
snapshot
goto chrI:11077078-11077266
snapshot
goto chrI:11100674-11100746
snapshot
goto chrI:11141313-11141373
snapshot
goto chrII:12635966-12635987
snapshot
goto chrII:12636749-12636857
snapshot
goto chrII:12638697-12639091
snapshot
goto chrII:12640249-12640363
snapshot
goto chrII:12647409-12647435
snapshot
goto chrII:12648860-12651963
snapshot
goto chrIV:16132880-16133330
snapshot
goto chrIV:16139857-16139879
snapshot
goto chrIV:16156718-16156745
snapshot
goto chrIV:16165366-16165468
snapshot
goto chrIV:16186430-16186496
snapshot
goto chrIV:16187759-16187916
snapshot
goto chrIV:16203511-16203633
snapshot
goto chrIV:16225780-16248687
snapshot
goto chrV:7795697-7795717
snapshot
goto chrV:7815340-7815379
snapshot
goto chrV:7885228-7885460
snapshot
goto chrV:7948325-7948369
snapshot
goto chrV:7948602-7948713
snapshot
goto chrX:14188407-14189353
snapshot
goto chrX:14193282-14193332
snapshot
goto chrX:14218719-14218753
snapshot
goto chrX:14224822-14224847
snapshot
'''