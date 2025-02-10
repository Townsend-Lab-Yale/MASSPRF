#!/usr/bin/python3

import sys
import os
import re
from os.path import isfile, join
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# get list of files in directory #
wd = os.getcwd()
file_list = [f for f in os.listdir(wd) if isfile(join(wd, f))]

# If necessary, make some directories for the processed results and plots

proc_dir = './processed_output/'
plot_dir = './processed_output/plots/'

if not os.path.exists(proc_dir):
	os.makedirs(proc_dir)
if not os.path.exists(plot_dir):
	os.makedirs(plot_dir)

# make lists to keep track of your genes 
boring_l = []
neg_l = []
pos_l = []
failed_l = []
str_pos_l = []

# Iterate through each file in the directory, checking to see if each run succeeded,
# and if so, whether it contains site showing positive or negative selection. 
for entry in file_list:
	if ".txt" not in entry:
		continue
	
	name = entry.split('_')
	gene_name = os.path.splitext(entry)[0]
	with open(entry) as f_in:
		if 'Mission accomplished' in f_in.read():
			# Start again from the beginning of the file
			f_in.seek(0)
			lines = f_in.read().splitlines()
			# Get scaling factor from MASS-PRF output to de-scale results
			for line in lines:
				if "computed scale factor" in line:
					#print(line)
					scaling_factor = int(line.split("computed scale factor is ")[-1].strip())
					#scaling_factor = int(re.search(r'(?<=computed scale factor is )\w+', line).group(0))
			# Pull out the part of the results file that we care about
			intro = '//Results based on model averaging of gamma using different models: '
			intro_index = lines.index(intro)
			end = 'Abbreviation: MS=Model Selection; MA=Model Averaging; CI=Confidence Interval; ps=Polymorphism Synonymous; pr=Polymorphism Replacement; ds=Divergence Synonymous; dr=Divergence Replacement; Gamma=N*s (Gamma: scaled selection coefficient (selection intensity); N: haploid effective population size; s: selection coefficient); NI=Neutrality Index (NI=(pr/dr)/(ps/ds), >1 Negative selection, <1 Positive selection); INF=Infinite; N-INF=Negative Infinite;'
			end_index = lines.index(end)
			chart_lines = lines[intro_index + 1: end_index -1]
			pd_l = []
			[pd_l.append(item.split()) for item in chart_lines]
			# Create a pd dataframe for the results
			scaled_sel_df = pd.DataFrame(pd_l[1:], columns =[i for i in pd_l[0]]) 


			scaled_sel_df['Lower_CI_Gamma'] = pd.to_numeric(scaled_sel_df['Lower_CI_Gamma'], errors='coerce')
			scaled_sel_df['Gamma'] = pd.to_numeric(scaled_sel_df['Gamma'], errors='coerce')
			scaled_sel_df['Upper_CI_Gamma'] = pd.to_numeric(scaled_sel_df['Upper_CI_Gamma'], errors='coerce')
			scaled_pos_df = scaled_sel_df[(scaled_sel_df.Gamma > 4) & (scaled_sel_df.Lower_CI_Gamma > 0)]
			scaled_str_pos_df = scaled_sel_df[(scaled_sel_df.Gamma > 4) & (scaled_sel_df.Lower_CI_Gamma > 4)]
			scaled_neg_df = scaled_sel_df[(scaled_sel_df.Gamma < -1) & (scaled_sel_df.Upper_CI_Gamma < 0)]
			if scaled_pos_df.empty:
				if scaled_neg_df.empty:
					boring_l.append(gene_name)
					print("Boring: " + gene_name)
				else:
					# The results need to be de-scaled to match with actual positions along a gene
					# Here we de-scale completely to the nucleotide position (result_position*scaling_factor),
					# though you can adjust this to scale to amino acid position ((result_positon*scaling_factor)/3)

					sel_df = pd.DataFrame(columns=scaled_sel_df.columns)
					posit = 0
					for index, row in scaled_sel_df.iterrows():
						for i in range(0, scaling_factor):
							sel_df.loc[posit] = row
							posit += 1
					sel_df.drop('Position', axis=1, inplace=True)
					position_list = list(range(1,len(sel_df)+1))
					sel_df.insert(0, 'Position', position_list)
					pos_df = sel_df[(sel_df.Gamma > 4) & (sel_df.Lower_CI_Gamma > 0)]
					str_pos_df = sel_df[(sel_df.Gamma > 4) & (sel_df.Lower_CI_Gamma > 4)]
					neg_df = sel_df[(sel_df.Gamma < -1) & (sel_df.Upper_CI_Gamma < 0)]

					## If we got here, pos_df and neg_df can't both be empty so

					if pos_df.empty:
						neg_l.append(name[0])
						neg_sites = list(neg_df['Position'])
						neg_sites_st = [str(int) for int in neg_sites]
						sel_df.to_csv(proc_dir+gene_name+".csv", index=False)
						neg_f = open(proc_dir+gene_name+"_neg_sites.txt", "w")
						neg_f.write('\n'.join(neg_sites_st) + '\n')
						neg_f.close()

					else:
						pos_l.append(name[0])
						pos_sites = list(pos_df['Position'])
						pos_sites_st = [str(int) for int in pos_sites]
						pos_f = open(proc_dir+gene_name+"_pos_sites.txt", "w")
						pos_f.write('\n'.join(pos_sites_st) + '\n')
						pos_f.close()
						if not str_pos_df.empty:
							str_pos_l.append(name[0])
							str_pos_sites = list(str_pos_df['Position'])
							str_pos_sites_st = [str(int) for int in str_pos_sites]
							str_pos_f = open(proc_dir+gene_name+"_str_pos_sites.txt", "w")
							str_pos_f.write('\n'.join(str_pos_sites_st) + '\n')
							str_pos_f.close()
						if neg_df.empty:
							sel_df.to_csv(proc_dir+gene_name+".csv", index=False)
						else:
							neg_l.append(name[0])
							neg_sites = list(neg_df['Position'])
							neg_sites_st = [str(int) for int in neg_sites]
							sel_df.to_csv(proc_dir+gene_name+".csv", index=False)
							neg_f = open(proc_dir+gene_name+"_neg_sites.txt", "w")
							neg_f.write('\n'.join(neg_sites_st) + '\n')
							
			else:
				# The results need to be de-scaled to match with actual positions along a gene
				# Here we de-scale completely to the nucleotide position (result_position*scaling_factor),
				# though you can adjust this to scale to amino acid position ((result_positon*scaling_factor)/3)

				sel_df = pd.DataFrame(columns=scaled_sel_df.columns)
				posit = 0
				for index, row in scaled_sel_df.iterrows():
					for i in range(0, scaling_factor):
						sel_df.loc[posit] = row
						posit += 1
				sel_df.drop('Position', axis=1, inplace=True)
				position_list = list(range(1,len(sel_df)+1))
				sel_df.insert(0, 'Position', position_list)
				pos_df = sel_df[(sel_df.Gamma > 4) & (sel_df.Lower_CI_Gamma > 0)]
				str_pos_df = sel_df[(sel_df.Gamma > 4) & (sel_df.Lower_CI_Gamma > 4)]
				neg_df = sel_df[(sel_df.Gamma < -1) & (sel_df.Upper_CI_Gamma < 0)]

				## If we got here, pos_df and neg_df can't both be empty so

				if pos_df.empty:
					neg_l.append(name[0])
					neg_sites = list(neg_df['Position'])
					neg_sites_st = [str(int) for int in neg_sites]
					sel_df.to_csv(proc_dir+gene_name+".csv", index=False)
					neg_f = open(proc_dir+gene_name+"_neg_sites.txt", "w")
					neg_f.write('\n'.join(neg_sites_st) + '\n')
					neg_f.close()

				else:
					pos_l.append(name[0])
					pos_sites = list(pos_df['Position'])
					pos_sites_st = [str(int) for int in pos_sites]
					pos_f = open(proc_dir+gene_name+"_pos_sites.txt", "w")
					pos_f.write('\n'.join(pos_sites_st) + '\n')
					pos_f.close()
					if not str_pos_df.empty:
						str_pos_l.append(name[0])
						str_pos_sites = list(str_pos_df['Position'])
						str_pos_sites_st = [str(int) for int in str_pos_sites]
						str_pos_f = open(proc_dir+gene_name+"_str_pos_sites.txt", "w")
						str_pos_f.write('\n'.join(str_pos_sites_st) + '\n')
						str_pos_f.close()
					if neg_df.empty:
						sel_df.to_csv(proc_dir+gene_name+".csv", index=False)
					else:
						neg_l.append(name[0])
						neg_sites = list(neg_df['Position'])
						neg_sites_st = [str(int) for int in neg_sites]
						sel_df.to_csv(proc_dir+gene_name+".csv", index=False)
						neg_f = open(proc_dir+gene_name+"_neg_sites.txt", "w")
						neg_f.write('\n'.join(neg_sites_st) + '\n')
		
		else:
			print("Failed: " + gene_name)
			failed_l.append(gene_name)

## Plot your significant results!

results_list = [f for f in os.listdir(proc_dir) if isfile(join(proc_dir, f))]

for entry in results_list:
	if entry.endswith(".csv"):
		name = entry.split('.')[0]
		data = pd.read_csv(proc_dir + entry)
		plt.plot(data['Position'],data['Gamma'])
		plt.fill_between(data['Position'],data['Lower_CI_Gamma'],data['Upper_CI_Gamma'], alpha=0.2)
		plt.axhline(y=0, color='black', linestyle='-')
		plt.xlabel('Nucleotide position')
		plt.ylabel('Scaled selection coefficient (gamma)')
		plt.savefig(plot_dir + name + '.pdf') 
		plt.close()


neg_f = open(proc_dir+"Negative_genes.txt", "w")
neg_f.write('\n'.join(neg_l) + '\n')
neg_f.close()

pos_f = open(proc_dir+"Positive_genes.txt", "w")
pos_f.write('\n'.join(pos_l) + '\n')
pos_f.close()

str_pos_f = open(proc_dir+"Strongly_positive_genes.txt", "w")
str_pos_f.write('\n'.join(str_pos_l) + '\n')
str_pos_f.close()

boring_f = open(proc_dir+"Boring_genes.txt", "w")
boring_f.write('\n'.join(boring_l) + '\n')
boring_f.close()

failed_f = open(proc_dir+"Failed_genes.txt", "w")
failed_f.write('\n'.join(failed_l) + '\n')
failed_f.close()


