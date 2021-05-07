import subprocess 
import argparse
import sys

def find_values(line):
	idv = line.split("IDV=")[1].split(";")[0]
	dp = line.split("DP=")[1].split(";")[0]
	dp4 = line.split("DP4=")[1].split(";")[0]
	ad = line.split("AD=")[1].split(";")[0]
	return idv,dp,dp4,ad


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-file', help="Provide pre_bcftools vcf.")
	
	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)
	file = args.file

	with open(file) as f:
		content = f.readlines()
	content = [x.strip() for x in content]

	sample_name = file.split("_pre_b")[0]

	new_vcf = open(sample_name+"_pre1.vcf","w+") 

	list_of_indel_positions=[]
	list_of_indels = []
	indel_replaced = []

	for index, line in enumerate(content):
		# Let's correct some indels
		if "#" not in line and "INDEL" in line:
			position=int(line.split("\t")[1])

			ref = line.split("\t")[3]
			ref_length=len(ref)
			indel = (line.split("\t")[4])
			indel_length=len(indel)

			# Fix deletions
			if(indel_length < ref_length and indel_length < 500):
				# Grab the last nearest depth of non-indel
				index_num=index-1
				while("INDEL" in content[index_num] or int(content[index_num].split("\t")[1])>position):
					index_num -= 1
				prev_dp=int((content[index_num-1]).split("DP=")[1].split(";")[0])

				# Grab the next nearest depth of non-indel
				index_num = index + 1
				while(int(content[index_num].split("\t")[1])<(position+ref_length) or "INDEL" in content[index_num]):
					index_num+=1
				next_dp=int((content[index_num+1]).split("DP=")[1].split(";")[0])

				max_dp = max(prev_dp, next_dp)

				idv,dp,dp4,ad = find_values(line)
				
				# Trust the AD for the indel depth
				indel_depth = ad.split(",")[1]

				# Compare to depth around the position
				replacement_ref = max_dp - int(indel_depth)

				# Check IDV compared to local depth
				if(float(idv)>=(max_dp/2)):
					# Check if we have an indel at the same position
					if(position not in list_of_indel_positions):
						list_of_indel_positions.append(position)
						list_of_indels.append(line)

					# Use AD for these these since IDV sucks with indels at same position
					else:
						index = list_of_indel_positions.index(position)
						comparison_indel = list_of_indels[index]
						comparison_ad = comparison_indel.split("AD=")[1].split(";")[0].split(",")[1]

						if(int(indel_depth) >= int(comparison_ad)):
							list_of_indels[index] = line
			else:
				new_vcf.write(line+"\n")
		else:
			new_vcf.write(line+"\n")

	for indel in list_of_indels:
		idv,dp,dp4,ad = find_values(indel)
		indel_depth_half = int(int(indel_depth) / 2)
		ref_depth_half = int(replacement_ref / 2)

		new_line = indel.replace(dp4,str(ref_depth_half)+","+str(ref_depth_half)+","+str(indel_depth_half)+","+str(indel_depth_half))
		new_line = indel.replace(ad,str(replacement_ref)+","+str(indel_depth))

		new_vcf.write(new_line+"\n")