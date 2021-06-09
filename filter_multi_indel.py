import vcf
import itertools
import numpy as np
import argparse
#import pdb

def read_vcf(file_path):
    return vcf.Reader(open(file_path, 'r'))

def check_for_indels(vcfs):
    return [(vcfs[i].POS, vcfs[i].ALT, vcfs[i].INFO) for i,v in enumerate(vcfs) if "IMF" in vcfs[i].INFO]

def find_min_multiIMF(indels):
    combs = list(itertools.combinations([indels[i][0] for i in range(len(indels))], 2))
    duplicates=np.unique([x[0] for x in combs if x[0]==x[1]])
    dup_dict={x:[] for x in duplicates}
    for i in indels:
        if i[0] in dup_dict:
            dup_dict[i[0]].append(i[2]['AD'])
    return {k:min(v) for k,v in dup_dict.items()}

def construct_filtered_vcf(vcfs, indels_to_filter):
    new_vcf = []    
    for record in vcfs:
        if record.POS not in indels_to_filter.keys():
            new_vcf.append(record)
        if record.POS in indels_to_filter.keys() and record.INFO['AD'] not in indels_to_filter.values():
            new_vcf.append(record)
        else:
            continue
    return new_vcf

def write_new_vcf(filtered_vcf, vcf_reader, vcf_file_path):
    #vcf_reader = vcf.Reader(filename=vcf_filename)
    vcf_writer = vcf.Writer(open(vcf_file_path, 'w'), vcf_reader)
    for record in filtered_vcf:
        #print(record)
        vcf_writer.write_record(record)

def lines_that_start_with(string, fp):
    return [line for line in fp if line.startswith(string)]

def add_metadata(vcf_file_path, indels_to_filter):
    a_file = open(vcf_file_path, "r")
    list_of_lines = a_file.readlines()
    metadata_header=len([line for line in lines_that_start_with("##", list_of_lines)])
    for indel_position, indel_ad in indels_to_filter.items():
        list_of_lines.insert(metadata_header, "##Python filtered an INDEL at POS=" + str(indel_position)+" with and AD="+ str(indel_ad) +"\n")
        metadata_header+=1

    with open(vcf_file_path, "w") as f:
        list_of_lines = "".join(list_of_lines)
        f.write(list_of_lines)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('vcf_file_to_filter', help='')
    args = parser.parse_args()  
    # Grabs name of file.
    vcf_file_path=str(args.vcf_file_to_filter)
    vcf_reader = read_vcf(vcf_file_path)
    vcf_list=[record for record in vcf_reader]
    indels = check_for_indels(vcf_list)
    if indels:
        indels_to_filter=find_min_multiIMF(indels)
        filtered_vcf=construct_filtered_vcf(vcf_list, indels_to_filter)
        write_new_vcf(filtered_vcf, vcf_reader, vcf_file_path)
        add_metadata(vcf_file_path, indels_to_filter)



