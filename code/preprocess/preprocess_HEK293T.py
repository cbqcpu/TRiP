# Pre-process of polysome sequencing data
# Cut adapter, mapping human genome, mapping yeast genome
# CBQ 2023.3.28

import os

def rename_file_and_cut_adapter(renamelist_file, input_dir, output_dir):
    f = open(renamelist_file, "r")
    renamelist = f.readlines()
    f.close()
    # renamelist format, csv:
    # first line: origin_name,new_name
    # origin_name is the file name from sequencing company with prefix, e.g. CBQ20230309_NC2h-1-1_1 prifix is *_NC2h-1-1
    # new_name format: Condition-replication-polysome frection
    os.system("mkdir -p "+output_dir+"_clean")
    os.system("mkdir -p "+output_dir)

    for i in renamelist[1:]:
        filenames = i.strip().split(",")
        #cat_cmd_R1 = "cat "+input_dir+"/"+filenames[0]+"_1.fq.gz >> "+output_dir+"/"+filenames[1]+"_1.fq.gz"
        #cat_cmd_R2 = "cat "+input_dir+"/"+filenames[0]+"_2.fq.gz >> "+output_dir+"/"+filenames[1]+"_2.fq.gz"
        #print(cat_cmd_R1)
        #print(cat_cmd_R2)
        #os.system(cat_cmd_R1)
        #os.system(cat_cmd_R2)
        filenames = i.strip().split(",")
        cutadapt_cmd = "cutadapt -j 2 -a CTGTCTCTTATA -A CTGTCTCTTATA -q 35 -m 15 -u 4 -u -4 -o "+output_dir+"_clean/"+filenames[1]+"_1.fq.gz "+"-p "+output_dir+"_clean/"+filenames[1]+"_2.fq.gz "+output_dir+"/"+filenames[1]+"_1.fq.gz "+output_dir+"/"+filenames[1]+"_2.fq.gz"
        print(cutadapt_cmd)
        os.system(cutadapt_cmd)
        merge_cmd = f"cat {output_dir}_clean/{filenames[1]}_*.fq.gz >> {output_dir}_clean/{filenames[1]}.merge.fq.gz"
        print(merge_cmd)
        os.system(merge_cmd)

def mapping_to_human_genome(renamelist_file, hg_file, anotation_bed, raw_clean, mapping_out):
    f = open(renamelist_file, "r")
    renamelist = f.readlines()
    f.close()
    f = open("data/data_list_slamdunk.csv", "w")
    for i in renamelist[1:]:
        filenames = i.strip().split(",")
        samp = raw_clean+"/"+filenames[1]+".merge.fq.gz,"+"SampleName"+",Condition"+",2\n"
        f.write(samp)
    f.close()
    mapping_cmd = "slamdunk all -r "+hg_file+" -b "+anotation_bed+" -o "+mapping_out+" -5 12 -n 100 -t 48 --max-read-length 160 "+"data/data_list_slamdunk.csv"
    print(mapping_cmd)
    os.system(mapping_cmd)
    
def mapping_to_yeast_genome(renamelist_file, yg_file, raw_clean, mapping_out):
    f = open(renamelist_file, "r")
    renamelist = f.readlines()
    f.close()
    f = open("data/datayeast_mapping_list.txt", "w")
    for i in renamelist[1:]:
        filenames = i.strip().split(",")
        f.write(filenames[1]+"\n")
    f.close()
    mapping_cmd = "code/mapping_to_yeast.sh "+raw_clean+" "+yg_file+" data/datayeast_mapping_list.txt "+mapping_out
    print(mapping_cmd)
    os.system(mapping_cmd)



import os
import subprocess

def count_reads_bam(bam_file_path, gtf_file_path, output_file_path):
    """
    This function uses featureCounts from the subread package to count reads in a BAM file.

    :param bam_file_path: str, full path to input BAM file
    :param gtf_file_path: str, full path to GTF annotation file
    :param output_file_path: str, full path to output file
    :return: None
    """

    # Check if BAM file exists
    if not os.path.isfile(bam_file_path):
        raise FileNotFoundError(f"BAM file not found: {bam_file_path}")

    # Check if GTF file exists
    if not os.path.isfile(gtf_file_path):
        raise FileNotFoundError(f"GTF file not found: {gtf_file_path}")

    # Run featureCounts
    cmd = ['featureCounts', '-a', gtf_file_path, '-T', '20', '-o', output_file_path, bam_file_path]
    subprocess.run(cmd, check=True)

# Example usage:
# count_reads_bam('/path/to/featureCounts', '/path/to/input.bam', '/path/to/annotation.gtf', '/path/to/output.txt')


import os
import pandas as pd

def process_bam_files(csv_file_path, bam_files_dir, gtf_file_path, output_dir):
    """
    This function processes multiple BAM files based on a CSV file with new names for each file.

    :param csv_file_path: str, full path to CSV file with old and new names for each BAM file
    :param bam_files_dir: str, directory where BAM files are located
    :param gtf_file_path: str, full path to GTF annotation file
    :param output_dir: str, directory where to save output files
    :return: None
    """

    # Load CSV file with pandas
    df = pd.read_csv(csv_file_path)

    # Initialize dataframe for concatenated results
    all_counts_df = pd.DataFrame()

    # Process each BAM file
    for _, row in df.iterrows():
        bam_file_name = row['new_name'] + ".sort.rmdup.bam"
        bam_file_path = os.path.join(bam_files_dir, bam_file_name)
        output_file_path = os.path.join(output_dir, row['new_name'] + "_counts.txt")

        # Count reads
        count_reads_bam(bam_file_path, gtf_file_path, output_file_path)

        # Load result into dataframe
        try:
            counts_df = pd.read_csv(output_file_path, sep="\t", skiprows=1, usecols=[0, 6], index_col=0, names=['GeneID', 'Counts'])
        except IndexError:
            print(f"Error: Unable to read file {output_file_path}, skipping.")
            continue

        # Rename column to sample name and concatenate
        counts_df.columns = [row['new_name']]
        all_counts_df = pd.concat([all_counts_df, counts_df], axis=1)

    # Save concatenated result
    all_counts_df.to_csv(os.path.join(output_dir, "all_counts.csv"))


# Example usage:
# process_bam_files('/path/to/bam_prefixes.csv', '/path/to/bam_files', '/path/to/annotation.gtf', '/path/to/output')


#30 min
rename_file_and_cut_adapter("data/renamelist_NC30min.csv", "/storage/cbq/sequencing_data/CBQ/Polysome_Sequencing/*/*/01.RawData/*/", "data/raw_NC30min")
mapping_to_human_genome("data/renamelist_NC30min.csv", "data/refs/GRCh38/GRCh38.primary_assembly.genome.fa", "data/refs/GRCh38/gencode.v35.annotation.bed","data/raw_NC30min_clean", "data/NC30min_mapping_hg")
mapping_to_yeast_genome("data/renamelist_NC30min.csv", "data/refs/saccharomyces_cerevisiae/R64_hisat2", "data/raw_NC30min_clean", "data/NC30min_mapping_yeast")
process_bam_files('data/renamelist_NC30min.csv', 'data/NC30min_mapping_yeast', 'data/refs/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.109.gtf', 'data/NC30min_mapping_yeast')

#1h
rename_file_and_cut_adapter("data/renamelist_NC1h.csv", "/storage/cbq/sequencing_data/CBQ/Polysome_Sequencing/*/*/01.RawData/*/", "data/raw_NC1h")
mapping_to_human_genome("data/renamelist_NC1h.csv", "data/refs/GRCh38/GRCh38.primary_assembly.genome.fa", "data/refs/GRCh38/gencode.v35.annotation.bed","data/raw_NC1h_clean", "data/NC1h_mapping_hg")
mapping_to_yeast_genome("data/renamelist_NC1h.csv", "data/refs/saccharomyces_cerevisiae/R64_hisat2", "data/raw_NC1h_clean", "data/NC1h_mapping_yeast")
process_bam_files('data/renamelist_NC1h.csv', 'data/NC1h_mapping_yeast', 'data/refs/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.109.gtf', 'data/NC1h_mapping_yeast')

#2h
rename_file_and_cut_adapter("data/renamelist_NC2h_lcd.csv", "/storage/cbq/sequencing_data/CBQ/Polysome_Sequencing/20230529_NC_2h/CP2023020200044/H101SC23030069/RSSQ00504/X101SC23030069-Z01/X101SC23030069-Z01-J117/01.RawData/*/", "data/raw_NC2h_lcd")
mapping_to_human_genome("data/renamelist_NC2h_lcd.csv", "data/refs/GRCh38/GRCh38.primary_assembly.genome.fa", "data/refs/GRCh38/gencode.v35.annotation.bed","data/raw_NC2h_lcd_clean", "data/NC2h_lcd_mapping_hg")
mapping_to_yeast_genome("data/renamelist_NC2h_lcd.csv", "data/refs/saccharomyces_cerevisiae/R64_hisat2", "data/raw_NC2h_lcd_clean", "data/NC2h_lcd_mapping_yeast")
process_bam_files('data/renamelist_NC2h_lcd.csv', 'data/NC2h_lcd_mapping_yeast', 'data/refs/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.109.gtf', 'data/NC2h_lcd_mapping_yeast')

#Puro 2h dose and labeling
rename_file_and_cut_adapter("data/renamelist_Purox2h.csv", "/storage/cbq/sequencing_data/CBQ/Polysome_Sequencing/20230328_STM2457_12h_puro_2h/CP2023020200044/H101SC23030069/RSSQ00504/X101SC23030069-Z01/X101SC23030069-Z01-J028/01.RawData/*/", "data/raw_Purox2h")
mapping_to_human_genome("data/renamelist_Purox2h.csv", "data/refs/GRCh38/GRCh38.primary_assembly.genome.fa", "data/refs/GRCh38/gencode.v35.annotation.bed","data/raw_Purox2h_clean", "data/Purox2h_mapping_hg")
mapping_to_yeast_genome("data/renamelist_Purox2h.csv", "data/refs/saccharomyces_cerevisiae/R64_hisat2", "data/raw_Purox2h_clean", "data/Purox2h_mapping_yeast")
process_bam_files('data/renamelist_Purox2h.csv', 'data/Purox2h_mapping_yeast', 'data/refs/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.109.gtf', 'data/Purox2h_mapping_yeast')

#STM2457 20uM 2h
rename_file_and_cut_adapter("data/renamelist_STM2457x2h.csv", "/storage/cbq/sequencing_data/CBQ/Polysome_Sequencing/20230408_WT_STM2457_20uM_12h_2h/CP2023020200044/H101SC23030069/RSSQ00504/X101SC23030069-Z01/X101SC23030069-Z01-J046/01.RawData/*/", "data/raw_STM2457x2h")
mapping_to_human_genome("data/renamelist_STM2457x2h.csv", "data/refs/GRCh38/GRCh38.primary_assembly.genome.fa", "data/refs/GRCh38/gencode.v35.annotation.bed","data/raw_STM2457x2h_clean", "data/STM2457x2h_mapping_hg")
mapping_to_yeast_genome("data/renamelist_STM2457x2h.csv", "data/refs/saccharomyces_cerevisiae/R64_hisat2", "data/raw_STM2457x2h_clean", "data/STM2457x2h_mapping_yeast")
process_bam_files('data/renamelist_STM2457x2h.csv', 'data/STM2457x2h_mapping_yeast', 'data/refs/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.109.gtf', 'data/STM2457x2h_mapping_yeast')

#YTHDFs
rename_file_and_cut_adapter("data/renamelist_YTHDFs.csv", "/storage/cbq/sequencing_data/CBQ/Polysome_Sequencing/20230606_YTHDF123/CP2022041300022/H101SC23023451/RSSQ00504/X101SC23023451-Z01/X101SC23023451-Z01-J131/01.RawData/*/", "data/raw_YTHDFs")
mapping_to_human_genome("data/renamelist_YTHDFs.csv", "data/refs/GRCh38/GRCh38.primary_assembly.genome.fa", "data/refs/GRCh38/gencode.v35.annotation.bed","data/raw_YTHDFs_clean", "data/YTHDFs_mapping_hg")
mapping_to_yeast_genome("data/renamelist_YTHDFs.csv", "data/refs/saccharomyces_cerevisiae/R64_hisat2", "data/raw_YTHDFs_clean", "data/YTHDFs_mapping_yeast")
process_bam_files('data/renamelist_YTHDFs.csv', 'data/YTHDFs_mapping_yeast', 'data/refs/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.109.gtf', 'data/YTHDFs_mapping_yeast')

