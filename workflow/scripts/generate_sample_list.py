from os import walk, path
from sys import argv

"""
This python script takes 1 argument from the command-line which is the absolute common path to all Fastq files. 
The output is a tabular file of sample names, lane, ILLUMINA sample number, ID, platform, PU, path of fq1, and fq2. 
"""
def get_fastq_path(common_dir):
    file_paths = []
    for folderName, subfolders, filenames in walk(common_dir):
        for filename in filenames:
            if "fastq.gz" in filename:
                file_paths.append(f"{folderName}/{filename}")
    return file_paths


def tabulate_data(files):
    fq1 = [f.replace("_1.fastq.gz", "") for f in files if "_1.fastq.gz" in f]
    fq2 = [f.replace("_2.fastq.gz", "") for f in files if "_2.fastq.gz" in f]

    samples = [f"sample\tlane\tsample_number\tid\tplatform\tpu\tfq1\tfq2"]
    for i in range(len(fq1)):
        lane = path.basename(fq1[i])[0]
        num_identifier = path.basename(fq1[i]).split("_")[1]
        h_identifier = path.basename(fq1[i]).split("_")[2]
        sample_id = path.basename(fq1[i]).split("_")[3]
        sample_num = path.basename(fq1[i]).split("_")[4]
        pu = f"{h_identifier}.{lane}.{sample_num}"
        platform = "ILLUMINA"
        id = f"{sample_id}_{num_identifier}_{h_identifier}_{sample_num}_lane{lane}"
        for j in range(len(fq2)):
            if fq1[i] == fq2[j]:
                samples.append(f"{sample_id}\t{lane}\t{sample_num}\t{id}\t{platform}\t{pu}\t{fq1[i]}_1.fastq.gz\t{fq2[j]}_2.fastq.gz")
    return samples


if __name__ == '__main__':
    samples = tabulate_data(get_fastq_path(argv[1]))
    with open("samples.tsv", "w") as tsv_file:
        for element in samples:
            tsv_file.write(f"{element}\n")
    if path.isfile("samples.tsv"):
        abs_path = path.abspath("samples.tsv")
        print(f"File written: {abs_path}")
    else:
        print("File didn't write")