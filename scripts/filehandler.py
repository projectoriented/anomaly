import os


class FamilyData(object):
    def __init__(self, common_dir):
        base_path = f"{common_dir}"
        family_id_list = get_subdirs(base_path)
        self.family = {}
        counter = 0
        for id in family_id_list:
            family_id_path = f"{common_dir}/{id}"
            sample_ids = get_subdirs(family_id_path)
            for sample_id in sample_ids:
                counter += 1
                if not sample_id == "dna_vcf":
                    sample_id_path = f"{family_id_path}/{sample_id}"
                    self.family[f"record_{counter}"] = {"family": id, "sample_id": sample_id, "fastq": f"{sample_id_path}/fastq"}

class FileHandler(object):
    def __init__(self, rna_data_dir):
        self.rna_data_dir = rna_data_dir
        self.dna_vcf_dir = None
        self.rna = FamilyData(rna_data_dir)


def get_subdirs(path):
    return os.listdir(path)


#if __name__ == '__main__':
#    filehandler = FileHandler("/home/proj/development/rare-disease/rna_data")
#    print(filehandler.rna.family)
