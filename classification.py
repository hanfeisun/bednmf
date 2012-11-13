from numpy import array, nonzero, empty, hstack, ones
from scipy.sparse import csr_matrix
from scipy.spatial.distance import cosine
from bed import read_bed, divide_bed, merge_bed, bed2vector, read_chromosome_len, coverage_bed, extract_vec
from nmf import factorization
from heapq import nsmallest

class BedNMF:
    def __init__(self, beds_path, chromosome_len_path, dhs_path, bin_size):
        self.beds_path = beds_path
        self.bin_size = bin_size
        self.dhs_path = dhs_path

        with open(chromosome_len_path) as clp:
            self.chromosome_len = read_chromosome_len(clp)
            divide_bed(self.chromosome_len, self.bin_size)

        self.genome_len = 0
        for chromosome in self.chromosome_len:
            self.genome_len += self.chromosome_len[chromosome][0]["end"] + 1

    def _get_vector(self, bp):
        bed = read_bed(bp)
        divide_bed(bed, self.bin_size)
        merge_bed(bed)
        try:
            vector = bed2vector(bed, self.chromosome_len)
        except:
            print bp
            raise
        return vector


    def _bed2extracted_vector(self, bp):
        return array(extract_vec(self._get_vector(bp), self.dhs_vector), "int8")

    def init_dhs_extraction(self):
        with open(self.dhs_path) as dp:
            dhs = read_bed(dp)
            divide_bed(dhs, self.bin_size)
            merge_bed(dhs)
            self.dhs_vector = bed2vector(dhs, self.chromosome_len)
            coverage = len(nonzero(self.dhs_vector)[0])
            print coverage,"coverage"

        ret = empty((coverage, len(self.beds_path)), dtype="int8")

        current_col = 0
        for bed_path in self.beds_path:
            with open(bed_path) as bp:
                ret[:,current_col] = self._bed2extracted_vector(bp)
                current_col += 1
        self.sparse = csr_matrix(ret,dtype="int16")
        print "init DHS extraction finished"
    def init_factorization(self, rank = 4):
        print "start factorization"
        self.factor = factorization(self.sparse, rank=rank)
        print "factorization finished"

    def classification(self, test_file):
        print "start classification"
        with open(test_file) as tf:
            test_vector = self._bed2extracted_vector(tf)
        test_feature = self.factor["projection"].dot(test_vector)
        nearest_distance = 100
        nearest_index = 0
        distance_list = []
        for i in range(self.factor["feature"].shape[1]):
            tmp_cosine = cosine(self.factor["feature"][:,i], test_feature)
            distance_list.append(tmp_cosine)
            if tmp_cosine < nearest_distance:
                nearest_distance = tmp_cosine
                nearest_index = i
        return {"nearest_distance":nearest_distance,
                "nearest_index": nearest_index,
                "distance_list" : distance_list}

def n_nearest(sparse_obj, names, n=3):
    cnt = len(sparse_obj["distance_list"])
    return nsmallest(n, zip(sparse_obj["distance_list"], [names[i] for i in range(cnt)]))



