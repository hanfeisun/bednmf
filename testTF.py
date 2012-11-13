from pprint import pprint
from bed import read_chromosome_len, divide_bed, read_bed, merge_bed
from classification import BedNMF, n_nearest
import cPickle
__author__ = 'hanfeisun'

import unittest
import glob



class TfTestCase(unittest.TestCase):
    def setUp(self):
        self.train_path = glob.glob("../TF_train/*.bed")
        self.len_path = "./static_data/hg19_len"
        self.dhs_path = "./static_data/DHS_hg19.bed"
        self.test_path = "./CTCF_test/5011_peaks.bed"
        self.test2_path = "./CTCF_train/1252_peaks.bed"
        self.dump_path = "./persist/all_TF_train.p"
        self.bin_size = 50000
        self.InitNMF()

    def InitNMF(self):
        try:
            self.bsm = cPickle.load(open(self.dump_path,"r"))
            print "Load from database"
        except:
            print "Init database"
            self.bsm = BedNMF(self.train_path, self.len_path, self.dhs_path, self.bin_size)
            self.bsm.init_dhs_extraction()
            self.bsm.init_factorization(rank=20)
            cPickle.dump(self.bsm, open(self.dump_path, 'w'))
    def gen_testClassification(self, test_path):
        result = self.bsm.classification(test_path)
        pprint(n_nearest(result, self.train_path))
        return result
    def testClassification1(self):
        self.gen_testClassification(self.test_path)
        self.gen_testClassification(self.test2_path)
    #    def testClassification2(self):
#










if __name__ == "__main__":
    unittest.main()



