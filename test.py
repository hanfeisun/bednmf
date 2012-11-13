import unittest
from scipy.sparse import csr_matrix

from scipy import array
from numpy import dot, zeros, testing, shape, alltrue, ravel
from scipy.spatial.distance import cosine
import time
from bed import read_bed, divide_bed, merge_bed, bed2vector, read_chromosome_len
from classification import BedNMF
from nmf import factorization
def read_vector(vector_file):

    """
    read estimated vector
    """
    ret = []
    dic = {}
    for vector in vector_file:
        vector = vector.split(" ")
        dic[vector[0]] = map(int, vector[1:])
    for chromosome in sorted(dic.keys()):
        ret.extend(dic[chromosome])
    return ret



class NmfTestCase(unittest.TestCase):
    def setUp(self):
        self.V = csr_matrix((array([1,1,1,1,1,1]), array([0,2,2,0,1,2]), array([0,2,3,6])), shape=(3,3))
    def tearDown(self):
        self.V = None
    def testNmfRun(self):
        self.factor = factorization(self.V, 4)
        self.estimate = dot(self.factor["basis"], self.factor["coef"])
        print self.estimate
        self.error = self.estimate - self.V


class BedTestCase(unittest.TestCase):
    def setUp(self):
        with open("./utest/test1.bed") as tb:
            self.bed_obj = read_bed(tb)
        self.maxDiff = None
    def tearDown(self):
        self.bed_obj = None
    def testDivide(self):
        divide_bed(self.bed_obj, 100)
        with open("./utest/test1_divided.bed") as dd:
            self.bed_div = read_bed(dd)
        self.assertEqual(self.bed_div, self.bed_obj)
    def testMerge(self):
        merge_bed(self.bed_obj)
        with open("./utest/test1_merged.bed") as md:
            self.bed_mer = read_bed(md)
        self.assertEqual(self.bed_obj, self.bed_mer)
    def testMerge2(self):
        divide_bed(self.bed_obj, 100)
        merge_bed(self.bed_obj)
        with open("./utest/test1_divided_merged.bed") as dmd:
            self.bed_dmer = read_bed(dmd)
        self.assertEqual(self.bed_obj, self.bed_dmer)

class VectorTestCase(unittest.TestCase):
    def setUp(self):
        with open("./utest/test1.len") as tl:
            self.len_obj = read_chromosome_len(tl)
            divide_bed(self.len_obj, 100)
    def tearDown(self):
        self.bed_obj = None
        self.len_obj = None


    def gen_testVector(self, bed, vec):
        with open(bed) as md:
            bed_obj = read_bed(md)
        with open(vec) as dmv:
            vector_dmv = array(read_vector(dmv),"int8")
        self.vec_obj = bed2vector(bed_obj, self.len_obj)
        print "\nestimate\t", vector_dmv
        print "result\t\t", self.vec_obj
        self.assertTrue(alltrue(self.vec_obj== vector_dmv))

    def testVector(self):
        self.gen_testVector("./utest/test1_divided_merged.bed", "./utest/test1_divided_merged.vec")

    def testVector2(self):
        """
        This dataset has a error peak
        """
        self.gen_testVector("./utest/test1_different.bed", "./utest/test1_different.vec")


class SparseTestCase(unittest.TestCase):
    def setUp(self):
        self.file_bed =  ["./utest/test1_divided_merged.bed",
                          "./utest/test1_different.bed"]
        self.file_len = "./utest/test1_divided.len"
        self.file_test = "./utest/test1_similar.bed"
        self.file_dhs = "./utest/test1_dhs.bed"

        self.estimate_extracted = ["./utest/test1_divided_merged_extracted.vec",
                                   "./utest/test1_different_extracted.vec"]
        self.bsm = BedNMF(self.file_bed,self.file_len,self.file_dhs, 1)


    def testDhsExtraction(self):
        self.bsm.init_dhs_extraction()

        cnt = 0
        for efile in self.estimate_extracted:
            with open(efile) as ef:
                vector = read_vector(ef)
            result = ravel(self.bsm.sparse[:,cnt].todense())
            print "estimate extracted",vector
            print "result extracted", result
            self.assertTrue(alltrue(vector == result ))
            cnt += 1
    def testNmfClassification(self):
        self.bsm.init_dhs_extraction()
        self.bsm.init_factorization()
        print self.bsm.classification(self.file_bed[0])
        print self.bsm.classification(self.file_bed[1])
        print self.bsm.classification(self.file_test)

if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(BedTestCase('testDivide'))
    suite.addTest(BedTestCase('testMerge'))
    suite.addTest(BedTestCase('testMerge2'))
    suite.addTest(VectorTestCase('testVector'))
    suite.addTest(VectorTestCase('testVector2'))
    suite.addTest(SparseTestCase('testDhsExtraction'))
    suite.addTest(SparseTestCase('testNmfClassification'))
    runner = unittest.TextTestRunner()
    runner.run(suite)


