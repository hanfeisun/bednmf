import numpy as np
from numpy.ma import nonzero

from scipy.sparse import csc_matrix, csr_matrix


def read_bed(bed_file):
    """
    Convert a bed file to a python object
    """
    ret = {}
    for i in bed_file:
        isp = i.split()
        try:
            chromosome, start, end = isp[0], int(isp[1]), int(isp[2])
        except:
            print isp,"isp"
            raise
        if chromosome not in ret:
            ret[chromosome] = []
        ret[chromosome].append({"start":start, "end":end})
    return ret


def divide_bed(bed_object, bin_size):
    """
    divide a bed object by bin_size
    """
    for chromosome in bed_object:
        for interval in bed_object[chromosome]:
            interval["start"] = interval["start"] // bin_size
            interval["end"] = (interval["end"] - 1) //bin_size + 1


def merge_bed(bed_object):
    """
    merge the sorted bed
    """
    for chromosome in bed_object.keys():
        tmp_chromosome = []
        tmp_interval = {"start": bed_object[chromosome][0]["start"],
                        "end": bed_object[chromosome][0]["end"]}
        for interval in bed_object[chromosome][1:]:
            if interval["start"] > tmp_interval["end"] + 1:
                tmp_chromosome.append(tmp_interval)
                tmp_interval = interval
            else:
                tmp_interval = {"start": min(interval["start"],tmp_interval["start"]),
                                "end": max(interval["end"], tmp_interval["end"])}
        tmp_chromosome.append(tmp_interval)

        bed_object[chromosome] = tmp_chromosome


def read_chromosome_len(len_file):
    """
    Convert a chromosome length file to a bed object
    """
    ret = {}
    for i in len_file:
        isp = i.split()
        chromosome, chromosome_length = isp[0], int(isp[1])
        assert chromosome not in ret  # Don't repeat
        ret[chromosome] = [{"start":0, "end":chromosome_length}]

    return ret


def coverage_bed(bed_obj):
    """
    Stat the coverage base-pairs in a bed object
    """
    coverage = 0
    for chromosome in bed_obj:
        for interval in bed_obj[chromosome]:
            coverage += interval["end"] - interval["start"]
    return coverage

def bed2vector(bed_main, bed_length):
    """
    Input bed object and chromosome length object

    """
    length_of = lambda chromosome: bed_length[chromosome][0]["end"] + 1
    ret = np.zeros(sum(map(length_of,bed_length.keys())), "int8")
    base = 0
    tmp = 0
    for chromosome in sorted(bed_length.keys()):
        for interval in bed_main.get(chromosome, [{"start": 0, "end": 0}]):
            # a dirty trick
            if interval["end"] > length_of(chromosome) - 1:
                interval["end"] = length_of(chromosome) - 1
                if interval["start"] > interval["end"]:
                    print "error"
                    print chromosome, "chromosome isher"
                    raise

            start = base + interval["start"]
            end = base + interval["end"]

            ret[tmp : start] = np.zeros(start - tmp, "int8")
            ret[start : end ] = np.ones(end - start, "int8")
            tmp = end
        base += length_of(chromosome)

    return ret

def extract_vec(target_vector, dhs_vector):
    dhs_index = nonzero(dhs_vector)[0]
    return target_vector[dhs_index]
