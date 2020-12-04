import numpy
import pandas
import os
import itertools
import argparse
from tqdm import tqdm
numpy.random.seed(0)


parser = argparse.ArgumentParser()
parser.add_argument('-a', '--assay', type=str, default=None,
        help="comma delimited list of assays for each sample; ex H3K27ac, H3K4me2...")
parser.add_argument('-f', '--file', type=str, default=None,
        help="tab delimited file, cell type, assay, name")
parser.add_argument('-c', '--chromsize', type=str, default="/home/CAM/ewentworth/cotney/genome/hg19/hg19.chrom.sizes.avocado", help='chromosome size file')


args = parser.parse_args()

assays, names = [], []
with open(args.file, "r") as infile:
        for line in infile:
                        line = line.strip("\r\n").split()
                        assay = line[0]
                        name = line[1]
                        assays.append(assay)
                        names.append(name)

chroms= []
with open(args.chromsize, "r") as infile:
        for line in infile:
                        chrom = line.strip("\r\n").split()
                        chroms.append(chrom)


def bedgraph_to_dense(filename, verbose=True):
        bedgraph = pandas.read_csv(filename, sep="\t", header=None)
        n = bedgraph[2].values[-1]
        k = bedgraph.shape[0]
        data = numpy.zeros(n)
        d = not verbose
        for i, (_, start, end, v) in tqdm(bedgraph.iterrows(), total=k, disable=d):
                data[start:end] = v
        return data



def decimate_vector(x, k=25, func=numpy.mean):
        m = x.shape[0] // k
        y = numpy.zeros(m)
        for i in range(m):
                y[i] = func(x[i*k:(i+1)*k])
        return y



for chrom, assay in itertools.product(chroms, assays):
        data = bedgraph_to_dense("{}.{}.bedGraph".format(assay, chrom[0]))
        print("{}.{}.bedGraph".format(assay, chrom[0]))
        data = decimate_vector(data)
        data_ = numpy.zeros(int("{}".format(chrom[1])))
        data_[:len(data)] = data
        data = data_
        data = numpy.arcsinh(data).astype('float32')
        numpy.savez_compressed("{}.{}.arcsinh.npz".format(assay, chrom[0]), data)

















