import numpy
import pandas
import os
os.environ['KERAS_BACKEND'] = 'theano'
import keras
import itertools
import theano
import argparse
from tqdm import tqdm
from avocado import Avocado

chroms = list(range(1, 23)) + ['X']

parser = argparse.ArgumentParser()
parser.add_argument('-a', '--assay', type=str, default=None,
	help="comma delimited list of assays for each sample; ex H3K27ac, H3K4me2...")
parser.add_argument('-f', '--file', type=str, default=None,
	help="tab delimited file, cell type, assay, name")
parser.add_argument('-c', '--chromsize', type=str, default="/home/CAM/ewentworth/cotney/genome/hg19/hg19.chrom.sizes.avocado", 
	help='chromosome size file')
		
		

args = parser.parse_args()

chroms= []
with open(args.chromsize, "r") as infile:
        for line in infile:
                        chrom = line.strip("\r\n").split()
                        chroms.append(chrom)



celltypes, assays = [], []
with open(args.file, "r") as infile:
	for line in infile:
			line = line.strip("\r\n").split()
			celltype = line[0]
			assay = line[1]
			assays.append(assay)
			celltypes.append(celltype)


data = {}
for celltype, assay, chrom in itertools.product(celltypes, assays, chroms):
	filename = '{}_{}.{}.arcsinh.npz'.format(celltype, assay, chrom[0])
	print(filename)
	data[(celltype, assay)] = numpy.load(filename)['arr_0']
	modelname = 'avocado-{}'.format(chrom[0])
	print(modelname)
	Model = Avocado.load(modelname)
	Model.fit_celltypes(data, n_epochs=50)
	Model.save("model.{}".format(chrom[0]))





