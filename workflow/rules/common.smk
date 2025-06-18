import glob
import os

fqdir = 'results/fastq'
read1 = [read for read in os.listdir(fqdir) if (read.endswith('R1.fastq.gz'))]
SAMPLES = [sample[:-12] for sample in read1]