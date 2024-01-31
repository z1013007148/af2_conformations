import multiprocessing

from af2_conformations.scripts import predict
from af2_conformations.scripts import util
from af2_conformations.scripts import mmseqs2
import shutil
import random
import os
import uuid
from absl import logging

os.environ["CUDA_VISIBLE_DEVICES"] = '3'

logging.set_verbosity(logging.DEBUG)
fasta_folder = "../FLEXpdb_database/91test_set/fasta"
a3m_folder = "../FLEXpdb_database/91test_set/msa"
# In this example we introduce alanine mutations into the sequence of MCT1
# MCT1 exclusively adopts the inward-facing conformation when no templates are used
# Alanine mutations were placed at the extracellular gate
pname = '1F3Y-17_A'
sequence = (
    "GPLGSMDSPPEGYRRNVGICLMNNDKKIFAASRLDIPDAWQMPQGGIDEGEDPRNAAIRELREETGVTSAEVIAEVPYWLTYDFPPKVREKLNIQWGSDWKGQAQKWFLFKFTGQDQEINLLGDGSEKPEFGEWSWVTPEQLIDLTVEFKKPVYKEVLSVFAPHL")

# The MMSeqs2Runner object submits the amino acid sequence to
# the MMSeqs2 server, generates a directory, and populates it with
# data retrieved from the server.
# mmseqs2_runner = mmseqs2.MMSeqs2Runner(jobname, sequence)
#
# # Fetch sequences and download data
# a3m_lines, _ = mmseqs2_runner.run_job()

a3m_path = os.path.join(a3m_folder, pname, pname + '.a3m')
# with open(a3m_path) as f:
#     a3m_lines = f.read()
#     a3m_list = a3m_lines.split('>')[1:]
#     if len(a3m_list)/2>1024:
#         a3m_random_list = random.sample(a3m_list, 1024)
#     a3m_random_list = '>'.join(a3m_random_list)

import multiprocessing
import time
# lock = multiprocessing.Lock()

def job(v, num,lock):
    with lock:
        for _ in range(5):
            print(v+num, end=",")

def init(l):
    global lock
    lock = l


def multicore():
    lock = multiprocessing.Manager().RLock()
    pool = multiprocessing.Pool(2)

    v =1
    pool.apply_async(job, (v, 1, lock,))
    pool.apply_async(job, (v, 3,lock,))

    pool.close()
    pool.join()

if __name__ == '__main__':

    multicore()
# with open(MSA_path, 'r') as f1:
#     line_1 = f1.readlines()
#     N = len(line_1)
#     if N > 128:
#         randomList = random.sample(range(1, N), 127)
#         for i, line_64 in enumerate(line_1):
#             if i == 0:
#                 with open(msa_128_path, 'a') as f2:
#                     f2.write('>')
#                     f2.write('\n')
#                     f2.write(line_64)
#             elif i in randomList:
#                 with open(msa_128_path, 'a') as f2:
#                     f2.write('>')
#                     f2.write('\n')
#                     f2.write(line_64)
#     else:
#         with open(msa_128_path, 'a') as f2:
#             for line_64 in line_1:
#                 f2.write('>')
#                 f2.write('\n')
#                 f2.write(line_64)

# muts = {x: random.sample("ACDEFGHIKLMNPQRSTVWY-", 1) for x in
#             [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100,
#              101, 121, 122, 123, 124, 125, 126, 127, 128, 129, 164]}



