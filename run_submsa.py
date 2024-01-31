import multiprocessing

from af2_conformations.scripts import predict
from af2_conformations.scripts import util
from af2_conformations.scripts import mmseqs2
import shutil
import random
import os
import time
from absl import logging


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
with open(a3m_path) as f:
    a3m_lines = f.read()


# muts = {x: random.sample("ACDEFGHIKLMNPQRSTVWY-", 1) for x in
#             [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100,
#              101, 121, 122, 123, 124, 125, 126, 127, 128, 129, 164]}

def gpu_queue(args, gpu_dic, lock):
    gpu_idx = 0
    use_gpu_flag = False
    while not use_gpu_flag:
        with lock:
            for i in gpu_dic.keys():
                if gpu_dic[i]:
                    gpu_dic[i] = False
                    gpu_idx = i
                    use_gpu_flag = True
                    break

        if not use_gpu_flag:
            print('没抢到gpu')
            time.sleep(5000)

    print('抢到gpu', gpu_idx)
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_idx)
    try:
        deal(args)
    except Exception as e:
        print(e)
    gpu_dic[gpu_idx] = True


def deal(args):
    jobname, threshold = args
    for n_model in range(20):
        a3m_list = a3m_lines.split('>')[1:]
        if len(a3m_list) > threshold:
            a3m_random_list = random.sample(a3m_list, threshold)
            random_a3m_lines = '>'.join(['']+a3m_random_list)
        else:
            random_a3m_lines = a3m_lines
        muts = {x: '-' for x in
                [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99,
                 100, 101, 121, 122, 123, 124, 125, 126, 127, 128, 129, 164]}
        random_muts = {}
        for idx in random.sample(muts.keys(), int(len(muts.keys()))):  # 随机取threshold%的残基
            random_muts[idx] = muts[idx]
        # Define the mutations and introduce into the sequence and MSA
        mutated_msa = util.mutate_msa(random_a3m_lines, random_muts)
        # mutated_seq = util.mutate_msa(sequence, muts)
        mutated_seq = sequence  # 序列不变

        with open(f"{jobname}/{jobname}_mutated_seq_{n_model}.txt", 'w') as f:
            f.write(str(mutated_seq))
        with open(f"{jobname}/{jobname}_mutated_msa_{n_model}.txt", 'w') as f:
            f.write(str(mutated_msa))

        # Specify the name of the output PDB
        outname = f"{jobname}_model_{n_model}.pdb"
        predict.predict_structure_no_templates(mutated_seq, outname, mutated_msa)
        shutil.copyfile(outname, os.path.join(jobname, outname))
        os.remove(outname)


if __name__ == '__main__':
    pool = multiprocessing.Pool(4)
    manager = multiprocessing.Manager()
    lock = manager.RLock()
    gpu_dic = manager.dict()
    for i in range(0, 4):
        gpu_dic[i] = True
    for threshold in [128, 256, 512, 1024, 2048]:
        jobname = f"{pname}_submsa_{threshold}"
        if not os.path.exists(jobname):
            os.mkdir(jobname)
        pool.apply_async(gpu_queue, ((jobname, threshold), gpu_dic, lock))
    pool.close()
    pool.join()
