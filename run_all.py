import multiprocessing

from af2_conformations.scripts import predict
from af2_conformations.scripts import util
from af2_conformations.scripts import mmseqs2
import shutil
import random
import os

from absl import logging
import time

logging.set_verbosity(logging.DEBUG)
fasta_folder = "../FLEXpdb_database/91test_set/fasta"
a3m_folder = "../FLEXpdb_database/91test_set/msa"
# In this example we introduce alanine mutations into the sequence of MCT1
# MCT1 exclusively adopts the inward-facing conformation when no templates are used
# Alanine mutations were placed at the extracellular gate
pair_dic = {}
with open("/storage/zhouqiangLab/xiayh/zpx/FLEXpdb_database/91test_set/apo_polo_pair.txt") as f:
    for line in f.readlines():
        lines = line.split(",")
        p1 = lines[0]
        p2 = lines[1]
        pair_dic[p1] = p2
        pair_dic[p2] = p1


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
    job_folder, pname, pname_flex_idx_dic = args
    residue_threshold = 2  # 灵活区域两边加2
    random_threshold = 0.6
    fasta_path = os.path.join(fasta_folder, pname + '.fasta')
    if not os.path.exists(fasta_path):
        fasta_path = os.path.join(fasta_folder, pair_dic[pname] + '.fasta')
        if not os.path.exists(fasta_path):
            return
    with open(fasta_path) as f:
        content = f.read().split("\n")
        sequence = "".join(content[1:])
        if len(sequence) < 2:
            return
    a3m_path = os.path.join(a3m_folder, pname, pname + '.a3m')
    if not os.path.exists(a3m_path):
        a3m_path = os.path.join(a3m_folder, pair_dic[pname], pair_dic[pname] + '.a3m')
        if not os.path.exists(a3m_path):
            return
    with open(a3m_path) as f:
        a3m_lines = f.read()
    seq_length = len(a3m_lines.split('\n')[1])

    flex_idxs = pname_flex_idx_dic[pname]

    muts = {x: random.sample("ACDEFGHIKLMNPQRSTVWY-", 1)[0] for x in flex_idxs}

    more_muts = {}
    for i in list(muts.keys()):
        for j in range(i - residue_threshold, i + residue_threshold + 1):
            if j < 0 or j >= seq_length:
                continue
            more_muts[j] = '-'

    for n_model in range(20):
        random_muts = {}
        for idx in random.sample(more_muts.keys(), int(len(more_muts.keys()) * random_threshold)):  # 随机取threshold%的残基
            random_muts[idx] = more_muts[idx]
        # Define the mutations and introduce into the sequence and MSA
        mutated_msa = util.mutate_msa(a3m_lines, random_muts)
        # mutated_seq = util.mutate_msa(sequence, muts)
        mutated_seq = sequence  # 序列不变

        with open(f"{job_folder}/{pname}_mutated_seq_{n_model}.txt", 'w') as f:
            f.write(str(mutated_seq))
        with open(f"{job_folder}/{pname}_mutated_msa_{n_model}.txt", 'w') as f:
            f.write(str(mutated_msa))

        # Specify the name of the output PDB
        outname = f"{pname}_model_{n_model}.pdb"
        predict.predict_structure_no_templates(mutated_seq, outname, mutated_msa)
        shutil.copyfile(outname, os.path.join(job_folder, outname))
        os.remove(outname)


if __name__ == '__main__':
    output_folder = "91*2testset_runall"
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    flex_txt = "/storage/zhouqiangLab/xiayh/zpx/af2_conformations_master/flex_data2.txt"
    pname_flex_idx_dic = {}
    with open(flex_txt) as f:
        for line in f.readlines():
            if len(line) < 2:
                continue
            lines = line.split(":")
            pname = lines[0]
            flexs = lines[1].strip("\n").strip("[").strip("]")
            flex_idxs = flexs.split(", ")
            flex_idxs = [int(i) for i in flex_idxs]
            pname_flex_idx_dic[pname] = flex_idxs

    pool = multiprocessing.Pool(4)
    manager = multiprocessing.Manager()
    lock = manager.RLock()
    gpu_dic = manager.dict()
    for i in range(0, 4):
        gpu_dic[i] = True
    # for pname in pname_flex_idx_dic:
    for pname in ["1WCW_A"]:  # 第二次跑，没跑完的几个蛋白
        job_folder = os.path.join(output_folder, pname)
        if not os.path.exists(job_folder):
            os.mkdir(job_folder)
        pool.apply_async(gpu_queue, ((job_folder, pname, pname_flex_idx_dic), gpu_dic, lock,))
    pool.close()
    pool.join()
