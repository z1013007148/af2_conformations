import os
import shutil

p_dic = {}
with open("/storage/zhouqiangLab/xiayh/zpx/FLEXpdb_database/91test_set/apo_polo_pair.txt") as f:
    for line in f.readlines():
        lines = line.split(',')
        p1 = lines[0]
        p2 = lines[1]
        p_dic[p1] = p2
        p_dic[p2] = p1
folder = "91*2testset_runall"
for pname in os.listdir(folder):
    dst_path1 = os.path.join(folder, pname, pname + '.pdb')
    dst_path2 = os.path.join(folder, pname, p_dic[pname] + '.pdb')
    source_path1 = os.path.join("/storage/zhouqiangLab/xiayh/zpx/FLEXpdb_database/91test_set/pdb", pname + '.pdb')
    source_path2 = os.path.join("/storage/zhouqiangLab/xiayh/zpx/FLEXpdb_database/91test_set/pdb",
                                p_dic[pname] + '.pdb')
    if not os.path.exists(source_path1) or not os.path.exists(source_path2):
        print(pname)
        continue
    shutil.copyfile(source_path1, dst_path1)
    shutil.copyfile(source_path2, dst_path2)
