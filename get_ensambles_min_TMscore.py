import multiprocessing
import os
import subprocess




def cal_TMscore(pdb1, pdb2, dic, lddt=False):
    if not lddt:
        out = subprocess.check_output(
            ['/storage/zhouqiangLab/xiayh/zpx/af2_conformations_master/TMscore', '-seq', pdb1, pdb2],
            # 改为执行当前目录下的TMscore可执行文件
            stderr=open('/dev/null', 'w'), cwd='./')
        start = out.find(b'RMSD')
        end = out.find(b'rotation')
        out = out[start:end]

        rmsd, _, tm, _, gdt_ts, gdt_ha, _, _ = out.split(b'\n')

        rmsd = float(rmsd.split(b'=')[-1])
        tm = float(tm.split(b'=')[1].split()[0])
        gdt_ts = float(gdt_ts.split(b'=')[1].split()[0])
        gdt_ha = float(gdt_ha.split(b'=')[1].split()[0])

        dic_key = os.path.split(pdb1)[-1] + '===' + os.path.split(pdb2)[-1]
        dic[dic_key] = tm
    else:
        out = subprocess.check_output(['/storage/zhouqiangLab/xiayh/zpx/af2_conformations_master/lddt', '-c', pdb1, pdb2],  # reference comes last  # 改为执行当前目录下的lddt可执行文件
            stderr=open('/dev/null', 'w'), cwd='./')
        for line in out.split(b'\n'):
            if b'Global LDDT score' in line:
                lddt = float(line.split(b':')[-1].strip())
                break
        dic_key = os.path.split(pdb1)[-1] + '===' + os.path.split(pdb2)[-1]
        dic[dic_key] = lddt


# ensemble_folder = "/storage/zhouqiangLab/xiayh/zpx/af2_conformations_master/1F3Y-17_A_randommsaaa_0.7"
if __name__ == '__main__':
    for ensemble_folder in [f"/storage/zhouqiangLab/xiayh/zpx/af2_conformations_master/1F3Y-17_A_randommsaaa_{i}"for i in [0.1,0.5,0.7,0.9]]:
        pool = multiprocessing.Pool(50)
        manager = multiprocessing.Manager()
        dic = manager.dict()
        pdb_list = []
        for file in os.listdir(ensemble_folder):
            if not file.endswith('.pdb'):
                continue
            file_path = os.path.join(ensemble_folder, file)
            pdb_list.append(file_path)

        for i in range(len(pdb_list)):
            for j in range(i + 1, len(pdb_list)):
                pool.apply_async(cal_TMscore, (pdb_list[i], pdb_list[j], dic, False,))
                # cal_TMscore(pdb_list[i], pdb_list[j], dic, True)
        pool.close()
        pool.join()
        # print(len(dic))
        sorted_list = sorted(dic.items(), key=lambda e: e[1])
        print('TMscore最小的结构对是', sorted_list[0])

