import os

folder = "91*2testset_runall2"
for sub_folder in os.listdir(folder):
    pdb_num = 0
    for file in os.listdir(os.path.join(folder, sub_folder)):
        if file.endswith(".pdb"):
            pdb_num += 1
    if pdb_num != 20:
        print(sub_folder)
