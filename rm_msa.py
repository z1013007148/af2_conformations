import os



folder = "91*2testset_runall2"
for sub_folder in os.listdir(folder):
    for file in os.listdir(os.path.join(folder, sub_folder)):
        if file.endswith(".txt"):
            file_path = os.path.join(folder, sub_folder, file)
            os.remove(file_path)
