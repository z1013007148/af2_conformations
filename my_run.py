from af2_conformations.scripts import predict
from af2_conformations.scripts import util
from af2_conformations.scripts import mmseqs2

import random
import os

from absl import logging
logging.set_verbosity(logging.DEBUG)

# In this example we introduce alanine mutations into the sequence of MCT1
# MCT1 exclusively adopts the inward-facing conformation when no templates are used
# Alanine mutations were placed at the extracellular gate
fasta_folder = "/home/data/user/guangxingh/zpx/FLEXpdb_database/91test_set/fasta"
a3m_folder = "/home/data/user/guangxingh/zpx/FLEXpdb_database/91test_set/msa"
for fasta in os.listdir(fasta_folder):
    pname = fasta.replace('.fasta', '')
    fasta_path = os.path.join(fasta_folder, fasta)
    with open(fasta_path) as f:
        content = f.read()
        sequence = ''.join(content.split('\n')[1:])

# jobname = 'MCT1'
# sequence = ("MPPAVGGPVGYTPPDGGWGWAVVIGAFISIGFSYAFPKSITVFFKEIEGIFHATTSEVSWISS"
#             "IMLAVMYGGGPISSILVNKYGSRIVMIVGGCLSGCGLIAASFCNTVQQLYVCIGVIGGLGLAF"
#             "NLNPALTMIGKYFYKRRPLANGLAMAGSPVFLCTLAPLNQVFFGIFGWRGSFLILGGLLLNCC"
#             "VAGALMRPIGPKPTKAGKDKSKASLEKAGKSGVKKDLHDANTDLIGRHPKQEKRSVFQTINQF"
#             "LDLTLFTHRGFLLYLSGNVIMFFGLFAPLVFLSSYGKSQHYSSEKSAFLLSILAFVDMVARPS"
#             "MGLVANTKPIRPRIQYFFAASVVANGVCHMLAPLSTTYVGFCVYAGFFGFAFGWLSSVLFETL"
#             "MDLVGPQRFSSAVGLVTIVECCPVLLGPPLLGRLNDMYGDYKYTYWACGVVLIISGIYLFIGM"
#             "GINYRLLAKEQKANEQKKESKEEETSIDVAGKPNEVTKAAESPDQKDTDGGPKEEESPV" )

# The MMSeqs2Runner object submits the amino acid sequence to
# the MMSeqs2 server, generates a directory, and populates it with
# data retrieved from the server.
#     mmseqs2_runner = mmseqs2.MMSeqs2Runner( jobname, sequence )

    # Fetch sequences and download data
    # a3m_lines, _ = mmseqs2_runner.run_job()
    a3m_path = os.path.join(a3m_folder, pname, pname+'.a3m')
    with open(a3m_path) as f:
        a3m_lines = f.read()
    # Define the mutations and introduce into the sequence and MSA
    muts = { x: "A" for x in [ 41,42,45,46,56,59,60,63,281,282,285,286,403,407 ] }

    mutated_msa = util.mutate_msa( a3m_lines, muts )
    mutated_seq = util.mutate_msa( sequence, muts )

    for n_model in range( 1 ):

      # Specify the name of the output PDB
      outname = f"model_{ n_model }.pdb"

      # predict.predict_structure_no_templates( mutated_seq, outname, mutated_msa )
      with open("mutated_seq.txt",'w') as f:
          f.write(str(mutated_seq))
      with open("mutated_msa.txt",'w') as f:
          f.write(str(mutated_msa))