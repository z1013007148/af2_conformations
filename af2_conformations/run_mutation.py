from scripts import predict
from scripts import util
from scripts import mmseqs2

import random
import os

from absl import logging
logging.set_verbosity(logging.DEBUG)

# In this example we introduce alanine mutations into the sequence of MCT1
# MCT1 exclusively adopts the inward-facing conformation when no templates are used
# Alanine mutations were placed at the extracellular gate
jobname = 'MCT1'
sequence = ("MPPAVGGPVGYTPPDGGWGWAVVIGAFISIGFSYAFPKSITVFFKEIEGIFHATTSEVSWISS"
            "IMLAVMYGGGPISSILVNKYGSRIVMIVGGCLSGCGLIAASFCNTVQQLYVCIGVIGGLGLAF"
            "NLNPALTMIGKYFYKRRPLANGLAMAGSPVFLCTLAPLNQVFFGIFGWRGSFLILGGLLLNCC"
            "VAGALMRPIGPKPTKAGKDKSKASLEKAGKSGVKKDLHDANTDLIGRHPKQEKRSVFQTINQF"
            "LDLTLFTHRGFLLYLSGNVIMFFGLFAPLVFLSSYGKSQHYSSEKSAFLLSILAFVDMVARPS"
            "MGLVANTKPIRPRIQYFFAASVVANGVCHMLAPLSTTYVGFCVYAGFFGFAFGWLSSVLFETL"
            "MDLVGPQRFSSAVGLVTIVECCPVLLGPPLLGRLNDMYGDYKYTYWACGVVLIISGIYLFIGM"
            "GINYRLLAKEQKANEQKKESKEEETSIDVAGKPNEVTKAAESPDQKDTDGGPKEEESPV" )

# The MMSeqs2Runner object submits the amino acid sequence to
# the MMSeqs2 server, generates a directory, and populates it with
# data retrieved from the server.
mmseqs2_runner = mmseqs2.MMSeqs2Runner( jobname, sequence )

# Fetch sequences and download data
a3m_lines, _ = mmseqs2_runner.run_job()

# Define the mutations and introduce into the sequence and MSA
muts = { x: "A" for x in [ 41,42,45,46,56,59,60,63,281,282,285,286,403,407 ] }

mutated_msa = util.mutate_msa( a3m_lines, muts )
mutated_seq = util.mutate_msa( sequence, muts )

for n_model in range( 5 ):

  # Specify the name of the output PDB
  outname = f"model_{ n_model }.pdb"

  predict.predict_structure_no_templates( mutated_seq, outname, mutated_msa )