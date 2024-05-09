from ftplib import FTP
import os
import shutil
import pandas as pd
import pickle
import random
import string
import time
import requests
from Bio import Entrez
from Bio import SeqIO
if os.path.exists('../GENOME_REPORTS/overview.txt') and os.path.isfile("../pickle/organism_df"):
    t = os.path.getmtime("../GENOME_REPORTS/overview.txt")


print("Suppression de l'arborescence actuelle...")
try: shutil.rmtree('../GENOME_REPORTS')
except: pass
print("terminé.")
os.mkdir("../GENOME_REPORTS")
os.chdir('../GENOME_REPORTS')
os.mkdir('IDS')

with FTP('ftp.ncbi.nlm.nih.gov') as ftp:
    ftp.login()  # connecter au FTP

    ftp.cwd('genomes/GENOME_REPORTS')

    print(" - overview.txt")
    ftp.retrbinary('RETR overview.txt', open("overview.txt",'wb').write)

    ftp.cwd('IDS')  # changer de répertoire courant
    n = len(ftp.nlst())
    for i,filename in enumerate(ftp.nlst()):
        print(f" - {filename}")
        ftp.retrbinary('RETR '+ filename, open("IDS/" + filename, 'wb').write)

