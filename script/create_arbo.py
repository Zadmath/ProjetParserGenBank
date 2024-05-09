from ftplib import FTP
import os
import shutil
import pandas as pd
import pickle
import random
import string
import time
import requests
import sys
import stat
import genericpath
from multiprocessing import Pool, Process, Manager, Queue,cpu_count, Semaphore, Value, Event

from genericpath import *
from Bio import Entrez
from Bio import SeqIO

# parse overview.txt
organism_names = []
organism_paths = []

with open('../GENOME_REPORTS/overview.txt') as f:
    print("traitement overview.txt...")
    first_row = True
    count_rows = 1
    for row in f:
        count_rows += 1
        if first_row:
            first_row=False
            continue
        parsed_row = row.split('\t')

        try :
            organism = parsed_row[0].replace(' ','_').replace('/','_')
            kingdom = parsed_row[1].replace(' ','_').replace('/','_')
            group = parsed_row[2].replace(' ','_').replace('/','_')
            subgroup = parsed_row[3].replace(' ','_').replace('/','_')
            path = '../Results/' + kingdom +'/' + group +'/' + subgroup +'/' + organism
            organism_names.append(parsed_row[0])
            organism_paths.append('../Results/' + kingdom +'/' + group +'/' + subgroup +'/')
        except IndexError :
            print(f"IndexError sur {parsed_row}")
            pass
    print("terminé traitement overview.txt")

# parse fichier ids
ids_files = os.listdir('../GENOME_REPORTS/IDS/')
organism_names_ids = []
organism_paths_ids = []
organism_NC_ids = []
i = 0

for ids in ids_files:
    print(f"Traitement ids : {ids}")
    i += 1

    with open('../GENOME_REPORTS/IDS/' + ids) as f:
        n_line = sum(1 for _ in f)

    with open('../GENOME_REPORTS/IDS/' + ids) as f:
        for i,row in enumerate(f):
            # update l'affichage de la progression

            parsed_row = row.replace('\n', '').split('\t')
            if (parsed_row[1][0:2] != 'NC'):
                continue

            try:
                index = organism_names.index(parsed_row[5])

            except ValueError:
                continue
            try:
                if parsed_row[1] not in organism_NC_ids[organism_names_ids.index(organism_names[index])]:
                    organism_NC_ids[organism_names_ids.index(organism_names[index])].append(parsed_row[1])
            except ValueError:

                organism_names_ids.append(organism_names[index])
                organism_paths_ids.append(organism_paths[index])
                organism_NC_ids.append([parsed_row[1]])
                name = organism_names[index].replace(" ", "_")
                name = name.replace("[", "_")
                name = name.replace("]", "_")
                name = name.replace(":", "_")
                name = name.replace("/", "_")
                path = organism_paths[index] + name + "/"
                if not os.path.exists(path):
                    os.makedirs(path)


#dataframe
organism_df = pd.DataFrame({
            "name":organism_names_ids,
            "path":organism_paths_ids,
            "NC":organism_NC_ids})

# crée le fichier pickle sauvegardant le dataframe
if not os.path.exists("../pickle"):
    os.makedirs("../pickle")

