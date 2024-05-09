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

df = pd.read_pickle("../pickle/organism_df")
df["nb_NC"] = df["NC"].apply(len)

homo = df[df['name'] == "Homo Sapiens"]


NC_list = homo["NC"].iloc[0]
NC_list = list(dict.fromkeys(NC_list))


def f(NC):

    Entrez.email = ''.join(random.choice(string.ascii_lowercase) for i in range(20))+'@'+''.join(random.choice(string.ascii_lowercase) for i in range(20))+ '.com'
    try :
        print(f"  request for {NC} text")
        handle_fasta = Entrez.efetch(db="nucleotide", id=NC, rettype="fasta", retmode="text")
    except Exception as e:
        print(f"1 for {NC} text",str(e))
        if str(e) == "HTTP Error 429: Too Many Requests":
            iter = 0
            iter_max = 10
            while True:
                time.sleep(10) # on se calme un peu
                try :
                    print(f"  request for {NC} text (rerun 1)")
                    handle_fasta = Entrez.efetch(db="nucleotide", id=NC, rettype="fasta", retmode="text")
                    break
                except Exception as e:
                    iter += 1
                    print(f"1.{iter} for {NC} text",str(e))
                    if str(e) == "HTTP Error 429: Too Many Requests":
                        continue
                    elif iter > iter_max:
                        return
        else:
            return
    try:
        # record_fasta = timeout(timeout=TIMEOUT_MAX)(SeqIO.read)(handle_fasta, "fasta")
        record_fasta = SeqIO.read(handle_fasta, "fasta")
    except Exception as e:
        print(f"2 for {NC} text",e)
        return
    handle_fasta.close()

def f2(number_region_found, number_region_already_found, semaphore, path, NC, name, selected_region):
    """
    renvois un couple (résultat, code)
    où résultat est un string décrivant le succès ou non de la requete
    et code est la couleur avec laquelle le résultat devrait être affiché
    """
    # if test_NC_downloaded(path, NC, name, selected_region):
    #     print(f"  already downloaded : {path}/{name}/{NC}")
    #     nb_region_already_found = number_region_already_found.get()
    #     nb_region_already_found += 1
    #     number_region_already_found.put(nb_region_already_found)
    #     return f"Déjà traité","green"

    name = name.replace(" ", "_")
    name = name.replace("[", "_")
    name = name.replace("]", "_")
    name = name.replace(":", "_")
 
    Entrez.email = ''.join(random.choice(string.ascii_lowercase) for i in range(20))+'@'+''.join(random.choice(string.ascii_lowercase) for i in range(20))+ '.com'
    # try :
    #     print(f"  request for {NC} text")
    #     handle_fasta = Entrez.efetch(db="nucleotide", id=NC, rettype="fasta", retmode="text")
    # except Exception as e:
    #     print(f"1 for {NC} text",str(e))
    #     if str(e) == "HTTP Error 429: Too Many Requests":
    #         iter = 0
    #         iter_max = 10
    #         while True:
    #             time.sleep(10) # on se calme un peu
    #             try :
    #                 print(f"  request for {NC} text (rerun 1)")
    #                 handle_fasta = Entrez.efetch(db="nucleotide", id=NC, rettype="fasta", retmode="text")
    #                 break
    #             except Exception as e:
    #                 iter += 1
    #                 print(f"1.{iter} for {NC} text",str(e))
    #                 if str(e) == "HTTP Error 429: Too Many Requests":
    #                     continue
    #                 elif iter > iter_max:
    #                     return
    #     else:
    #         return
    # try:
    #     # record_fasta = timeout(timeout=TIMEOUT_MAX)(SeqIO.read)(handle_fasta, "fasta")
    #     record_fasta = SeqIO.read(handle_fasta, "fasta")
    # except Exception as e:
    #     print(f"2 for {NC} text",e)
    #     return
    # handle_fasta.close()

    iter = 0
    while True:
        try:
            iter += 1
            print(f"  request for {NC} xml attempt {iter}")
            # handle = Entrez.efetch(db="nuccore", id=NC, rettype="gbwithparts", retmode="text")
            handle = Entrez.efetch(db="nuccore", id=NC, rettype="gbwithparts", retmode="xml")
            # handle = Entrez.efetch(db="nuccore", id=NC, retmode="xml")
            break
        except Exception as e:
            print(f"  error for request {NC} xml attempt {iter}: {str(e)}")
            if '429' in str(e):
                print('  -> waiting 4 sec and retrying')
                time.sleep(4)
            else:
                print("  unhandeld error : panic exit")
                return f"Erreur non gérée lors de la requête à l'API","red"
    try:
        # record = SeqIO.read(handle, 'gb')
        record = Entrez.read(handle)
        handle.close()
    except Exception as e:
        print(f"  error for request {NC} (reading) : {str(e)}")
        return "Erreur de lecture du fichier brut (gbwithparts)","red"
    iter=0
    while True:
        try:
            iter += 1
            print(f"  request for {NC} fasta attempt {iter}")
            handle_fasta = Entrez.efetch(db="nuccore", id=NC, rettype="fasta", retmode="text")
            break
        except Exception as e:
            print(f"  error for request {NC} fasta attempt {iter}: {str(e)}")
            if '429' in str(e):
                print('  -> waiting 4 sec and retrying')
                time.sleep(4)
            else:
                print("  unhandeld error : panic exit")
                return "Erreur non gérée lors de la requête à l'API","red"
    try:
        # record = SeqIO.read(handle, 'gb')
        record_fasta = SeqIO.read(handle_fasta,"fasta")
        handle_fasta.close()
    except Exception as e:
        print(f"  error for request {NC} fasta (reading) : {str(e)}")
        return "Erreur de lecture du fichier brut (fasta)","red"



    # print(f"  request for {NC} xml")
    # handle_text = Entrez.efetch(db="nucleotide", id=NC, retmode="xml")
    # try:
    #     record = Entrez.read(handle_text)
    # except Exception as e:
    #     print(f"3 for {NC} xml",e)
    #     return
    # handle_text.close()
    print(f"  Traitement {NC}")
    list_file = []
    prev_buffer_ecriture = ""
    buffer_ecriture = ""
    found_anything = False
    set_feature_key = set()
    set_expected_feature = set(["Aucun","CDS","centromere","intron","mobile_element","ncRNA","rRNA","telomere","tRNA","3'UTR","5'UTR"])
    for i in range(len(record[0]["GBSeq_feature-table"])):
        feature_location = record[0]["GBSeq_feature-table"][i]["GBFeature_location"]
        feature_key = record[0]["GBSeq_feature-table"][i]["GBFeature_key"]
        set_feature_key.add(feature_key)
        if feature_key != selected_region and not (selected_region == "intron" and feature_key == "CDS"):
            continue

        nb_region_found = number_region_found.get()
        nb_region_found += 1
        number_region_found.put(nb_region_found)
        NC_filename = selected_region + "_" + str(name) + "_" + str(NC) + ".txt"
        
        if len(list_file) != 0 :
            if NC_filename not in list_file:
                os.remove(path + name + "/" + NC_filename)
                list_file.append(NC_filename)
        else :
            if os.path.isfile(path + name + "/" + NC_filename):
                os.remove(path + name + "/" + NC_filename)
            list_file.append(NC_filename)
        
        # Parsing
        buffer_ecriture = ""
        
        header_str = selected_region + ' ' + name + ' ' + NC + ': '
        
        semaphore.acquire()
        if feature_location.find("complement") != -1 and feature_location.find("join") != -1:
            print("complement join")
            # buffer_ecriture = join(buffer_ecriture, header_str, selected_region, feature_location, record_fasta, True)

        elif feature_location.find("complement") != -1:
            print("complement")
            # buffer_ecriture = extract(buffer_ecriture, header_str, selected_region, feature_location, record_fasta, True)

        elif feature_location.find("join") != -1:
            print("join")
            # buffer_ecriture = join(buffer_ecriture, header_str, selected_region, feature_location, record_fasta, False)

        else:
            print("normal")
            # buffer_ecriture = extract(buffer_ecriture, header_str, selected_region, feature_location, record_fasta, False)
        semaphore.release()

        if NC == "NC_071382":
            truc = "aaaaaaaaaaaaaaaa"
            truc = 0

        # Si on doit écrire quelque chose on créé le fichier correspondant
        if buffer_ecriture != "" and buffer_ecriture != prev_buffer_ecriture:
            with open(path + name + "/" + NC_filename, 'a+') as out:
                out.write(buffer_ecriture)
            prev_buffer_ecriture = buffer_ecriture
            found_anything = True
        else:
            nb_region_found = number_region_found.get()
            nb_region_found -= 1
            number_region_found.put(nb_region_found)
    print(f"  Fin traitement {NC}")
    intersec = sorted(list(set_expected_feature.intersection(set_feature_key)))
    if intersec and False:
        return (f"fichier {NC_filename} créé","green") if found_anything else (f'[{selected_region}] non trouvée (présentes dans ce NC : {", ".join(intersec)})',"orange")
    else:
        return (f"fichier {NC_filename} créé","green") if found_anything else (f'[{selected_region}] non trouvée (aucune région connue dans ce NC)',"orange")


a=f2(NC_list)