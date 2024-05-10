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
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

from genericpath import *
from Bio import Entrez
from Bio import SeqIO

def load_dataframe():
    #on ouvre le fichier pickle contenant les informations sur les organismes
    organism_df = pd.read_pickle("../pickle/organism_df")
    return organism_df

organism_df = pd.read_pickle("../pickle/organism_df")
df = pd.read_pickle("../pickle/organism_df")
org="Campylobacter_coli"
line = organism_df[organism_df["name"].str.replace(" ","_").replace("[","_").replace("]","_").replace(":","_").replace("/","_",regex=True) == org]
df["nb_NC"] = df["NC"].apply(len)




# Split a path in a drive specification (a drive letter followed by a
# colon) and the path specification.
# It is always true that drivespec + pathspec == p
def splitdrive(p):
    """Split a pathname into drive/UNC sharepoint and relative path specifiers.
    Returns a 2-tuple (drive_or_unc, path); either part may be empty.

    If you assign
        result = splitdrive(p)
    It is always true that:
        result[0] + result[1] == p

    If the path contained a drive letter, drive_or_unc will contain everything
    up to and including the colon.  e.g. splitdrive("c:/dir") returns ("c:", "/dir")

    If the path contained a UNC path, the drive_or_unc will contain the host name
    and share up to but not including the fourth directory separator character.
    e.g. splitdrive("//host/computer/dir") returns ("//host/computer", "/dir")

    Paths cannot contain both a drive letter and a UNC path.

    """
    p = os.fspath(p)
    if len(p) >= 2:
        if isinstance(p, bytes):
            sep = b'\\'
            altsep = b'/'
            colon = b':'
            unc_prefix = b'\\\\?\\UNC\\'
        else:
            sep = '\\'
            altsep = '/'
            colon = ':'
            unc_prefix = '\\\\?\\UNC\\'
        normp = p.replace(altsep, sep)
        if normp[0:2] == sep * 2:
            # UNC drives, e.g. \\server\share or \\?\UNC\server\share
            # Device drives, e.g. \\.\device or \\?\device
            start = 8 if normp[:8].upper() == unc_prefix else 2
            index = normp.find(sep, start)
            if index == -1:
                return p, p[:0]
            index2 = normp.find(sep, index + 1)
            if index2 == -1:
                return p, p[:0]
            return p[:index2], p[index2:]
        if normp[1:2] == colon:
            # Drive-letter drives, e.g. X:
            return p[:2], p[2:]
    return p[:0], p


# Split a path in head (everything up to the last '/') and tail (the
# rest).  After the trailing '/' is stripped, the invariant
# join(head, tail) == p holds.
# The resulting head won't end in '/' unless it is the root.

def join(path, *paths):
    path = os.fspath(path)
    if isinstance(path, bytes):
        sep = b'\\'
        seps = b'\\/'
        colon = b':'
    else:
        sep = '\\'
        seps = '\\/'
        colon = ':'
    try:
        if not paths:
            path[:0] + sep  #23780: Ensure compatible data type even if p is null.
        result_drive, result_path = splitdrive(path)
        for p in map(os.fspath, paths):
            p_drive, p_path = splitdrive(p)
            if p_path and p_path[0] in seps:
                # Second path is absolute
                if p_drive or not result_drive:
                    result_drive = p_drive
                result_path = p_path
                continue
            elif p_drive and p_drive != result_drive:
                if p_drive.lower() != result_drive.lower():
                    # Different drives => ignore the first path entirely
                    result_drive = p_drive
                    result_path = p_path
                    continue
                # Same drive in different case
                result_drive = p_drive
            # Second path is relative to the first
            if result_path and result_path[-1] not in seps:
                result_path = result_path + sep
            result_path = result_path + p_path
        ## add separator between UNC and non-absolute path
        if (result_path and result_path[0] not in seps and
            result_drive and result_drive[-1:] != colon):
            return result_drive + sep + result_path
        return result_drive + result_path
    except (TypeError, AttributeError, BytesWarning):
        genericpath._check_arg_types('join', path, *paths)
        raise

def check_inf_sup(inf,sup):
    """Fonction simple pour savoir qui est inferieur a qui

    Args:
        inf: value
        sup: value

    Returns:
        true if inf is inferior or equal to sup
        false if not
    """
    if(inf <= sup):
        return True
    else:
        return False


def extract(buffer_ecriture, header_str, selected_region, feature_location, record_fasta, is_complement):
    if selected_region == "intron":
        return buffer_ecriture
    if is_complement:
        feature_location = feature_location[11:-1]
    x = feature_location.split("..")
    try:
        (int(x[0]),int(x[1]))
    except ValueError:
        return buffer_ecriture
    else:
        if(check_inf_sup(int(x[0]),int(x[1])) == False):
            return buffer_ecriture
        if is_complement: #TODO changer formattage
            buffer_ecriture = header_str + f"complement( {x[0]}..{x[1]} )\n"
        else:
            buffer_ecriture = header_str + f"{x[0]}..{x[1]}\n" 
        f = SeqFeature(FeatureLocation(int(x[0])-1, int(x[1])), type="domain")
        if is_complement:
            buffer_ecriture += str(f.extract(record_fasta.seq).reverse_complement())
        else:
            buffer_ecriture += str(f.extract(record_fasta.seq))
    buffer_ecriture += "\n"
    return buffer_ecriture

def f2(number_region_found, number_region_already_found, path, NC, name, selected_region):
    name = name.replace(" ", "_")
    name = name.replace("[", "_")
    name = name.replace("]", "_")
    name = name.replace(":", "_")
    Entrez.email = ''.join(random.choice(string.ascii_lowercase) for i in range(20))+'@'+''.join(random.choice(string.ascii_lowercase) for i in range(20))+ '.com'

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

        # nb_region_found = number_region_found.get()
        # nb_region_found += 1
        # number_region_found.put(nb_region_found)
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
        
        # semaphore.acquire()
        if feature_location.find("complement") != -1 and feature_location.find("join") != -1:
            buffer_ecriture = join(buffer_ecriture, header_str, selected_region, feature_location, record_fasta, True)

        elif feature_location.find("complement") != -1:
            buffer_ecriture = extract(buffer_ecriture, header_str, selected_region, feature_location, record_fasta, True)

        elif feature_location.find("join") != -1:
            buffer_ecriture = join(buffer_ecriture, header_str, selected_region, feature_location, record_fasta, False)

        else:
            buffer_ecriture = extract(buffer_ecriture, header_str, selected_region, feature_location, record_fasta, False)
        # semaphore.release()

        if NC == "NC_071382":
            truc = "aaaaaaaaaaaaaaaa"
            truc = 0

        # Si on doit écrire quelque chose on créé le fichier correspondant
        if buffer_ecriture != "" and buffer_ecriture != prev_buffer_ecriture:
            print("path + name + '/' + NC_filename",path + name + "/" + NC_filename)
            with open(path + name + "/" + NC_filename, 'a+') as out:
                out.write(buffer_ecriture)
            prev_buffer_ecriture = buffer_ecriture
            found_anything = True
        # else:
            # nb_region_found = number_region_found.get()
            # nb_region_found -= 1
            # number_region_found.put(nb_region_found)
    print(f"  Fin traitement {NC}")
    intersec = sorted(list(set_expected_feature.intersection(set_feature_key)))
    if intersec and False:
        return (f"fichier {NC_filename} créé","green") if found_anything else (f'[{selected_region}] non trouvée (présentes dans ce NC : {", ".join(intersec)})',"orange")
    else:
        return (f"fichier {NC_filename} créé","green") if found_anything else (f'[{selected_region}] non trouvée (aucune région connue dans ce NC)',"orange")

# nb_region_found = 0
# nb_region_already_downloaded = 0
# number_region_found = 0
# number_region_already_found = 0
# selected_region = "CDS"
# name, path, NC_list = line.iloc[0].values
# print(name)
# print(path)
# print(NC_list)

# nb_region_found = 0
# nb_region_already_downloaded = 0
# for NC in NC_list:
#     a=f2(number_region_found, number_region_already_found, path, NC, name, selected_region)
#     print(a)

