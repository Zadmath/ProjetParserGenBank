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
from window import FENETRE
from genericpath import *
from Bio import Entrez
from Bio import SeqIO

def load_dataframe():
    #on ouvre le fichier pickle contenant les informations sur les organismes
    organism_df = pd.read_pickle("../pickle/organism_df")
    return organism_df

# organism_df = pd.read_pickle("../pickle/organism_df")
# df = pd.read_pickle("../pickle/organism_df")
# org="Campylobacter_coli"
# line = organism_df[organism_df["name"].str.replace(" ","_").replace("[","_").replace("]","_").replace(":","_").replace("/","_",regex=True) == org]
# df["nb_NC"] = df["NC"].apply(len)




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

def join(buffer_ecriture, header_str, selected_region, feature_location, record_fasta, is_complement):
    # if selected_region != "intron":
    #     buffer_ecriture += header_str + feature_location
    
    if is_complement:
        feature_location = feature_location[16:-2]
    else:
        feature_location = feature_location[5:-1]
    
    x = feature_location.split(",")

    is_valid = True

    indexes = []
    for xi in x:
        xi = xi.split("..")
        try:
            indexes.append([int(xi[0]),int(xi[1])])
        except:
            is_valid = False

    # check si les index sont dans le bon ordre
    for i in range(len(indexes) - 1):
        if indexes[i] > indexes[i + 1]:
            return ""
    indexes.sort(key=lambda r:r[0])

    # reconstruction de la localisation, correctement mise en forme
    # TODO changer formattage
    str_loc = " , ".join( f'{deb}..{fin}' for deb,fin in indexes )
    if len(indexes) > 1:
        str_loc = f"join( {str_loc} )"
    if is_complement:
        str_loc = f"complement( {str_loc} )"
    feature_location = str_loc

    if selected_region == "intron":
        fn = []
        for i in range(len(indexes) - 1):
            if(check_inf_sup(indexes[i][1],indexes[i+1][0]) == False):
                is_valid = False
            try :
                fn.append(FeatureLocation(indexes[i][1], indexes[i+1][0]-1))
                # fn.append(FeatureLocation(indexes[i][1]-1, indexes[i+1][0]))
            except :
                return ""
        if not is_valid:
            return ""
    else:
        buffer_ecriture += header_str + feature_location
        fn = []
        for xi in indexes:
            if(check_inf_sup(xi[0],xi[1]) == False):
                is_valid = False
            else :
                fn.append(FeatureLocation(xi[0]-1, xi[1]))
        if not is_valid:
            return ""
    
    if len(fn) > 1:
        f = CompoundLocation(fn)
    else:
        f = SeqFeature(fn[0], type="domain")
    
    if selected_region != "intron":
        if is_complement:
            buffer_ecriture += "\n" +str(f.extract(record_fasta.seq).reverse_complement())
        else:
            buffer_ecriture += "\n" +str(f.extract(record_fasta.seq))
        buffer_ecriture += "\n"
    
    if len(fn) >= 1:
        for i in range(len(fn)):
            type_seq = " Exon "
            if selected_region == "intron":
                type_seq = " Intron "
            if is_complement:
                if selected_region=="intron":
                    buffer_ecriture += header_str + feature_location + type_seq+ str(i+1)+ "\n" + str(SeqFeature(fn[len(fn)-i-1], type="domain").extract(record_fasta.seq).reverse_complement())
                else:
                    buffer_ecriture += header_str + feature_location + type_seq+ str(i+1)+ "\n" + str(SeqFeature(fn[len(fn)-i-1], type="domain").extract(record_fasta.seq).reverse_complement())
            else:
                if selected_region=="intron":
                    buffer_ecriture += header_str + feature_location + type_seq+ str(i+1)+ "\n" + str(SeqFeature(fn[i], type="domain").extract(record_fasta.seq))
                else:
                    buffer_ecriture += header_str + feature_location + type_seq+ str(i+1)+ "\n" + str(SeqFeature(fn[i], type="domain").extract(record_fasta.seq))
            buffer_ecriture += "\n"
    return buffer_ecriture

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
    print("DANS LA FONCTION")
    print(selected_region)
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
        print("selected_region",selected_region)
        print("feature_key",feature_key)
        if feature_key != selected_region and not (selected_region == "intron" and feature_key == "CDS"):
            continue
        print("TRAITEMENT DE",feature_key) 
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
    
        print(type(selected_region))
        header_str = selected_region + ' ' + name + ' ' + NC + ': '
        
        if feature_location.find("complement") != -1 and feature_location.find("join") != -1:
            buffer_ecriture = join(buffer_ecriture, header_str, selected_region, feature_location, record_fasta, True)

        elif feature_location.find("complement") != -1:
            buffer_ecriture = extract(buffer_ecriture, header_str, selected_region, feature_location, record_fasta, True)

        elif feature_location.find("join") != -1:
            buffer_ecriture = join(buffer_ecriture, header_str, selected_region, feature_location, record_fasta, False)

        else:
            buffer_ecriture = extract(buffer_ecriture, header_str, selected_region, feature_location, record_fasta, False)

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

