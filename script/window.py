#!/usr/bin/env python3
from tkinter import *
from tkinter import ttk
import tkinter.font as tkFont
from ttkwidgets import CheckboxTreeview
# from ttkthemes import ThemedTk
import pandas as pd

import os
import sys
import time
import ctypes

from PIL import Image as Img

import create_arbo as arbo

try:
    ctypes.windll.shcore.SetProcessDpiAwareness(True) # améliore la netteté de l'app
except:
    pass #ne marche pas sous linux

#Couleurs
BG_COLOR = "#292a3d"
TXT_COLOR= "#E0FFE0" # "white"
DISABLED_SCROLLBAR= "#191a2d"
HOVER_COLOR = "#595a6d"
STOP_COLOR = "#ff3030"
HOVER_STOP_COLOR = "#ff0000"
VERT_FONCE = "#00ff00"
ANTHRACITE = "#191a2d"


def collect_data() :
    #Création de la fenetre de bienvenue
    root = Tk()
    root.title("Parser GENBANK")

    #Affichage texte
    lbl = Label(root, text="Traitement des fichiers...",foreground=TXT_COLOR,background=BG_COLOR)
    lbl.pack(side = LEFT, anchor=  W, padx = 15, pady = 15)
    #Configuration du style de la progressbar
    style = ttk.Style(root)
    style.configure("start.Horizontal.TProgressbar", troughcolor=ANTHRACITE, background=TXT_COLOR, thickness=20)
    root.configure(bg=ANTHRACITE )

    
    #Création progressbar
    progressbar = ttk.Progressbar(root,style="start.Horizontal.TProgressbar", orient = HORIZONTAL, length = 100, mode = 'determinate')
    progressbar.pack( padx = 15, pady = 15)
    root.after(200,  lambda : arbo.collect_ids(progressbar,root, lbl = lbl))
    root.mainloop()

class Affichage:
    def __init__(self):
        collect_data()

if __name__ == "__main__":
    App = Affichage()
