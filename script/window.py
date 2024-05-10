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


class Affichage:
    def __init__(self):
        #Création de la fenetre de bienvenue
        self.root = Tk()
        self.root.title("Parser GENBANK")
        
        #Affichage texte
        self.lbl = Label(self.root, text="Traitement des fichiers...",foreground=TXT_COLOR,background=BG_COLOR)
        self.lbl.pack(side = LEFT, anchor=  W, padx = 15, pady = 15)
        #Configuration du style de la progressbar
        self.style = ttk.Style(self.root)
        self.style.configure("start.Horizontal.TProgressbar", troughcolor=ANTHRACITE, background=TXT_COLOR, thickness=20)
        self.root.configure(bg=ANTHRACITE )
        
        #Création progressbar
        self.progressbar = ttk.Progressbar(self.root,style="start.Horizontal.TProgressbar", orient = HORIZONTAL, length = 100, mode = 'determinate')
        self.progressbar.pack( padx = 15, pady = 15)
        self.root.after(200,  lambda : arbo.collect_ids(self.progressbar,self.root, lbl = self.lbl))
        self.root.mainloop()
        
    def get_root(self):
        return self.root
    
    def get_progressbar(self):
        return self.progressbar
    
    def get_style(self):
        return self.style

if __name__ == "__main__":
    App = Affichage()
