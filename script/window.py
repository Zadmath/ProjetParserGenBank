#!/usr/bin/env python3
from tkinter import *
from tkinter import ttk
import tkinter.font as tkFont
from ttkwidgets import CheckboxTreeview
import traitement_ids as fetch
# from ttkthemes import ThemedTk
import pandas as pd
from multiprocessing import Pool, Process, Manager, Queue,cpu_count, Semaphore, Value, Event
import traceback
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
cwd = os.path.dirname(os.path.realpath(__file__))
#Couleurs
BG_COLOR = "#727387"
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



class FENETRE:
    def __init__(self):
        # attribute
        #si le fichier pickle dans /pickle existe, on le charge
        #is.path.exists detecte un fichier supprimé
        #on affiche les fichiers dans le repertoire courant 
        for file in os.listdir("."):
            print(file)
        if os.path.exists("../Results/"):
            print("Fichier pickle trouvé")
        else:
            print("Fichier pickle non trouvé")
            Affichage()
        if not os.path.exists("../Results/"):
            Affichage()
        self.organism_df : pd.DataFrame = fetch.load_dataframe()

        # Creation de la fenetre
        self.window = Tk()
        # self.window = ThemedTk(theme="azure")
        self.window.geometry("1010x550")
        self.window.resizable(width=False,height=False)
        self.window.title("Parseur GenBank")
        self.window.configure(background="#292a3d")

        # self.window.tk.call("source","azure.tcl")
        # self.window.tk.call("source",r"C:\Users\sim20\AppData\Local\Programs\Python\Python311\Lib\site-packages\ttkthemes\gif\adapta\adapta.tcl")
        # self.window.tk.call("set_theme","dark")
        #Configuration des styles
        style = ttk.Style(self.window)
        
        style.theme_use("clam")

        self.indent_treeview = 10
        self.font_treeview = (None,9)
        style.configure("Treeview", 
                        font=self.font_treeview,
                        background=BG_COLOR, 
                        fieldbackground=BG_COLOR,
                        foreground=TXT_COLOR,
                        indent=self.indent_treeview)
        style.map("Treeview.Heading",
                background = [('active', BG_COLOR)])
        style.configure('Treeview.Heading',
                        font=(None,11),
                        background=BG_COLOR,
                        foreground=TXT_COLOR,
                        # activebackground=BG_COLOR,
                        relief='flat')
        # style.configure("Treeview.Cell",
        #                 # indicatormargins=1,
        #                 padding=10
        #                 )

        style.configure("start.Horizontal.TProgressbar",
                        troughcolor=BG_COLOR)
        style.map("TScrollbar", 
                background =[("!disabled",TXT_COLOR),("disabled", DISABLED_SCROLLBAR)],
                darkcolor  =[("!disabled",BG_COLOR ),("disabled", BG_COLOR)],
                lightcolor =[("!disabled",TXT_COLOR),("disabled", DISABLED_SCROLLBAR)],
                troughcolor=[("!disabled",BG_COLOR ),("disabled", BG_COLOR)],
                bordercolor=[("!disabled",BG_COLOR ),("disabled", BG_COLOR)],
                arrowcolor =[("!disabled",BG_COLOR ),("disabled", BG_COLOR)],
                gripcount  =[("!disabled",0),("disabled",0)])

        style.map("TMenubutton",
                        background=[("!active",BG_COLOR),("active",HOVER_COLOR)],
                        foreground=[("!active",TXT_COLOR),("active",TXT_COLOR)],
                        arrowcolor=[("!active",TXT_COLOR),("active",TXT_COLOR)],
                        # borderwidth=[("!active",100), ("active",100)]
                        )
        # style.configure("TMenubutton",
        #                 bordercolor="red")
        self.window.option_add("*Menu",BG_COLOR)
        self.window.option_add("*Menu.background",BG_COLOR)
        self.window.option_add("*Menu.foreground",TXT_COLOR)
        self.window.option_add("*Menu.selectColor",TXT_COLOR)
        self.window.option_add("*Menu.activeBackground",HOVER_COLOR)
        self.window.option_add("*Menu.relief","flat")

        print(style.element_options("Treeview.Heading"))

        self.tree_array = []
        self.org_selected = []
        # Zone 1 : Liste des fichiers
        (self.treeview, self.scrollbar, self.scrollbar_h) = self.create_tree()
        
        self.treeview.bind("<Button-1>", self.box_click, True)
        self.treeview.bind("<ButtonRelease-1>", self.after_click, True)
        self.treeview.bind("<Double-Button-1>", self.open_file, True)

        # Zone 2 : Informations
        (self.Info, self.Info_lab1, self.label_progress, self.label_current_search, self.progressbar, self.button_arret) = self.create_info()
        # self.label_progress.config(text="Aucune recherche en cours")
        # Zone 3 : Logs
        (self.Logs, self.Logs_lab, self.log_text, self.log_scroll) = self.create_log()

        # Menu de selection des regions
        (self.selected_region, self.menu_menu, self.run_search, self.OptionList_region) = self.create_region_menu()
        print(self.selected_region)
        list = Frame(self.window)
        list.grid(row=3,column=1)
        # Lancement

        self.window.mainloop()


    # Methode de création de l'arbre
    def tree_exist(self, tree, name) :
        if tree.exists(name) :
            name += 'AZERTY'
            # print(f"tree_exist {name = }")
            return self.tree_exist(tree, name)
        return name

    def create_node(self, tree : CheckboxTreeview, node_path, node, depth=0) :
        name_nodes = os.listdir(node_path)
        name_nodes.sort()


        total = 0
        for name in name_nodes :
            name_path = os.path.join(node_path, name)
            if os.path.isdir(name_path):
                # iid = self.tree_exist(tree, name)
                iid = name_path
                # print(iid)
                try :
                    tree.insert(node, '1000000',iid=iid, text=name)
                    self.tree_array.append(name)
                    # name_path += '/'
                    self.create_node(tree, name_path, iid, depth=depth+1)
                except Exception as e:
                    # print(str(e))
                    pass
            else:
                pass 
                # print(f"{name_path} n'est pas une dir (?!)")

    # Méthode de création de l'abre graphique
    def create_tree(self):
        list = Frame(self.window)
        list.place(x=5, y=5, anchor="nw", width=270, height=540)
        
        # Création de la scrollbar
        # scrollbar = Scrollbar(list,background="#545577",activebackground="#464768")
        scrollbar = ttk.Scrollbar(list)
        scrollbar.pack( side = RIGHT, fill = Y )

        scrollbar_h = ttk.Scrollbar(list,orient="horizontal")
        scrollbar_h.pack( side = BOTTOM, fill = X )


        # Création de la fenetre treeview
        treeview = CheckboxTreeview(list)
        treeview.column("#0", stretch=False, minwidth=250)
        treeview.heading('#0', text='Arborescence des fichiers', anchor="w")
        treeview.configure(yscrollcommand=scrollbar.set)
        treeview.configure(xscrollcommand=scrollbar_h.set)

        # os.chdir(os.path.dirname(os.path.realpath(__file__)))
        cwd = os.getcwd()
        im_unchecked = Img.open(os.path.join(cwd,"img_theme","unchecked.png"))
        im_tristate = Img.open(os.path.join(cwd,"img_theme","tristate.png"))
        im_checked = Img.open(os.path.join(cwd,"img_theme","checked.png"))

        for im in [im_unchecked, im_tristate, im_checked]:
            pixels = im.load()
            old_color = pixels[0,0]
            new_color = tuple([ int(x,16) for x in (BG_COLOR[1:3],BG_COLOR[3:5],BG_COLOR[5:7]) ])
            # print(new_color)
            # Parcourir tous les pixels de l'image
            for i in range(im.width):
                for j in range(im.height):
                    # Vérifier si le pixel a la couleur à remplacer
                    if pixels[i, j] == old_color:
                        # Remplacer la couleur du pixel par la nouvelle couleur
                        pixels[i, j] = new_color

        treeview.im_unchecked.paste(im_unchecked)
        treeview.im_tristate.paste(im_tristate)
        treeview.im_checked.paste(im_checked)
        treeview.pack(fill="both", expand="yes")

        scrollbar.configure(command=treeview.yview)
        scrollbar_h.configure(command=treeview.xview)
        #Generation de l'arbre
        root_path = "../Results/"
        item = treeview.insert('', '0', text='Results', iid='Results')

        # treeview.change_state(item, "disabled") # retire la checkbox
        
        self.tree_array.append('Results')


        self.create_node(treeview, root_path, 'Results')

        # suppression des checkboxs autres que les leafs
        self.remove_boxes_for_nonleaf(treeview)

        

        treeview.column("#0", width=460)

        return treeview, scrollbar, scrollbar_h

    def remove_boxes_for_nonleaf(self, treeview=None):
        return
        if not treeview:
            treeview = self.treeview
        for node in self.tree_array:
            if len(treeview.get_children(node)):
                treeview.change_state(node, "disabled")

    def box_click(self, event):
        """ check or uncheck box when clicked """
        return
        x, y, widget = event.x, event.y, event.widget
        elem = widget.identify("element", x, y)
        if "image" in elem:
            # a box was clicked
            item = self.treeview.identify_row(y)
            tags = self.treeview.item(item, "tags")
            if ("unchecked" in tags) or ("tristate" in tags):
                try : 
                    self.org_selected.remove(item)
                except : 
                    pass
            else:
                self.org_selected.append(item)
        self.remove_boxes_for_nonleaf()
        print(self.org_selected)

    def open_file(self, event):
        # print("test")
        x, y, widget = event.x, event.y, event.widget
        elem = widget.identify("element", x, y)
        # print(elem)
        if "text" in elem:
            # a box was clicked
            item_id = self.treeview.identify_row(y)
            item_str = self.treeview.item(item_id, "text")
            tags = self.treeview.item(item_id, "tags")
            # print(tags)
            if ("checked" in tags):
                # print("OPEN CESAME")
                # print(os.getcwd())
                parents = []
                current_id = self.treeview.parent(item_id)
                if current_id == '':
                    fullpath = os.path.join( os.getcwd(), "..", "Results" )
                else:
                    current_str = self.treeview.item(current_id, "text")
                    while current_id != self.treeview.get_children()[0]:
                        parents.append(current_str)
                        current_id = self.treeview.parent(current_id)
                        current_str = self.treeview.item(current_id, "text")
                    fullpath = os.path.join( os.getcwd(), "..", current_str, *list(reversed(parents)), item_str )
                print(fullpath)
                if sys.platform=='win32':
                    os.startfile(fullpath)
                else:
                    os.system(f"xdg-open \"{fullpath}\"")
    

    def after_click(self,evt=None):

        return
        font = tkFont.Font(size=self.font_treeview[1])
        maxwidth = -1
        maxitem = None
        def explore_shown(item, depth=0):
            nonlocal maxwidth
            nonlocal maxitem
            width = font.measure(str(item)) + depth*self.indent_treeview + 18 + (16 if len(self.treeview.get_children(item)) == 0 else 0) + 5
            # print(" "*depth*2 + str(item) + f"{width}")
            if width > maxwidth:
                maxwidth = width
                maxitem = str(item)
            # print(self.treeview.item(item))
            if self.treeview.item(item)["open"]:
                if depth != 3:
                    for c in self.treeview.get_children(item):
                        explore_shown(c, depth+1)
                else:
                    maxwidth = max(maxwidth, self.item_lengths[str(item)])
        
        explore_shown(self.treeview.get_children()[0])
        if maxwidth != -1:
            print(f'{maxwidth} {maxitem}')
            self.treeview.column("#0", width=max(252,maxwidth))


    # INFO CREATION METHODS
    def create_info(self):
        Info = Frame(self.window)
        Info.place(x=275, y=5, anchor="nw", width=730, height=95)

        Info_lab1 = LabelFrame(Info, text="Selectionner une région",font=(None,11), padx=10, pady=10,foreground=TXT_COLOR,background=BG_COLOR)
        Info_lab1.pack(fill="both", expand="yes")

        label_progress=Label(Info_lab1,text="0/0",font=(None,9),foreground=TXT_COLOR,background=BG_COLOR)
        label_progress.grid(row=1,column=1,sticky=W,pady=5,padx=10)

        label_current_search = Label(Info_lab1, text="Aucune recherche en cours",font=(None,9), foreground=TXT_COLOR,background=BG_COLOR)
        # label_current_search = Label(Info_lab1, text="CDS de 'Homo Sapiens'",font=(None,9), foreground="white",background="#292a3d")
        label_current_search.grid(row=1, column=2, sticky=W, pady=5, padx=10)

        Label(Info_lab1, justify = LEFT, text="Région fonctionnelle choisie : ",font=(None,9),foreground=TXT_COLOR,background=BG_COLOR).grid(row = 0, column = 0, sticky = NW,pady=5)
        #Création progressbar
        progressbar = ttk.Progressbar(Info_lab1,style="start.Horizontal.TProgressbar", orient = HORIZONTAL, length = 200, mode = 'determinate')
        progressbar.grid(row=1,column=0,sticky=W,pady=5)


        # bouton d'arrêt
        button_arret = Button(self.window, text ="Arrêter la recherche",font=(None,9), command = self.interrupt_search, borderwidth=1,foreground=TXT_COLOR,background=STOP_COLOR,activebackground=HOVER_STOP_COLOR)
        
        # button_arret.grid(row = 0, column = 3,sticky=NW,padx=10)      
        button_arret.place(x=840,y=30,anchor="nw")

        button_arret.config(state=DISABLED)
        button_arret.config(background=HOVER_COLOR)


        return Info, Info_lab1, label_progress, label_current_search, progressbar, button_arret

    def interrupt_search(self):
        self.interrupt_signal = True
        self.print_on_window("Interruption !","warning")

    # LOG CREATION METHODS
    def create_log(self):
        Logs = Frame(self.window)
        Logs.place(x=275, y=100, anchor="nw", width=730, height=445)
        Logs_lab = LabelFrame(Logs, text="Logs",font=(None,11),foreground=TXT_COLOR,background=BG_COLOR)
        Logs_lab.pack(fill="both", expand="yes")


        log_text = Text(Logs_lab, wrap=WORD,font=(None,9), foreground=TXT_COLOR,background="#001000")
        log_text.pack(fill="both", expand="yes",side=LEFT)
        log_text.tag_config('warning',foreground="red")
        log_text.tag_config('red', foreground="#e69138")
        log_text.tag_config('green', foreground="green")
        log_text.tag_config('white', foreground="white")
        log_text.tag_config('orange', foreground="orange")


        # log_scroll = Scrollbar(Logs_lab, command = log_text.yview, background="#545577",activebackground="#464768")
        log_scroll = ttk.Scrollbar(Logs_lab, command = log_text.yview)
        log_scroll.pack( side = RIGHT, fill = Y )


        log_text.configure(yscrollcommand=log_scroll.set)

        def on_mousewheel(event):
            log_text.yview_scroll(-1*(event.delta//120), "units")

        log_text.bind("<MouseWheel>", on_mousewheel)

        log_text.bindtags((str(log_text), str(Logs_lab), "all"))
        return Logs, Logs_lab, log_text, log_scroll

    # REGION MENU CREATION METHODS
    def create_region_menu(self):
        OptionList = [
        "Aucun",
        "CDS",
        "centromere",
        "intron",
        "mobile_element",
        "ncRNA",
        "rRNA",
        "telomere",
        "tRNA",
        "3'UTR",
        "5'UTR"
        ]

        variable = StringVar(self.Info_lab1)
        variable.set(OptionList[1]) #par défaut : CDS

        menu = ttk.OptionMenu(self.Info_lab1, variable, variable.get(), *OptionList)


        menu.grid(row = 0, column = 1, sticky = NW)
        variable.trace("w", self.callback)

        run_search = Button(self.Info_lab1, text ="Recherche",font=(None,9), command = self.search_button_callback, borderwidth=1,foreground=TXT_COLOR,background=BG_COLOR,activebackground=HOVER_COLOR)
        run_search.grid(row = 0, column = 2,sticky=NW,padx=10)      
        
        return variable, menu, run_search, OptionList

    # FEATURES METHODS

    def update_tree_tags(self):
        for node in self.tree_array:
            current_path = self.get_path(node).replace("Other1", "Other").replace("unclassified1", "unclassified").replace("uncultured_bacterium1", "uncultured_bacterium")
            is_leaf = False
            for (index, name, path, NC_list) in self.organism_df.itertuples():
                path_full = path + name.replace(" ", "_").replace("[", "_").replace("]", "_").replace(":", "_") + '/'
                if path_full == current_path:
                    is_leaf = True
                    if os.listdir(path_full) == []:
                        self.treeview.item(node)
                    else:
                        self.treeview.item(node)
                    break
            if is_leaf:
                continue

        self.print_on_window("Arborescence mise a jour",'white')

    def print_on_window(self, t, couleur): #affiche t dans les logs
        time_string = time.strftime('%H:%M:%S')
        self.log_text.insert(INSERT, time_string + ' : ' + str(t) + "\n",couleur)
        self.log_text.yview(END)

    def callback(self, *args): # fonction pour executer du code pour le menu, a changer
        self.print_on_window("Région choisie : " + self.selected_region.get(),"white")
        self.menu_menu["text"] = self.selected_region.get()
    def get_path(self, current_item):
        path = '/' + current_item + '/'
        while len(self.treeview.parent(current_item)) != 0 :
            parent = self.treeview.parent(current_item)
            current_item = parent
            path = '/' + parent + path
        return ('..' + path)

    def search_button_callback(self): # Fonction boutton
        
        
        self.org_selected = [self.treeview.item(x)["text"] for x in self.treeview.get_checked()]


        self.print_on_window("------------------Début de la recherche-------------------","white")
        self.window.update()


        if self.org_selected == [] :
            self.print_on_window("Aucun organisme choisi","warning")
            self.window.update()
            return
        elif self.selected_region.get() == 'Aucun':
            self.print_on_window("Aucune région choisie","warning")
            self.window.update()
            return
        else :
            self.run_search.config(state=DISABLED)
            self.run_search.config(text="Recherche en cours...")
            self.button_arret.config(state=NORMAL)
            self.button_arret.config(background=STOP_COLOR)

            self.interrupt_signal = False

            self.progressbar['value']=0
            self.label_progress['text']="0/"+str(len(self.org_selected))
            self.window.update_idletasks()
            org_done = []
            c = 0
            count=0
            for i,org in enumerate(self.org_selected) :
                if self.interrupt_signal:
                    break
                
                nb_region_found = 0
                nb_region_found_1 = 0
                nb_region_found_2 = 0
                nb_region_already_downloaded_2 = 0
                nb_region_already_downloaded_1 = 0
                line = self.organism_df[self.organism_df["name"].str.replace(" ","_").replace("[","_").replace("]","_").replace(":","_").replace("/","_",regex=True) == org]
                try:
                    name, path, NC_list = line.iloc[0].values
                    index = line.index[0]

                    self.label_current_search.config(text=f"{self.selected_region.get()} de '{name}' ({i+1}/{len(self.org_selected)} organisme{'s' if len(self.org_selected) > 1 else ''})")
                    self.progressbar.config(value=0)
                    self.label_progress.config(text=f"0/{len(NC_list)}")
                    self.print_on_window(f"Organisme {i+1}/{len(self.org_selected)}: {name} ({len(NC_list)} NC à traiter)", "white")
                    

                    nb_new_region_found = -1

                    print(f"1 {index} {name} {path} {NC_list} {len(NC_list)}")
                    c += 1
                    # if(self.selected_region.get() == 'CDS'):
                    #     nb_region_found_1, nb_region_already_downloaded_1 = fetch.load_data_from_NC(index, name, path, NC_list, 'intron')
                    self.window.update()
                    #on affiche le type de self.selected_region
                    print(self.selected_region.get())
                    nb_region_found = 0
                    nb_region_already_downloaded = 0
                    # Retire les duplicatas
                    NC_list = list(dict.fromkeys(NC_list))

                    myStack = Queue()
                    number_region_found = Queue()
                    number_region_found.put(nb_region_found)
                    number_region_already_found = Queue()
                    number_region_already_found.put(nb_region_found)
                    nb_NC_done = Value('i', 0)
                    to_log = Queue()
                    interrupt_event = Event()


                    for NC in NC_list:
                        myStack.put(NC)          
                    print("STACK")   
                    for NC in NC_list:               
                        nb_region_found_2, nb_region_already_downloaded_2 = fetch.f2(number_region_found, number_region_already_found, path, NC, name, self.selected_region.get())
                except:
                    #on affiche l'erreur 
                    error_message = traceback.format_exc()
                    print(error_message)
                    

                    self.print_on_window(f"Erreur : l'organisme {i+1}/{len(self.org_selected)} [{org}] n'a aucun NC à traiter !","warning")

                self.window.update()

                self.print_on_window("","white")
            count+=1
            self.progressbar['value']=(float)(count/len(self.org_selected))*100
            
            self.label_progress['text']=str(count)+"/"+str(len(self.org_selected))
            self.window.update_idletasks()
            self.window.update()


        # self.print_on_window(str(c) + " nouveau item(s) téléchargé",'green')
        self.window.update()
        self.print_on_window("--------------------Recherche terminée--------------------","white")
        self.print_on_window("","white")
        
        self.run_search.config(state=NORMAL)
        self.run_search.config(text="Recherche")
        self.button_arret.config(state=DISABLED)
        self.button_arret.config(background=HOVER_COLOR)

        self.progressbar.config(value=0)
        self.label_progress.config(text="0/0")
        self.label_current_search.config(text="Aucune recherche en cours")


if __name__ == "__main__":
    App = FENETRE()