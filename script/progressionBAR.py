def init_bar(progressbar, root):
    progressbar['value'] = 0
    root.update_idletasks()
    
def update_bar(progressbar,root,n_line,i):
    if i%(max(1,n_line//100)) == 0:
        progressbar['value'] = i/n_line*100
        root.update_idletasks()
        root.update()