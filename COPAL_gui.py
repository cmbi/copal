"""
Module that creates GUI of COPAL tool. Run this module to run GUI version of the COPAL tool

part of: COPAL -- COmplexome Profile ALignment Tool
Copyright (C) 2018  Radboud universitair medisch centrum
    for full notice, reference readme.md
"""

# import statements
import Tkinter as tk
import tkFileDialog
import COPAL_main
import os
import threading

# ---------- INPUT PROCESSING AND HELPER FUNCTIONS ---------- #

def get_type(type):
    """ takes file type from tk widget, transforms it to fit COPAL_main module  """
    if type == "excel":
        return "excel"
    elif type == "csv del ';' dec ','":
        return (';', ',')
    elif type == "csv del ',' dec '.'":  
        return (',', '.')
    elif type == "tsv del '\\t' dec ','":
        return ('\t', ',')
    elif type == "tsv del '\\t' dec '.'":
        return ('\t', '.')

def get_sample_columns(first, last):
    """ takes input from samplecolumn widgets, combines into list of tuples for each 
    sample """
    return [(x,y) for x,y in zip(first,last)]

def get_groups(group1, group2, samplenames):
    """takes sample names, transforms to numbers to fit main program"""
    group1_index = [samplenames.index(i) + 1 for i in group1]
    group2_index = [samplenames.index(i) + 1 for i in group2] 
    return [group1_index, group2_index]

def get_normalisation(normtype, normcol, normfile):
    """
    takes normtype and normcol input from tk wigdets, assign norm variables to fit COPAL_main

    Keyword arguments: input from tkinter widgets
    Returns:
        norm_check -- boolean, wether to perform normalisation
        normcol -- string, header of input datacolumn that contains normalisation booleans
        normfile -- string, path/name of file containing protein IDs to be used for 
                    normalisation
    """
    if normtype == "None":
        norm_check = False
        normcol = None
        normfile = None
    elif normtype == "Using all Proteins":
        norm_check = True
        normcol = None
        normfile = None
    elif normtype == "Using subset from Column":
        norm_check = True
        normfile = None
    elif normtype == "Using subset from File":
        norm_check = True
        normcol = None
        normfile = normfile
    return (norm_check, normcol, normfile)
    
def save_input_settings():
    """ saves first input file settings to global input variables """
    global input
    input['filename'] = [file_name.get()]
    input['skiprows'] = [skip_rows.get()]
    input['input_type'] = [get_type(file_type.get())]
    sample_columns_first = sample_columns_box.get("1.0", tk.END).splitlines()
    sample_columns_last = sample_columns_second_box.get("1.0", tk.END).splitlines()
    input['samplecolumns'] = [get_sample_columns(sample_columns_first, sample_columns_last)]
    input['identifier'] = [ident_col.get()]
    input['sheetname'] = [sheet_name.get()]
    input['samplenames'] = sample_names_box.get("1.0", tk.END).splitlines()

def clear_input_vars():
    """clears input variables, for use by new input frame"""
    file_name.set("")
    skip_rows.set(0)
    sample_names_box.delete('1.0', tk.END)
    sample_columns_box.delete('1.0', tk.END)
    sample_columns_second_box.delete('1.0', tk.END)
    ident_col.set("")
    sheet_name.set("")
    
def append_extra_input():
    """ appends extra file input to existing global variable lists"""
    global input
    input['filename'].append(file_name.get())
    input['skiprows'].append(skip_rows.get())
    input['input_type'].append(get_type(file_type.get()))
    sample_columns_first = sample_columns_box.get("1.0", tk.END).splitlines()
    sample_columns_last = sample_columns_second_box.get("1.0", tk.END).splitlines()
    input['samplecolumns'].append(get_sample_columns(sample_columns_first, sample_columns_last))
    input['samplenames'] += sample_names_box.get("1.0", tk.END).splitlines()
    input['identifier'].append(ident_col.get())
    input['sheetname'].append(sheet_name.get())
    
def save_output_settings():
    """saves output settings to global variables"""
    global input
    input['analysis_name'] = job_name.get()
    normalisation_parameters = get_normalisation(normalisation_type.get(), norm_col.get(), norm_file.get())
    input['norm_check'] = normalisation_parameters[0]
    input['normcol'] = normalisation_parameters[1]
    input['normfile'] = normalisation_parameters[2]
    input['hausdorff_scoring'] = score_check.get()
    input['hausd_factor'] = hausd_factor.get()
    input['gsea_output'] = gsea_check.get()
    if input['gsea_output']:
        input['GSEA_rank_column'] = gsea_col.get()
    else:
        input['GSEA_rank_column'] = None
    if input['hausdorff_scoring']:
        group1 = group1_text.get("1.0", tk.END).splitlines()
        group2 = group2_text.get("1.0", tk.END).splitlines()
        input['groups'] = get_groups(group1, group2, input['samplenames'])
    else:
        input['groups'] = None
        
def check_input():
    """
    checks input details for possible errors in input data entry by user

    Returns:
        if error found: set message in input status label and return False
        if no error: return True
    """
    if file_name.get() == '':
        input_status_var.set("no filename entered! please enter filename")
        return False
    if ident_col.get() == '':
        input_status_var.set("no identifier column entered! please enter column header")
        return False
    if file_type.get() == 'excel':
        if sheet_name.get() == '':
            input_status_var.set("no sheetname entered! please enter a sheetname when using an excel file")
            return False
    try:
        skip_rows.get()
    except:
        input_status_var.set("skip rows entry not valid. enter whole number")
        return False
    
    sample_names = sample_names_box.get("1.0", tk.END).splitlines()
    sample_columns = sample_columns_box.get("1.0", tk.END).splitlines()
    sample_columns_second = sample_columns_second_box.get("1.0", tk.END).splitlines()
    
    if sample_names == ['']:
        input_status_var.set("no sample names entered! please enter samplenames")
        return False
    if sample_columns == [''] or sample_columns_second == ['']:
        input_status_var.set("sample columns missing! please check input")
        return False
    if len(sample_names) != len(sample_columns) or len(sample_names) != len(sample_columns_second):
        input_status_var.set("number of samples not equal number of sample column pairs! please check input")
        return False
    return True

def check_output():
    """
    checks output details for possible errors by user
    Returns:
        if error: set message in status label, return false
        if no error: return True
    """
    print input['samplenames']
    if job_name.get() == '':
        status_var.set("no job name entered!")
        return False
    if normalisation_type.get() == "Using subset from Column":
        if norm_col.get() == '':
            status_var.set("no normalisation column header entered. required with this data normalisation option")
            return False
    if normalisation_type.get() == "Using subset from File":
        if norm_file.get() == '':
            status_var.set("no normalisation file entered. required with this normalisation option")
            return False
    
    try:
        hausd_factor.get()
    except:
        status_var.set("hausdorff factor entry not valid. enter numeric value")
        return False

    if score_check.get():
        group1 = group1_text.get("1.0", tk.END).splitlines()
        group2 = group2_text.get("1.0", tk.END).splitlines()
        if group1 == [''] or group2 == ['']:
            status_var.set("a group is empty. fields required when scoring option is selected")
            return False
        
        groups = group1 + group2
        for sample in groups:
            if sample not in input['samplenames']:
                status_var.set("a sample entered in one of the groups is not in sample names. please check input")
                return False
    if gsea_check.get():
        if gsea_col.get() == '':
            status_var.set("no GSEA column header entered. required field if GSEA output option is selected")
            return False
        
    return True
            
def print_input():
    """prints out global input variables to stdout"""
    print "file name: ", input['filename']
    print "file type: ", input['input_type']
    print "skiprows: ", input['skiprows']
    print "identifier: ", input['identifier']
    print "sheetname : ", input['sheetname']
    print "sample columns: ", input['samplecolumns']
    print "sample names: ", input['samplenames']

def print_output():
    """prints out global output variables to stdout"""
    print "job name: ", input['analysis_name']
    print "norm_check: ", input['norm_check']
    print "normcol: ", input['normcol']
    print "normfile: ", input['normfile']
    print "score analysis: ", input['hausdorff_scoring']
    print "hausdorff factor: ", input['hausd_factor']
    print "gsea output: ", input['gsea_output']
    print "gsea column: ", input['GSEA_rank_column']
    print "groups: ", input['groups']

def run_analysis():
    """
    runs main function in COPAL_main module, running COPAL analysis
    
    global input dictionary is provided as input argument to COPAL_main.main()
    """
    save_settings_button.unbind("<Button-1>")
    try:
        COPAL_main.main(input = input)
        status_var.set("analysis complete! your output is in: " + output_folder.get() + "/" + input['analysis_name'] + "_results")  
    except Exception as e:
        status_var.set("Error occured. Check and re-enter input.  " + str(e))
        print "something went wrong! please check and re-enter your input"
        print "error message: ", e
    save_settings_button.bind("<Button-1>", save_output_and_run_handler) 

# ---------- EVENT HANDLERS ---------- #

def back_to_input_handler(event = None):
    """clears input and goes back to first GUI frame"""
    clear_input_vars()
    first_input_frame()

def quit_app(event = None):
    """quit app button handler: quits app"""
    root.destroy()

def get_file(event = None):
    """
    choose file button handler, opens menu to browse for file
    
    assigns path/filename to global file_name variable
    """
    file_name.set(tkFileDialog.askopenfilename(title = 'Choose a file'))

def get_norm_file(event = None):
    """
    choose file button event handler, opens menu to browse for a file
    
    assigns path/gilename to global file_name variable
    """
    norm_file.set(tkFileDialog.askopenfilename(title = 'Choose a file'))

def get_output_folder(event = None):
    """
    choose folder button event handler, opens menu to browse for folder,
    
    stores path to stringvar variable output_folder
    """
    output_folder.set(tkFileDialog.askdirectory(title = 'Choose folder'))

def save_proceed_handler(event = None):
    """
    save and proceed button event handler (for first input frame)
    
    checks if input seems ok, does not proceed and sets warning if not
    """
    if check_input():
        save_input_settings()
        output_frame = tk.Frame(root)
        output_frame.grid(row = 0, column = 0, sticky = "news")
        output_options_frame(output_frame)
        output_frame.tkraise()
    
def append_proceed_handler(event = None):
    """
    append and proceed button event handler (for later input frames)
    
    appends input, creates output frame, adds widgets to output frame
    checks if input seems okm does not proceed and sets warning if not
    """
    if check_input():
        append_extra_input()
        print_input()
        output_frame = tk.Frame(root)
        output_frame.grid(row = 0, column = 0, sticky = "news")
        output_options_frame(output_frame)
        output_frame.tkraise()
    
def save_extra_input_handler(event = None):
    """
    save and extra input button event handler (for first input frame)
    
    saves input to globals, clears input variables, creates new input frame
    checks if input seems ok, waits and sets warning otherwise
    """
    if check_input():   
        save_input_settings()       
        clear_input_vars()     
        extra_frame = tk.Frame(root)
        extra_frame.grid(row = 0, column = 0, sticky = "news")
        input_frame(extra_frame)
        extra_input_buttons(extra_frame)
        extra_frame.tkraise()

def append_extra_input_handler(event = None):
    """
    append_extra input button event handler (for later input frames)
    
    appends extra input to globals, clears input vars, creates new input frame
    checks if input seems ok, waits and sets warning otherwise
    """
    if check_input():
        append_extra_input()
        clear_input_vars()
        extra_frame = tk.Frame(root)
        extra_frame.grid(row = 0, column = 0, sticky = "news")
        input_frame(extra_frame)
        extra_input_buttons(extra_frame)
        extra_frame.tkraise()

def save_output_and_run_handler(event = None):
    """
    save output options and run button event handler (for output options frame)
    
    saves output to globals, set status label stringvar, tries main
    in case of except: set status label with error message
    """
    if check_output():
        try:
            os.chdir(output_folder.get())
        except:
            status_var.set("Output folder does not exist! re-select folder and try again")
            root.update_idletasks()
            return
        
        status_var.set("Running analysis. This might take a while...")
        root.update_idletasks()                # ensures the status label gets updated in a timely fashion
        save_output_settings()
        print_input()
        print_output()
        analysis_thread = threading.Thread(target = run_analysis)
        analysis_thread.daemon = True
        analysis_thread.start()

# ---------- GUI WIDGET CREATING FUNCTIONS ---------- #
def input_frame(master):
    """
    creates input widgets for file input frames, assigns text widget content to globals
    """
    global sample_names_box, sample_columns_box, sample_columns_second_box
    tk.Label(master, text = "COPAL -- COmplexome Profile ALignment\n", font = ("Helvetica", 18), pady = 20).grid(columnspan = 6)
    
    # filename label and entry
    filename_label = tk.Label(master, text = "filename:").grid(column = 1, sticky = tk.E)
    filename_entry = tk.Entry(master, textvariable = file_name, width = 50)         # create entry in master, set variable to strVar
    filename_entry.grid(row = 1, column = 2, columnspan = 2, sticky = tk.W, pady = 4)                            # add to frame with pack 

    # choose filename button
    get_file_button = tk.Button(master, text = "Choose file")
    get_file_button.bind("<Button-1>", get_file)
    get_file_button.grid(row = 1, column = 4, padx = 4, sticky = tk.W)

    # choose filetype dropdown menu
    file_type_dropdown = tk.OptionMenu(master, file_type,"excel", "csv del ';' dec ','",
                                         "csv del ',' dec '.'", "tsv del '\\t' dec ','",
                                         "tsv del '\\t' dec '.'")
    file_type_dropdown.grid(row = 1, column = 5, sticky = tk.EW)

    # skip rows label and entry
    skip_rows_label = tk.Label(master, text = "skip rows:").grid(sticky = tk.E)
    skip_rows_entry = tk.Entry(master, textvariable = skip_rows, width = 3).grid(row = 2, column = 1)

    #identifier label and entry
    identifier_label = tk.Label(master, text = "protein identifier column:").grid(row = 2, column = 2, sticky = tk.E)
    identifier_entry = tk.Entry(master, textvariable = ident_col, width = 25).grid(row = 2, column = 3, sticky = tk.W)

    # sheetname label and entry
    sheetname_label = tk.Label(master, text = "sheetname:").grid(row = 2, column = 4, sticky = tk.E)
    sheetname_entry = tk.Entry(master, textvariable = sheet_name, width = 25).grid(row = 2, column = 5, padx = 4)       

    # sample names label and text box
    sample_names_label = tk.Label(master, text = "sample names").grid(row = 3, column = 0, columnspan = 2, pady = 5)
    sample_names_box = tk.Text(master, height = 15, width = 30, background = "gray", wrap = tk.NONE)
    sample_names_box.grid(row = 4, column = 0, columnspan = 2, padx = 10)

    # sample columns labels and text boxes
    sample_columns_label = tk.Label(master, text = "first column").grid(row = 3, column = 2)
    sample_columns_box = tk.Text(master, height = 15, width = 30, background = "gray", wrap = tk.NONE)
    sample_columns_box.grid(row = 4, column = 2)
    sample_columns_second_label = tk.Label(master, text = "last column").grid(row = 3, column = 3)
    sample_columns_second_box = tk.Text(master, height = 15, width = 30, background = "gray", wrap = tk.NONE)
    sample_columns_second_box.grid(row = 4, column = 3)
    
    # quit application button
    quit_app_button = tk.Button(master, text = "quit     ")
    quit_app_button.bind("<Button-1>", quit_app)
    quit_app_button.grid(row = 5, column = 0, sticky = tk.W)    
    
    # input status bar
    status_label = tk.Label(master, textvariable = input_status_var).grid(row = 5, column = 1, columnspan = 3)

    
def first_input_buttons(master):
    """adds button widgets, for first input frame"""

    #add another input file button
    new_frame_button = tk.Button(master, text = "add another file")
    new_frame_button.bind("<Button-1>", save_extra_input_handler)
    new_frame_button.grid(row = 5, column = 4, sticky = tk.E)
    
    # save and proceed to output frame
    proceed_button = tk.Button(master, text = "proceed to output")
    proceed_button.bind("<Button-1>", save_proceed_handler)
    proceed_button.grid(row = 5, column = 5, sticky = tk.E)

def extra_input_buttons(master):
    """
    adds button widgets, for extra input frames 
    
    add another input file button
    """
    new_frame_button = tk.Button(master, text = "add another file")
    new_frame_button.bind("<Button-1>", append_extra_input_handler)
    new_frame_button.grid(row = 5, column = 4, sticky = tk.E)
    
    # append input and proceed button
    append_input_button = tk.Button(master, text = "proceed to output")
    append_input_button.bind("<Button-1>", append_proceed_handler)
    append_input_button.grid(row = 5, column = 5, sticky = tk.E)

def output_options_frame(master):
    """
    creates output widgets for output info frame
    
    asigns text box widgets' content to global variables
    """
    global group1_text, group2_text, save_settings_button
    # top label
    tk.Label(master, text = "COPAL -- COmplexome Profile ALignment\n", font = ("Helvetica", 18), pady = 20).grid(columnspan = 6)

    # job_name label and entry
    job_name_label = tk.Label(master, text = "job name:").grid(column = 0, sticky = tk.E)
    job_name_entry = tk.Entry(master, textvariable = job_name, width = 50)         # create entry in master, set variable to strVar
    job_name_entry.grid(row = 1, column = 1, columnspan = 2, sticky = tk.W, pady = 4)                            # add to frame with pack 

    # select output location folder
    select_folder_label = tk.Label(master, text = "select output folder:").grid(column = 0, sticky = tk.E)
    select_folder_entry = tk.Entry(master, textvariable = output_folder, width = 25)
    select_folder_entry.grid(row = 2, column = 1, sticky = tk.W, pady = 4)
    
    # choose folder button
    choose_folder_button = tk.Button(master, text = "Choose folder")
    choose_folder_button.bind("<Button-1>", get_output_folder)
    choose_folder_button.grid(row = 2, column = 2, sticky = tk.W)
    
    # data normalisation drop-down
    normalisation_label = tk.Label(master, text = "data normalisation:").grid(sticky = tk.E, pady = 5)
    normalisation_dropdown = tk.OptionMenu(master, normalisation_type, "None", "Using all Proteins", "Using subset from Column", "Using subset from File")
    normalisation_dropdown.grid(row = 3, column = 1, columnspan = 2, sticky = tk.EW)
    
    # norm_col label and entry
    norm_col_label = tk.Label(master, text = "if from column:").grid(sticky = tk.E)
    norm_col_entry = tk.Entry(master, textvariable = norm_col, width = 25).grid(row = 4, column = 1)

    # norm_file label and entry
    norm_file_label = tk.Label(master, text = "if from file:").grid(row = 5, sticky = tk.E)
    norm_file_entry = tk.Entry(master, textvariable = norm_file, width = 25).grid(row = 5, column = 1, sticky = tk.W)
    
    # choose norm_file button
    norm_file_button = tk.Button(master, text = "Choose file")
    norm_file_button.bind("<Button-1>", get_norm_file)
    norm_file_button.grid(row = 5, column = 2, padx = 4, sticky = tk.W)

    # score analysis check button
    score_analysis_check = tk.Checkbutton(master, text = "perform score analysis", variable = score_check)
    score_analysis_check.grid(sticky = tk.W)        # create checkbutton, set text, assign to boolVar variable

    # score analysis group text boxes
    group1_label = tk.Label(master, text = "group 1 samples").grid(row = 6, column = 1)
    group2_label = tk.Label(master, text = "group 2 samples").grid(row = 6, column = 2)
    group1_text = tk.Text(master, height = 10, width = 25, background = "gray", wrap = tk.NONE)
    group1_text.grid(column = 1)
    group2_text = tk.Text(master, height = 10, width = 25, background = "gray", wrap = tk.NONE)
    group2_text.grid(row = 7, column = 2)
    
    # hausdorff factor label and entry
    hausd_factor_label = tk.Label(master, text = 'hausdorff factor:').grid(row = 6, column = 3)
    hausd_factor_entry = tk.Entry(master, textvariable = hausd_factor, width = 10).grid(row = 7, column = 3, sticky = tk.N)

    # GSEA check button and column entry
    GSEA_check = tk.Checkbutton(master, text = "provide rank ordered protein list", variable = gsea_check).grid(columnspan = 2, sticky = tk.W, pady = 10)
    GSEA_col_label = tk.Label(master, text = "ranked list identifier column:").grid(row = 8, column = 2, sticky = tk.E)
    GSEA_col_entry = tk.Entry(master, textvariable = gsea_col, width = 25).grid(row = 8, column = 3, sticky = tk.W) 

    # status label
    status_label = tk.Label(master, textvariable = status_var).grid(row = 9, column = 1, columnspan = 4)
        
    # save and run button
    save_settings_button = tk.Button(master, text = "save and run")
    save_settings_button.bind("<Button-1>", save_output_and_run_handler)
    save_settings_button.grid(row = 9, column = 5, sticky = tk.E)
    
    # quit application button
    quit_app_button = tk.Button(master, text = "quit  ")
    quit_app_button.bind("<Button-1>", quit_app)
    quit_app_button.grid(row = 9, column = 0, sticky = tk.W)    
    
    # back to input button
    back_to_input_button = tk.Button(master, text = "back to start")
    back_to_input_button.bind("<Button-1>", back_to_input_handler)
    back_to_input_button.grid(row = 9, column = 0, sticky = tk.E)

    # sample names message box
    sample_names_message = "Sample names:\n\n" + '\n'.join(input['samplenames'])
    names_message_box = tk.Message(master, text = sample_names_message, width = 250)
    names_message_box.grid(row = 1, column = 4, rowspan = 6, sticky = tk.NE, padx = 20)


# FIRST INPUT FRAME CREATER
def first_input_frame():
    """initialises and creates first frame of GUI. later frames are created on the fly"""
    first_frame = tk.Frame(root)
    first_frame.grid(row = 0, column = 0, sticky = "news")
    input_frame(first_frame)
    first_input_buttons(first_frame)

if __name__ == "__main__":
    # initialise root
    root = tk.Tk()
    root.wm_title("master")
    root.resizable(0,0)

    # initialize tkinter variables
    job_name = tk.StringVar()
    file_name = tk.StringVar()
    skip_rows = tk.IntVar()
    file_type = tk.StringVar()
    ident_col = tk.StringVar()
    sheet_name = tk.StringVar()
    normalisation_type = tk.StringVar()
    norm_col = tk.StringVar()
    norm_file = tk.StringVar()
    score_check = tk.BooleanVar()
    hausd_factor = tk.DoubleVar()
    hausd_factor.set(1.0)
    gsea_check = tk.BooleanVar()
    gsea_col = tk.StringVar()
    status_var = tk.StringVar()
    input_status_var = tk.StringVar()
    output_folder = tk.StringVar()

    # initialize global dictionary containing input for analysis
    input = {}

    # set variable default values
    file_type.set("excel")
    normalisation_type.set("None")
    score_check.set(True)
    gsea_check.set(True)
    status_var.set("press save and run to start analysis")
    input_status_var.set("add another file or proceed to ouput details")

    # open first frame, start mainloop
    first_input_frame()
    root.mainloop()