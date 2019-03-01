"""
Module that creates GUI of COPAL tool. Run this module to run GUI version of the COPAL tool

part of: COPAL -- COmplexome Profile ALignment Tool
Copyright (C) 2018  Radboud universitair medisch centrum
    for full notice, reference readme.md
"""

# import external packages
import tkinter as tk
from tkinter import filedialog
from tkinter import ttk
from tkinter import messagebox
import os
import threading
from pandas.errors import ParserError

#import copal
import copal
from copal.dataprep import multi_dataload

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

def clear_last_input():
    """
    removes last added set of input specs from global input dict 
    """
    added_samples = len(input['samplecolumns'][-1])
    print('clearing wrong input..')
    print('number of added samples: {}'.format(added_samples))
    del input['identifier'][-1]
    del input['filename'][-1]
    del input['sheetname'][-1]
    del input['skiprows'][-1]
    del input['input_type'][-1]
    del input['samplecolumns'][-1]
    del input['samplenames'][-1*added_samples:]

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
    input['warp_method'] = warp_method.get()
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
        warn_msg = "no filename entered! please enter filename"
        input_status_var.set(warn_msg)
        messagebox.showwarning('input error',warn_msg)
        return False

    if ident_col.get() == '':
        warn_msg = "no identifier column entered! please enter column header"
        input_status_var.set(warn_msg)
        messagebox.showwarning('input error',warn_msg)
        return False

    if file_type.get() == 'excel':
        if sheet_name.get() == '':
            warn_msg = "no sheetname entered! please enter a sheetname when using an excel file"
            input_status_var.set(warn_msg)
            messagebox.showwarning('input error',warn_msg)
            return False
    try:
        skip_rows.get()
    except:
        warn_msg = "skip rows entry not valid. enter whole number"
        input_status_var.set(warn_msg)
        messagebox.showwarning("input error",warn_msg)
        return False

    sample_names = sample_names_box.get("1.0", tk.END).splitlines()
    sample_columns = sample_columns_box.get("1.0", tk.END).splitlines()
    sample_columns_second = sample_columns_second_box.get("1.0", tk.END).splitlines()

    if sample_names == ['']:
        warn_msg ="no sample names entered! please enter samplenames"
        input_status_var.set(warn_msg)
        messagebox.showwarning("input error",warn_msg)
        return False

    if sample_columns == [''] or sample_columns_second == ['']:
        warn_msg = "sample columns missing! please check input"
        input_status_var.set(warn_msg)
        messagebox.showwarning("input error",warn_msg)
        return False

    if len(sample_names) != len(sample_columns) or len(sample_names) != len(sample_columns_second):
        warn_msg = "number of samples not equal number of sample column pairs! please check input"
        input_status_var.set(warn_msg)
        messagebox.showwarning("input error",warn_msg)
        return False

    return True

def check_input_data():
    """
    check if input data is valid and matches input specs
    """
    identifier = input['identifier'][-1]
    filename = input['filename'][-1]
    sheet = input['sheetname'][-1]
    skip_rows = input['skiprows'][-1]
    file_type = input['input_type'][-1]
    samplecolumns = input['samplecolumns'][-1]
    
    # try to load data
    try:
        data = multi_dataload(identifier, filename, sheet, skip_rows, file_type)
    except FileNotFoundError as error:
        warn_msg = "file not found! check input file name and path"
        input_status_var.set(warn_msg)
        messagebox.showwarning("input error", warn_msg)
        clear_last_input()
        return False

    except ParserError as error:
        warn_msg = "problem parsing input file. is file type correct? error: {}".format(error)
        input_status_var.set('loading data failed. check file type')
        messagebox.showwarning("input error", warn_msg)
        clear_last_input()
        return False

    except ValueError as error:
        warn_msg = "something went wrong with loading data file. please check file input details. error: {}".format(error)
        input_status_var.set("loading data failed. please check input details.")
        messagebox.showwarning("input error", warn_msg)
        clear_last_input()
        return False

    except Exception as error:
        warn_msg = "something went wrong with loading data. check input details. error: {}".format(error)
        input_status_var.set("problem loading data. check input details")
        messagebox.showwarning("input error",warn_msg)
        clear_last_input()
        return False

    # check if samplecolumns exist
    dataset_columns = list(data.columns)
    for col_pair in samplecolumns:
        if not col_pair[0] in dataset_columns:
            warn_msg = "start column header not in dataset: {}. check input".format(col_pair[0])
            input_status_var.set(warn_msg)
            messagebox.showwarning("input error",warn_msg)
            clear_last_input()
            return False
        if not col_pair[1] in dataset_columns:
            warn_msg = "end column header not in dataset: {}. check input".format(col_pair[1])
            input_status_var.set(warn_msg)
            messagebox.showwarning("input error",warn_msg)
            clear_last_input()
            return False

        # check for samples with 0 or 1 slices
        try:
            sample_data = data.loc[:,col_pair[0]:col_pair[1]]
        except Exception as error:
            warn_msg = "sample could not be extracted from data. error: {}".format(error)
            input_status_var.set("problem extracting sample from data. check input")
            messagebox.showwarning("input error",warn_msg)
            clear_last_input()
            return False            

        if sample_data.shape[1] < 2:
            warn_msg = "found sample with 0 or 1 columns: {}:{}. ".format(col_pair[0],col_pair[1])
            input_status_var.set(warn_msg)
            messagebox.showwarning("input error", warn_msg)
            clear_last_input()
            return False

        # check for nan values in sample
        if sample_data.isnull().values.any():
            warn_msg = "sample contains missing (nan) values: {}:{}. ".format(col_pair[0],col_pair[1])
            input_status_var.set(warn_msg)
            messagebox.showwarning("input error", warn_msg)
            clear_last_input()
            return False
    
    # store columns of each dataset to check outputspecs
    input['dataset_columns'].append(dataset_columns)
    return True 

def check_output():
    """
    checks output details for possible errors by user
    Returns:
        if error: set message in status label, return false
        if no error: return True
    """
    if job_name.get() == '':
        warn_msg = "no job name entered!"
        status_var.set(warn_msg)
        messagebox.showwarning("input error",warn_msg)
        return False
    if normalisation_type.get() == "Using subset from Column":
        if norm_col.get() == '':
            warn_msg = "no normalisation column header entered. required with this data normalisation option"
            status_var.set(warn_msg)
            messagebox.showwarning("input error",warn_msg)
            return False
    if normalisation_type.get() == "Using subset from File":
        if norm_file.get() == '':
            warn_msg = "no normalisation file entered. required with this normalisation option"
            status_var.set(warn_msg)
            messagebox.showwarning("input error",warn_msg)
            return False

    try:
        hausd_factor.get()
    except:
        warn_msg = "hausdorff factor entry not valid. enter numeric value"
        status_var.set(warn_msg)
        messagebox.showwarning("input error",warn_msg)
        return False

    if score_check.get():
        group1 = group1_text.get("1.0", tk.END).splitlines()
        group2 = group2_text.get("1.0", tk.END).splitlines()
        if group1 == [''] or group2 == ['']:
            warn_msg = "a group is empty. fields required when scoring option is selected"
            status_var.set(warn_msg)
            messagebox.showwarning("input error",warn_msg)
            return False

        groups = group1 + group2
        for sample in groups:
            if sample not in input['samplenames']:
                warn_msg = "a sample entered in one of the groups is not in sample names. please check input"
                status_var.set(warn_msg)
                messagebox.showwarning("input error",warn_msg)
                return False
    if gsea_check.get():
        if gsea_col.get() == '':
            warn_msg = "no GSEA column header entered. required field if GSEA output option is selected"
            status_var.set(warn_msg)
            messagebox.showwarning("input error",warn_msg)
            return False

    return True

def check_output_data():
    """
    checks if output data is valid
    """
    input_columnsets = input['dataset_columns']
    normcol = input['normcol']
    normfile = input['normfile']
    gseacol = input['GSEA_rank_column']

    # check normcol 
    if normcol:
        for dataset in input_columnsets:
            if not normcol in dataset:
                warn_msg = "normalisation column header not present in one of the input datasets!"
                status_var.set(warn_msg)
                messagebox.showwarning("input error",warn_msg)
                return False
    # check gsea_col
    if gseacol:
        for dataset in input_columnsets:
            if not gseacol in dataset:
                warn_msg = "gsea column header not present in one of the input datasets!"
                status_var.set(warn_msg)
                messagebox.showwarning("input error",warn_msg)
                return False
    
    # check normfile
    if normfile:
        if not os.path.isfile(normfile):
            warn_msg = "normfile does not exist! check input."
            status_var.set(warn_msg)
            messagebox.showwarning("input error",warn_msg)
            return False
    
    return True

def print_input():
    """prints out global input variables to stdout"""
    print("file name: ", input['filename'])
    print("file type: ", input['input_type'])
    print("skiprows: ", input['skiprows'])
    print("identifier: ", input['identifier'])
    print("sheetname : ", input['sheetname'])
    print("sample columns: ", input['samplecolumns'])
    print("sample names: ", input['samplenames'])

def print_output():
    """prints out global output variables to stdout"""
    print("job name: ", input['analysis_name'])
    print("norm_check: ", input['norm_check'])
    print("normcol: ", input['normcol'])
    print("normfile: ", input['normfile'])
    print("score analysis: ", input['hausdorff_scoring'])
    print("hausdorff factor: ", input['hausd_factor'])
    print("gsea output: ", input['gsea_output'])
    print("gsea column: ", input['GSEA_rank_column'])
    print("groups: ", input['groups'])

def run_analysis():
    """
    runs main function in COPAL_main module, running COPAL analysis

    global input dictionary is provided as input argument to COPAL_main.main()
    """
    save_settings_button.unbind("<Button-1>")
    try:
        copal.main(input = input)
        status_var.set("analysis complete! your output is in: " + output_folder.get() + "/" + input['analysis_name'] + "_results")
    except Exception as e:
        status_var.set("Error occured. Check and re-enter input.  " + str(e))
        print("something went wrong! please check and re-enter your input")
        print("error message: ", e)
    save_settings_button.bind("<Button-1>", save_output_and_run_handler)

# ---------- EVENT HANDLERS ---------- #

def back_to_input_handler(event = None):
    """clears input and goes back to first GUI frame"""
    # clear dataset columns when re-entering input
    input['dataset_columns'] = []
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
    file_name.set(filedialog.askopenfilename(title = 'Choose a file'))

def get_norm_file(event = None):
    """
    choose file button event handler, opens menu to browse for a file

    assigns path/gilename to global file_name variable
    """
    norm_file.set(filedialog.askopenfilename(title = 'Choose a file'))

def get_output_folder(event = None):
    """
    choose folder button event handler, opens menu to browse for folder,

    stores path to stringvar variable output_folder
    """
    output_folder.set(filedialog.askdirectory(title = 'Choose folder'))

def save_proceed_handler(event = None):
    """
    save and proceed button event handler (for first input frame)

    checks if input seems ok, does not proceed and sets warning if not
    """
    if check_input():
        save_input_settings()
        if check_input_data():
            if len(input['samplenames']) < 2:
                warn_msg = ('only 1 sample specified! need at least two samples for alignment. '
                                     'add another sample from this file or add another file')
                status_var.set(warn_msg)
                messagebox.showwarning("input error",warn_msg)
                return 
            print_input()
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
        if check_input_data():
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
        if check_input_data():
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
        if check_input_data():
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

        save_output_settings()
        if check_output_data():
            status_var.set("Running analysis. This might take a while...")
            root.update_idletasks()                # ensures the status label gets updated in a timely fashion
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
    filename_entry = ttk.Entry(master, textvariable = file_name, width = 50)         # create entry in master, set variable to strVar
    filename_entry.grid(row = 1, column = 2, columnspan = 2, sticky = tk.W, pady = 4)                            # add to frame with pack

    # choose filename button
    get_file_button = ttk.Button(master, text = "Choose file")
    get_file_button.bind("<Button-1>", get_file)
    get_file_button.grid(row = 1, column = 4, padx = 4, sticky = tk.W)

    # choose filetype dropdown menu
    file_type_dropdown = ttk.OptionMenu(master, file_type,"csv del ',' dec '.'","excel", "csv del ';' dec ','",
                                         "csv del ',' dec '.'", "tsv del '\\t' dec ','",
                                         "tsv del '\\t' dec '.'")
    file_type_dropdown.grid(row = 1, column = 5, sticky = tk.EW)

    # skip rows label and entry
    skip_rows_label = tk.Label(master, text = "skip rows:").grid(sticky = tk.E)
    skip_rows_entry = ttk.Entry(master, textvariable = skip_rows, width = 3).grid(row = 2, column = 1)

    #identifier label and entry
    identifier_label = tk.Label(master, text = "protein identifier column:").grid(row = 2, column = 2, sticky = tk.E)
    identifier_entry = ttk.Entry(master, textvariable = ident_col, width = 25).grid(row = 2, column = 3, sticky = tk.W)

    # sheetname label and entry
    sheetname_label = tk.Label(master, text = "sheetname:").grid(row = 2, column = 4, sticky = tk.E)
    sheetname_entry = ttk.Entry(master, textvariable = sheet_name, width = 25).grid(row = 2, column = 5, padx = 4)

    # sample names label and text box
    sample_names_label = tk.Label(master, text = "sample names").grid(row = 3, column = 0, columnspan = 2, pady = 5)
    sample_names_box = tk.Text(master, height = 15, width = 30, background = "white", wrap = tk.NONE,cursor='arrow')
    sample_names_box.grid(row = 4, column = 0, columnspan = 2, padx = 10)

    # sample columns labels and text boxes
    sample_columns_label = tk.Label(master, text = "first column").grid(row = 3, column = 2)
    sample_columns_box = tk.Text(master, height = 15, width = 30, background = "white", wrap = tk.NONE, cursor='arrow')
    sample_columns_box.grid(row = 4, column = 2)
    sample_columns_second_label = tk.Label(master, text = "last column").grid(row = 3, column = 3)
    sample_columns_second_box = tk.Text(master, height = 15, width = 30, background = "white", wrap = tk.NONE, cursor='arrow')
    sample_columns_second_box.grid(row = 4, column = 3)

    # quit application button
    quit_app_button = ttk.Button(master, text = "quit     ")
    quit_app_button.bind("<Button-1>", quit_app)
    quit_app_button.grid(row = 5, column = 0, sticky = tk.W,pady=5,padx=5)

    # input status bar
    status_label = tk.Label(master, textvariable = input_status_var).grid(row = 5, column = 1, columnspan = 3)


def first_input_buttons(master):
    """adds button widgets, for first input frame"""

    #add another input file button
    new_frame_button = ttk.Button(master, text = "add another file")
    new_frame_button.bind("<Button-1>", save_extra_input_handler)
    new_frame_button.grid(row = 5, column = 4, sticky = tk.E)

    # save and proceed to output frame
    proceed_button = ttk.Button(master, text = "proceed to output")
    proceed_button.bind("<Button-1>", save_proceed_handler)
    proceed_button.grid(row = 5, column = 5, sticky = tk.E)

def extra_input_buttons(master):
    """
    adds button widgets, for extra input frames

    add another input file button
    """
    new_frame_button = ttk.Button(master, text = "add another file")
    new_frame_button.bind("<Button-1>", append_extra_input_handler)
    new_frame_button.grid(row = 5, column = 4, sticky = tk.E)

    # append input and proceed button
    append_input_button = ttk.Button(master, text = "proceed to output")
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
    job_name_entry = ttk.Entry(master, textvariable = job_name, width = 50)         # create entry in master, set variable to strVar
    job_name_entry.grid(row = 1, column = 1, columnspan = 2, sticky = tk.W, pady = 4)                            # add to frame with pack

    # select output location folder
    select_folder_label = tk.Label(master, text = "select output folder:").grid(column = 0, sticky = tk.E)
    select_folder_entry = ttk.Entry(master, textvariable = output_folder, width = 25)
    select_folder_entry.grid(row = 2, column = 1, sticky = tk.W, pady = 4)

    # choose folder button
    choose_folder_button = ttk.Button(master, text = "Choose folder")
    choose_folder_button.bind("<Button-1>", get_output_folder)
    choose_folder_button.grid(row = 2, column = 2, sticky = tk.W)

    # choose warping type drop-down
    warptype_label = tk.Label(master, text = 'warping type:').grid(sticky = tk.E, pady = 5)
    warptype_dropdown = ttk.OptionMenu(master, warp_method, "interpolate","interpolate","repeat")
    warptype_dropdown.grid(row = 3, column = 1, columnspan = 2, sticky = tk.EW)

    # data normalisation drop-down
    normalisation_label = tk.Label(master, text = "data normalisation:").grid(sticky = tk.E, pady = 5)
    normalisation_dropdown = ttk.OptionMenu(master, normalisation_type, "None","None", "Using all Proteins", "Using subset from Column", "Using subset from File")
    normalisation_dropdown.grid(row = 4, column = 1, columnspan = 2, sticky = tk.EW)

    # norm_col label and entry
    norm_col_label = tk.Label(master, text = "if from column:").grid(sticky = tk.E)
    norm_col_entry = ttk.Entry(master, textvariable = norm_col, width = 25).grid(row = 5, column = 1,sticky=tk.W)

    # norm_file label and entry
    norm_file_label = tk.Label(master, text = "if from file:").grid(row = 6, sticky = tk.E)
    norm_file_entry = ttk.Entry(master, textvariable = norm_file, width = 25).grid(row = 6, column = 1, sticky = tk.W)

    # choose norm_file button
    norm_file_button = ttk.Button(master, text = "Choose file")
    norm_file_button.bind("<Button-1>", get_norm_file)
    norm_file_button.grid(row = 6, column = 2, padx = 4, sticky = tk.W)

    # score analysis check button
    score_analysis_check = ttk.Checkbutton(master, text = "perform score analysis", variable = score_check)
    score_analysis_check.grid(sticky = tk.W)        # create checkbutton, set text, assign to boolVar variable

    # score analysis group text boxes
    group1_label = tk.Label(master, text = "group 1 samples").grid(row = 7, column = 1)
    group2_label = tk.Label(master, text = "group 2 samples").grid(row = 7, column = 2)
    group1_text = tk.Text(master, height = 10, width = 25, background = "white", wrap = tk.NONE,cursor = 'arrow')
    group1_text.grid(column = 1)
    group2_text = tk.Text(master, height = 10, width = 25, background = "white", wrap = tk.NONE,cursor='arrow')
    group2_text.grid(row = 8, column = 2)

    # hausdorff factor label and entry
    hausd_factor_label = tk.Label(master, text = 'hausdorff factor:').grid(row = 7, column = 3)
    hausd_factor_entry = ttk.Entry(master, textvariable = hausd_factor, width = 10).grid(row = 8, column = 3, sticky = tk.N)

    # GSEA check button and column entry
    GSEA_check = ttk.Checkbutton(master, text = "provide rank ordered protein list", variable = gsea_check).grid(columnspan = 2, sticky = tk.W, pady = 10)
    GSEA_col_label = tk.Label(master, text = "ranked list identifier column:").grid(row = 9, column = 2, sticky = tk.E)
    GSEA_col_entry = ttk.Entry(master, textvariable = gsea_col, width = 25).grid(row = 9, column = 3, sticky = tk.W)

    # status label
    status_label = tk.Label(master, textvariable = status_var).grid(row = 10, column = 1, columnspan = 4)

    # save and run button
    save_settings_button = ttk.Button(master, text = "save and run")
    save_settings_button.bind("<Button-1>", save_output_and_run_handler)
    save_settings_button.grid(row = 10, column = 5, sticky = tk.E)

    # quit application button
    quit_app_button = ttk.Button(master, text = "quit  ")
    quit_app_button.bind("<Button-1>", quit_app)
    quit_app_button.grid(row = 10, column = 0, sticky = tk.W,pady=5,padx=5)

    # back to input button
    back_to_input_button = ttk.Button(master, text = "back to start")
    back_to_input_button.bind("<Button-1>", back_to_input_handler)
    back_to_input_button.grid(row = 10, column = 1, sticky = tk.W,pady=5,padx=5)

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

    print('Loading COPAL..')

    # initialise root
    root = tk.Tk()
    root.wm_title("COPAL")
    root.geometry("1080x520")
    #root.resizable(0,0)        # prevents window from being resized

    # try to set application icon (does not work on all platforms)
    try:
        root.iconbitmap(os.path.join(os.getcwd(),'copal/static/Copal.ico'))
    except Exception as e: 
        print('failed to load copal icon')
        print(e)        

    # initialize tkinter variables  
    job_name = tk.StringVar()
    file_name = tk.StringVar()

    skip_rows = tk.IntVar()
    file_type = tk.StringVar()
    ident_col = tk.StringVar()
    sheet_name = tk.StringVar()
    normalisation_type = tk.StringVar()
    warp_method = tk.StringVar()
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
    # create empty list. to enter column-list of each dataset into.
    # used to check if output frame specs match input data
    input['dataset_columns'] = []

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
