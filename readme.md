# COPAL -- COmplexome Profile ALignment
## Synopsis

COPAL is a tool that aligns complexome profile samples using an adaptation on the dynamic time warping algorithm.
After alignment it can quantify differences in protein migration patterns between sample groups using the hausdorff effect size

## Installation

COPAL can be installed by cloning this repository.
To run the tool from python, run the maindtw.py module.
To run the tool with the GUI wrapper from python, run the gui.py module.
To run the packaged tool with GUI, execute the executable single file on your preferred platform (Windows, Linux or Mac).

### requirements

* python 3  (for python2 switch to py2_version branch)
* numpy
* pandas
* tkinter
* xlsxwriter
* xlrd

## Test Example

*test input files:*

* Gel1\_unaligned\_samples\_500\_rows.csv
* Gel2\_unaligned\_samples\_500\_rows.csv
* Gel3\_unaligned\_samples\_500\_rows.csv

To run COPAL on test data, run example_run.py
For documentation on how to run test data using the GUI, refer to the user guide

### Modules:

**COPAL main modules**

* *maindtw* -- Contains higher level functions running steps of the analysis specified below. run this module to run COPAL analysis without GUI wrapper
* *gui* -- Contains code that wraps COPAL analysis with a GUI. run this module to run COPAL with GUI wrapper

**Progressive alignment modules**

* *localdist* -- Contains functions that calculate local distances between samples. These distances are then used in dynamic time warping  
* *timewarp* -- Contains core dynamic time warping function. This function performs dynamic time warping between two samples or two alignments.
* *multipletimewarp* -- Contains functions that perform progressive alignment procedure. Pairwise alignments to determine order, and multiple alignment to get final alignment

**hausdorff scoring modules**

* *hausdorff* -- Computes the hausdorff distance between two protein migration patterns
* *findshift* -- Calculates hausdorff distances between sample pairs for all proteins using hausdorff module
* *shiftscore* -- Calculates hausdorff effect sizes between sample groups

**data preparation modules**

* *dataprep* -- Prepares data for complexome alignment, complementing the provided datasets and normalising data
* *complement_frames* -- Complements uneven datasets, generating rows for missing proteins

**output modules**

* *datatoexcel* -- Produces excel output of COPAL results
* *txtoutput* -- Produces text output with information on COPAL process. writes to file and to stdout

## License

    COPAL -- COmplexome Profile ALignment Tool
    Copyright (C) 2018  Radboud universitair medisch centrum

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
	along with this program.  If not, see <https://www.gnu.org/licenses/>.
	
## Contact

Prof. Dr. Martijn Huynen

email: Martijn.Huijnen@radboudumc.nl