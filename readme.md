# COPAL -- COmplexome Profile ALignment
## Synopsis

COPAL is a tool that aligns complexome profiling samples using a multi-dimensional adaptation of the dynamic time warping algorithm.
After alignment it can quantify differences in protein migration patterns between sample groups using the hausdorff effect size.

## Installation

COPAL source code can be accessed by cloning this repository.
To run the tool with python, run the copal/COPAL_main.py module.
To run the tool with the Graphical user interface (GUI) wrapper with python, run the COPAL_gui.py module.
To use the packaged tool (with GUI), download the [latest release](https://github.com/cmbi/copal/releases), and run the executable single file (Windows or Mac).

### requirements

* python 3  (for python2.7 switch to py2_version branch)
* numpy
* pandas
* xlsxwriter
* xlrd

requirements.txt: can be used to create (conda) virtual environment with all necessary requirements

example command:  conda create --name *env_name*  --file requirements.txt

## Test Example

*test data:*

* Gel1\_unaligned\_samples\_500\_rows.csv
* Gel2\_unaligned\_samples\_500\_rows.csv
* Gel3\_unaligned\_samples\_500\_rows.csv

To run COPAL on test data, run example_run.py
For documentation on how to run an analysis using the GUI, refer to the user guide

### Modules:

**COPAL main modules**

* *COPAL_main* -- Contains higher level functions running steps of the analysis specified below. run this module to run COPAL analysis without GUI wrapper
* *COPAL_gui* -- Contains code that wraps COPAL analysis with a GUI. run this module to run COPAL with GUI wrapper

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

Joeri van Strien Msc

email: Joeri.vanStrien@radboudumc.nl

## Citing COPAL

Van Strien, J., Guerrero-Castillo, S., Chatzispyrou, I. A., Houtkooper, R. H., Brandt, U., & Huynen, M. A. 
[COmplexome Profiling ALignment (COPAL) reveals remodeling of mitochondrial protein complexes in Barth syndrome.](https://doi-org.ru.idm.oclc.org/10.1093/bioinformatics/btz025)
Bioinformatics, 2019, 1–9
