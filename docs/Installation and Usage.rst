Installation and Usage
======================

As of now (December 2020), this tool works on macOS 10.15.7 (not tested yet 
on big sur) and should work on other 
\*nux operating systems. It uses Python 3.5+ as well as C++ 11. In order to 
build the C++ modules, you must have the g++ compiler. Here is the link to 
install Python_ on your device and here for the `g++`_. 

.. _Python: https://www.python.org/downloads/
.. _`g++`: https://www.cs.odu.edu/~zeil/cs250PreTest/latest/Public/installingACompiler/

Installation
^^^^^^^^^^^^

First install the repo: 

.. code-block:: bash 

    $> git clone https://github.com/zmcgrath96/hypedsearch.git


Then run the setup script that will install dependencies and build the C++ code. 

.. code-block:: bash 

    $> cd hypedsearch
    $hypedsearch> ./setup.sh

If you get a permissions error, try the following: 

.. code-block:: bash 

    $> cd hypedsearch
    $hypedsearch> chmod u+x setup.sh
    $hypedsearch> ./setup.sh

Usage 
^^^^^

There are two ways to use hypedsearch: command line arguments or the param file

command line arguments
""""""""""""""""""""""

In order to see the arguments, run 

.. code-block:: bash 

    $> python3 -m src.main --help

    usage: main.py [-h] [--spectra-folder SPECTRA_FOLDER]
               [--database-file DATABASE_FILE] [--output-dir OUTPUT_DIR]
               [--params PARAMS] [--min-peptide-len MIN_PEPTIDE_LEN]
               [--max-peptide-len MAX_PEPTIDE_LEN] [--tolerance TOLERANCE]
               [--precursor-tolerance PRECURSOR_TOLERANCE]
               [--peak-filter PEAK_FILTER]
               [--abundance-filter REL_ABUND_FILTER] [--digest DIGEST]
               [--verbose VERBOSE] [--cores CORES] [--n N]

    Tool for identifying proteins, both hybrid and non hybrid from MS/MS data

    optional arguments:
    -h, --help            show this help message and exit
    --spectra-folder SPECTRA_FOLDER
                            Path to folder containing spectra files.
    --database-file DATABASE_FILE
                            Path to .fasta file containing proteins
    --output-dir OUTPUT_DIR
                            Directory to save all figures. Default=~/
    --params PARAMS       Use the params.py file adjacent to main.py instead of
                            using command line arguments. Default=False
    --min-peptide-len MIN_PEPTIDE_LEN
                            Minimum peptide length to consider. Default=5
    --max-peptide-len MAX_PEPTIDE_LEN
                            Maximum peptide length to consider. Default=20
    --tolerance TOLERANCE
                            ppm tolerance to allow in search. Deafult=20
    --precursor-tolerance PRECURSOR_TOLERANCE
                            ppm tolerance to accept when matching precursor
                            masses. Default=10
    --peak-filter PEAK_FILTER
                            The number of peaks to take from a spectrum. The most
                            abundant peaks will be taken. Leave blank if you want
                            no filter or to use relative abundance filter.
                            Defualt=0
    --abundance-filter REL_ABUND_FILTER
                            Take only peaks from a spectrum where the abundance of
                            the peak is >= the percentage give. Leave blank if you
                            want no filter or to use peak filter. Default=0.0
    --digest DIGEST       The digest performed. Default=None
    --verbose VERBOSE     Extra printing to console during run. Default=True
    --cores CORES         The number of cores allowed to use when searching.
                            Uses at least 1 and at most the number of available
                            cores. Default=1
    --n N                 The number of alignments to keep per spectrum.
                            Default=5

Just run the :code:`python3 -m src.main` with your spectra folder, database 
file, output folder and any other parameters you want to include and hit enter.

.. note::
   
   If you choose to use command line arguments, do ot set the :code:`params` 
   flag to true. This will read from the params file. This is discussed in 
   the next section.

param file
""""""""""

If you open up the :code:`params.py` file in the :code:`src` directory, you 
will be presented with a many different variables and descriptions. Read the 
descriptions and fill in your own parameters. Once you have done this, save the 
params file. Finally, in order to run using this params file, run the following:

.. code-block:: bash 

    $hypedsearch> python3 -m src.main --params True

