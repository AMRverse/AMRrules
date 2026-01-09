**************************
Installation
**************************

Dependencies
=============

AMRrules requires Python 3.12 or higher, and needs pip to be be installed.


Download and install AMRrules
=============================

Currently, AMRrules is only available for installation via source. We recommend you set up a new conda environment while the tool is under development.

::
    conda create -n amrrules_beta -c bioconda python=3.12

    conda activate amrrules_beta

    conda install pip

    git clone https://github.com/AMRverse/AMRrules

    cd AMRrules

    make dev

After installation, you you must download the required AMRFinderPlus resource files. Run::
    
    amrrules download_resources

This will download and cache the necessary files for AMRrules to function. You only need to run this **once** after installation, or when updating resources (eg a new AMRFinderPlus database has been released).