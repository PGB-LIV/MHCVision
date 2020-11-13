# FDR_estimation
FDR/PEP estimation for MHC-peptide binding prediction
### **Introduction**
what is the name of this tool?

### **How to install?**
The model requires Python 3 ( >= 3.7) and the following python packages:

- pandas (>= 1.1.2)

- numpy (>= 1.19.1)

- scipy (>=1.5.2)

Note that we cannot guarantee whether MHCSeqNet will work with older versions of these packages.

To install the model:

1. Clone this repository
```
git clone https://github.com/Phorutai/FDR_estimation
```
For other methods for cloning a GitHub repository, please see  [here](https://help.github.com/articles/cloning-a-repository/)

2. Install the latest version of 'pip' and 'setuptools' packages for Python 3 if your system does not already have them
```
python -m ensurepip --default-pip
pip install setuptools
```
For more information, please see [here](https://packaging.python.org/tutorials/installing-packages/#install-pip-setuptools-and-wheel)

3.  then what??? setup.py how to do that?

### **Usage**
```
usage: noname.py [options] input_file.csv
options:
-a, --allele   REQUIRED: type the allele name i.e. HLA-A0101, which are supported in the "supplied_alleles.txt"
-i, --input    REQUIRED: specify the input filename, the input file must be in ".CSV" format (comma-separated values), the column headers must contain 'Peptide', 'IC50','%Rank'
-o, --output   specify the output filename (optional)
-h, --help     Print the usage information'
```

### **Sample scripts**
You can use input_sample.csv as the input file
```
python noname.py -a HLA-A0201 -i input_sample.csv
```
