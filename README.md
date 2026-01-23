# Overview

[Catfish](https://github.com/Kingsford-Group/catfish) is an efficient heuristic algorithm for decomposing a given flow into a set of minimum number of paths.
**Catfish-LP** incorporates a lightweight linear programming (LP) model to address the core limitations of catfish and largely improves the decomposition quality.
Please refer to [our manuscript](https://www.biorxiv.org/content/10.64898/2025.12.11.693570v1) for more details of the model.


# Installation
Download catfish-LP:
```
git clone https://github.com/Shao-Group/catfish-LP.git
```
Catfish requires Boost and catfish-LP requires Gurobi.

## Install Boost
Boost can be downloaded [here](https://www.boost.org/releases/latest/).
Decompress it and note down the directory (compiling and installing are not necessary).

## Install Gurobi
Gurobi Optimizer can be downloaded [here](https://www.gurobi.com/downloads/gurobi-software/); instructions for installation available [here](https://support.gurobi.com/hc/en-us/articles/14799677517585-Getting-Started-with-Gurobi-Optimizer).
Once completed, note down the library version at `$GUROBI_HOME/lib/libgurobiXXX.so` and modify [Makefile.am in src](src/Makefile.am) by replacing `-lgurobi110` with the installed version.


## Compile Catfish-LP
Use the following to compile and install catfish-LP:
```
aclocal && autoreconf --install
./configure --prefix=/path/to/your/installation/directory --with-boost=/path/to/your/boost/directory
make -j
make install
```
The `--prefix` argument for `configure` specifies the directory where you would put the binary of catfish-LP.
It is optional and the default is `/usr/local/bin` for most Linux distributions.
If Boost has been installed in your system, the `--with-boost` argument for `configure` can also be omitted.


# Usage
Catfish-LP accepts the same command line arguments as [catfish](https://github.com/Kingsford-Group/catfish):
```
catfish -a full -i <input-file> -o <output-file> [-h] [-v]
```

# Reproduce Results in the Manuscript
The benchmark datasets with perfect splice graphs introduced by catfish can be found at [catfishtest](https://github.com/Kingsford-Group/catfishtest).
The tool `catcompare` is used to compare decomposition results with the ground truths.

Additional simulated graphs with various complexity can be found in [data](data/).


In comparisons, the results of *greedy-width* are obtained using the original [catfish](https://github.com/Kingsford-Group/catfish):
```
catfish -a greedy -i <input-file> -o <output-file>
```
The results of *catfish* are obtained using the original [catfish](https://github.com/Kingsford-Group/catfish):
```
catfish -a full -i <input-file> -o <output-file>
```
The results of *catfish-LP* are obtained with this repository:
```
catfish -a full -i <input-file> -o <output-file>
```
The results of *ILP* are obtained using the [accelerated ILP model](https://drops.dagstuhl.de/storage/00lipics/lipics-vol301-sea2024/LIPIcs.SEA.2024.14/LIPIcs.SEA.2024.14.pdf):
```
python mfd_optimization.py -i <input-file> -o <output-file> --heuristic <greedy-solution>
```
where the required greedy solution can be transformed from the output file of *greedy-width* with the provided [script](scripts).
