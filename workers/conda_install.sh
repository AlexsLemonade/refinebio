bash miniconda.sh -b -p $HOME/miniconda
PATH="$HOME/miniconda/bin:$PATH"

conda install -c bioconda bioconductor-affy=1.50.0
conda install -c bioconda bioconductor-affyio=1.42.0

# rpy2 only uses miniconda/lib/R/lib, not miniconda/lib
# Since we can't change what rpy2 uses, just move the lib files.
cp miniconda/lib/lib* miniconda/lib/R/lib
