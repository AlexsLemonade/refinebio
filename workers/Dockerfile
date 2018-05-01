FROM ccdl/data_refinery_worker_base

COPY common/dist/data-refinery-common-* common/

COPY workers/requirements.txt .
# The base image does not have Python 2.X on it at all, so all calls
# to pip or python are by default calls to pip3 or python3
RUN pip install -r requirements.txt

# Get the latest version from the dist directory.
RUN pip3 install common/$(ls common -1 | sort --version-sort | tail -1)

# We need JRE for QC tools
RUN apt-get update -y 
RUN apt-get install -y unzip default-jre --fix-missing

USER user

# Install FastQC
RUN wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
RUN unzip fastqc_v0.11.7.zip
RUN chmod +x ./FastQC/fastqc

# Install Aspera
RUN wget -q https://download.asperasoft.com/download/sw/cli/3.7.7/aspera-cli-3.7.7.608.927cce8-linux-64-release.sh
RUN sh aspera-cli-3.7.7.608.927cce8-linux-64-release.sh

COPY workers/ .

ENTRYPOINT []
