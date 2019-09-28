FROM ubuntu:16.04

# Avoid warnings by switching to noninteractive
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -qq
RUN apt-get install -y software-properties-common
RUN add-apt-repository ppa:apt-fast/stable
RUN apt-get update -qq
RUN apt-get -y install apt-fast

# The packages related to R are somewhat weird, see the README for more details.

COPY workers/CRAN.gpg /var/tmp
RUN \
  apt-fast update -qq && \
  apt-get install -y apt-transport-https && \
  apt-fast install -y lsb-release && \
  echo "deb https://cran.revolutionanalytics.com/bin/linux/ubuntu $(lsb_release -sc)/" \
      >> /etc/apt/sources.list.d/added_repos.list && \
  apt-key add /var/tmp/CRAN.gpg && \
  apt-fast update -qq && \
  apt-fast install -y \
  ed \
  git \
  mercurial \
  libcairo-dev \
  libedit-dev \
  lsb-release \
  python3 \
  python3-pip \
  r-base-core=3.4.2-1xenial1 \
  r-base-dev=3.4.2-1xenial1 \
  libpq-dev \
  libxml2-dev \
  libssl-dev \
  libcurl4-openssl-dev \
  curl \
  wget

# Create a folder to host the project
RUN mkdir /workspace
WORKDIR /workspace

RUN apt-get update \
    && apt-get -y install procps unzip jq iproute2
RUN pip3 install --upgrade pip
RUN pip3 install -U pip-tools setuptools pylint==2.3.1 pylint-django==2.0.11 pytest==5.0.1 pytest-django pytest-pythonpath pytest-env autopep8

# install git-crypt
RUN curl -L https://github.com/AGWA/git-crypt/archive/0.6.0.tar.gz | tar zxv -C /var/tmp
RUN cd /var/tmp/git-crypt-0.6.0 && make && make install PREFIX=/usr/local

# Install R packages (these are required to run all common tests)
ENV R_LIBS "/usr/local/lib/R/site-library"
COPY common/install_devtools.R /var/tmp
COPY workers/install_affy_only.R /var/tmp
RUN Rscript /var/tmp/install_devtools.R \
    && Rscript /var/tmp/install_affy_only.R

# Install this rpy2 here instead of via requirements.txt because
# pip-compile throws an error for it (got this from common/dockerfiles/Dockerfile.common_tests)
# rpy2 is required to run that code
RUN pip3 install rpy2==2.9.5

# install terraform
# RUN wget --quiet https://releases.hashicorp.com/terraform/0.12.9/terraform_0.12.9_darwin_amd64.zip \
#   && unzip terraform_0.12.9_darwin_amd64.zip \
#   && mv terraform /usr/bin \
#   && rm terraform_0.12.9_darwin_amd64.zip

# install Nomad (running it might be the issue)
# RUN wget --quiet "https://releases.hashicorp.com/nomad/0.8.3/nomad_0.8.3_linux_amd64.zip"\
#   && unzip -d /usr/bin "nomad_0.8.3_linux_amd64.zip" \
#   && chmod a+rx /usr/bin/nomad \
#   && rm "nomad_0.8.3_linux_amd64.zip"

# Install API project requirements
COPY api/requirements.txt /var/tmp
RUN pip3 install -r /var/tmp/requirements.txt

# Clean up
RUN apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

# Switch back to dialog for any ad-hoc use of apt-get
ENV DEBIAN_FRONTEND=