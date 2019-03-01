FROM ubuntu:16.04

RUN apt-get update -qq
RUN apt-get install -y software-properties-common
RUN add-apt-repository ppa:apt-fast/stable
RUN apt-get update -qq
RUN apt-get -y install apt-fast

# The packages related to R are somewhat weird, see the README for more details.

COPY workers/CRAN.gpg .
RUN \
  apt-fast update -qq && \
  apt-get install -y apt-transport-https && \
  apt-fast install -y lsb-release && \
  echo "deb https://cran.revolutionanalytics.com/bin/linux/ubuntu $(lsb_release -sc)/" \
      >> /etc/apt/sources.list.d/added_repos.list && \
  apt-key add CRAN.gpg && \
  apt-fast update -qq && \
  apt-fast install -y \
  ed \
  git \
  python3 \
  python3-pip \
  r-base-core=3.4.2-1xenial1 \
  r-base-dev=3.4.2-1xenial1 \
  libpq-dev \
  libxml2-dev \
  libssl-dev \
  libcurl4-openssl-dev \
  libpq-dev \
  curl \
  wget && \
  rm -rf /var/lib/apt/lists/*
RUN rm CRAN.gpg

RUN groupadd user && useradd --create-home --home-dir /home/user -g user user
WORKDIR /home/user

ENV R_LIBS "/usr/local/lib/R/site-library"

COPY common/install_devtools.R .

RUN Rscript install_devtools.R

COPY workers/install_affy_only.R .

RUN Rscript install_affy_only.R

RUN pip3 install --upgrade pip

COPY config config

COPY common/requirements.txt .

RUN pip3 install -r requirements.txt

# Install this rpy2 here instead of via requirements.txt because
# pip-compile throws an error for it.
RUN pip3 install rpy2==2.9.5

COPY common/ .

ARG SYSTEM_VERSION

ENV SYSTEM_VERSION $SYSTEM_VERSION

USER user

ENTRYPOINT [""]
