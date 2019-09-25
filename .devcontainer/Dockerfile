FROM python:3.6.1

# Avoid warnings by switching to noninteractive
# ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update
RUN apt-get -y install git procps lsb-release wget curl unzip python3-pip jq iproute2 r-base r-base-dev
RUN pip install --upgrade pip
RUN pip install -U pip-tools setuptools pylint pylint-django pytest==5.0.1 pytest-django pytest-pythonpath pytest-env

# install git-crypt
RUN curl -L https://github.com/AGWA/git-crypt/archive/0.6.0.tar.gz | tar zxv -C /var/tmp
RUN cd /var/tmp/git-crypt-0.6.0 && make && make install PREFIX=/usr/local

RUN mkdir /workspace
WORKDIR /workspace

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
COPY api/requirements.txt .
RUN pip install -r requirements.txt && rm requirements.txt

# Get the latest version from the dist directory.
# COPY common/dist/data-refinery-common-* common/
# RUN pip install \
#     common/$(ls common -1 | sort --version-sort | tail -1) \
#     && rm -r common/

# Clean up
RUN apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*


# Switch back to dialog for any ad-hoc use of apt-get
ENV DEBIAN_FRONTEND=