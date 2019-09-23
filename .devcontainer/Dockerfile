FROM python:3.6.1-slim

RUN mkdir /workspace
WORKDIR /workspace

RUN apt-get update
RUN apt-get -y install git procps wget unzip python3-pip jq iproute2 r-base r-base-dev sudo
RUN pip install --upgrade pip
RUN pip install -U pip-tools setuptools pylint pytest pytest-django pytest-pythonpath pytest-env

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


