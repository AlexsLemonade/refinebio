FROM ccdl/data_refinery_worker_base

# It's annoying that this can only be installed via git.
RUN git clone https://github.com/deweylab/RSEM.git

RUN cd RSEM && make install

RUN rm -r RSEM

# Starting here, code is based on:
# https://github.com/COMBINE-lab/salmon/blob/master/docker/Dockerfile
ENV PACKAGES gcc make g++ cmake libboost-all-dev liblzma-dev libbz2-dev \
    ca-certificates zlib1g-dev unzip autoconf

# Version 0.8.2 has an outdated sha hash which breaks its make command.
# I fixed this in the lastest commit, but that hasn't been cut to a
# new version yet.
ENV SALMON_COMMIT_HASH 0b9f1eb8a0a321adc6320f490b015e75ebf3d3bb

# salmon binary will be installed in /home/salmon/bin/salmon

### don't modify things below here for version updates etc.

WORKDIR /home

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN curl -k -L https://github.com/COMBINE-lab/salmon/archive/${SALMON_COMMIT_HASH}.zip \
         -o salmon-${SALMON_COMMIT_HASH}.zip && \
    unzip salmon-${SALMON_COMMIT_HASH}.zip && \
    cd salmon-${SALMON_COMMIT_HASH} && \
    mkdir build && \
    cd build && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local && make && make install

ENV PATH /home/salmon-${SALMON_COMMIT_HASH}/bin:${PATH}
# End code from salmon repo.

COPY common/dist/data-refinery-common-* common/

# Get the latest version from the dist directory.
RUN pip3 install common/$(ls common -1 | sort --version-sort | tail -1)

COPY workers/ .

USER user

CMD ["celery", "-l", "info", "-A", "data_refinery_workers", "-c", "2", "worker", "-Ofair"]
