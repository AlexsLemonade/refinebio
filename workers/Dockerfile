FROM ccdl/data_refinery_worker_base

# It's annoying that this can only be installed via git.
RUN git clone https://github.com/deweylab/RSEM.git

RUN cd RSEM && make install

RUN rm -r RSEM

COPY common/dist/data-refinery-common-* common/

# Get the latest version from the dist directory.
RUN pip3 install common/$(ls common -1 | sort --version-sort | tail -1)

COPY workers/ .

USER user

CMD ["celery", "-l", "info", "-A", "data_refinery_workers", "-c", "2", "worker", "-Ofair"]
