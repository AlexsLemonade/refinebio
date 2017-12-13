FROM ccdl/data_refinery_worker_base

COPY common/dist/data-refinery-common-* common/

# Get the latest version from the dist directory.
RUN pip3 install common/$(ls common -1 | sort --version-sort | tail -1)

RUN pip3 install python-nomad

COPY workers/ .

USER user

CMD ["celery", "-l", "info", "-A", "data_refinery_workers", "-c", "2", "worker", "-Ofair"]
