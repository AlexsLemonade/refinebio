FROM ccdl/data_refinery_worker_base

COPY data_models/dist/data-refinery-models-* data_models/

# Get the latest version from the dist directory.
RUN pip3 install data_models/$(ls data_models -1 | sort --version-sort | tail -1)

COPY common/dist/data-refinery-common-* common/

# Get the latest version from the dist directory.
RUN pip3 install common/$(ls common -1 | sort --version-sort | tail -1)

COPY workers/ .

USER user

ENTRYPOINT ["python3", "manage.py"]
