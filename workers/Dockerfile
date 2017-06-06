FROM continuumio/miniconda3:4.3.11

RUN groupadd user && useradd --create-home --home-dir /home/user -g user user
WORKDIR /home/user

RUN apt-get install libcurl4-openssl-dev

COPY workers/environment.yml .

RUN conda env create -f environment.yml

ENV PATH "/opt/conda/envs/workers/bin:$PATH"

COPY data_models/dist/data-refinery-models-* data_models/

# Get the latest version from the dist directory.
RUN pip install data_models/$(ls data_models -1 | tail -1)

COPY common/dist/data-refinery-common-* common/

# Get the latest version from the dist directory.
RUN pip install common/$(ls common -1 | tail -1)

COPY workers/ .

USER user

CMD ["celery", "-l", "info", "-A", "data_refinery_workers", "-c", "1", "worker", "-Ofair"]
