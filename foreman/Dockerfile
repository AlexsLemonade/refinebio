FROM python:3.6.1-slim

RUN groupadd user && useradd --create-home --home-dir /home/user -g user user
WORKDIR /home/user

COPY foreman/requirements.txt .

RUN pip install -r requirements.txt

COPY data_models/dist/data-refinery-models-* data_models/

# Get the latest version from the dist directory.
RUN pip install \
    data_models/$(ls data_models -1 | sort --version-sort | tail -1)

COPY common/dist/data-refinery-common-* common/

# Get the latest version from the dist directory.
RUN pip install \
    common/$(ls common -1 | sort --version-sort | tail -1)

COPY foreman/ .

USER user

ENTRYPOINT ["python3.6", "manage.py"]

CMD ["survey"]
