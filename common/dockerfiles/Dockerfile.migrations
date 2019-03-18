FROM python:3.6.1-slim

RUN groupadd user && useradd --create-home --home-dir /home/user -g user user
WORKDIR /home/user

# Needed for migration 0006
COPY config/ config

COPY common/requirements.txt .

RUN pip install -r requirements.txt

COPY common/ .

ARG SYSTEM_VERSION

ENV SYSTEM_VERSION $SYSTEM_VERSION

USER user

ENTRYPOINT [""]
