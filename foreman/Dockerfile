FROM python:3.6.1-slim

RUN groupadd user && useradd --create-home --home-dir /home/user -g user user
WORKDIR /home/user

COPY foreman/requirements.txt .

# The base image does not have Python 2.X on it at all, so all calls
# to pip or python are by default calls to pip3 or python3
RUN pip install -r requirements.txt

COPY common/dist/data-refinery-common-* common/

# Get the latest version from the dist directory.
RUN pip install \
    common/$(ls common -1 | sort --version-sort | tail -1)

COPY foreman/ .

USER user

ENTRYPOINT ["python3.6", "manage.py"]
