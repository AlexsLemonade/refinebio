FROM python:3.6.1-slim

RUN groupadd user && useradd --create-home --home-dir /home/user -g user user
WORKDIR /home/user

COPY common/requirements.txt .

RUN pip install -r requirements.txt

COPY common/ .

USER user

ENTRYPOINT [""] 
