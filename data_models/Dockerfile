FROM python:3.6.1-slim

RUN groupadd user && useradd --create-home --home-dir /home/user -g user user
WORKDIR /home/user

COPY requirements.txt .

RUN pip install -r requirements.txt

COPY manage.py .

USER user

RUN ls

ENTRYPOINT ["python3.6", "manage.py"]

CMD ["migrate"]

# CMD sleep 10000
