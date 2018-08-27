FROM python:3.6.1

RUN groupadd user && useradd --create-home --home-dir /home/user -g user user
WORKDIR /home/user

COPY config/ config/

COPY api/requirements.txt .

RUN pip install -r requirements.txt

COPY common/dist/data-refinery-common-* common/

# Get the latest version from the dist directory.
RUN pip install \
    common/$(ls common -1 | sort --version-sort | tail -1)

COPY api/ .

# uWSGI
RUN pip install uwsgi
RUN chmod +x /home/user/collect_and_run_uwsgi.sh

RUN mkdir -p /tmp/www/static
RUN chown user /tmp/www/static

ARG SYSTEM_VERSION

ENV SYSTEM_VERSION $SYSTEM_VERSION

USER user

# We collect Django's static files and expose them
# as a volume so that Nginx can serve them directly.
VOLUME /tmp/www/static

EXPOSE 8081

CMD ["/home/user/collect_and_run_uwsgi.sh"]
