To initialize the database for the Data Refinery project run:

```
sudo -u postgres install.sh
```

To create tables within that database run:

```
python manage.py migrate
```

## Developing

After making changes to any model files within this directory two things should
happen before you can utilize those changes in other projects that rely on the
model. First you should make and run migrations:

```
./make_migrations.sh
```

Next you should generate the package:

```
python setup.py sdist
```

The Dockerfiles know how to find the package and install it in the container.
