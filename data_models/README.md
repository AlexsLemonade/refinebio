To initialize the database for the Bioinformatics Mill project run:

```
sudo -u postgres install.sh
```

To create tables within that database run:

```
python manage.py migrate
```

## Developing

After making changes to any model files within this directory two things should happen before you can utilize those changes in other projects that rely on the model. First you should make and run migrations:

```
python manage.py makemigrations
python manage.py migrate
```

Next you should install the package locally:

```
python setup.py sdist
pip install --user dist/bioinformatics-mill-models-0.1.tar.gz
```

After merging your changes into the master branch be sure to release the package with:

```
twine upload dist/*
```

(If you are not already set up within pypi.org to be able to modify this package please contact @kurtwheeler.)
