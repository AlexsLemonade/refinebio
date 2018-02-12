# Data Refinery Common

This sub-project contains code common to other sub-projects.

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

The Dockerfiles for other projects know how to find the package and install it
in the container.
