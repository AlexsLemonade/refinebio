from django.db import models


class Contributor(models.Model):
    """This model represents an external contributer and stores their name so
    we can use it as a foreign key on information (e.g. keywords and
    attributes) supplied by that contributor
    """

    name = models.TextField(unique=True)
