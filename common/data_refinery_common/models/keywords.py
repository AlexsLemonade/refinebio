from django.db import models


class SampleKeyword(models.Model):
    """An ontology term associated with a sample in our database"""

    name = models.ForeignKey("OntologyTerm", on_delete=models.CASCADE, related_name="+")
    sample = models.ForeignKey("Sample", on_delete=models.CASCADE, related_name="keywords")
    source = models.ForeignKey("Contribution", on_delete=models.CASCADE)
