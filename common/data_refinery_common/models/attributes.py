from django.db import models

from data_refinery_common.models.ontology_term import OntologyTerm

BOOL = "bool"
INT = "int"
FLOAT = "float"
ONTOLOGY_TERM = "ont"


class AbstractAttribute(models.Model):
    """This is an abstract class that defines all of the properties of a single
    attribute on either a sample or an experiment. We then subclass this to
    associate attributes with either an experiment or a sample."""

    name = models.ForeignKey("OntologyTerm", on_delete=models.CASCADE, related_name="+")
    unit = models.ForeignKey("OntologyTerm", on_delete=models.CASCADE, related_name="+", null=True)
    probability = models.FloatField(null=True)
    source = models.ForeignKey("Contribution", on_delete=models.CASCADE)

    value_type = models.TextField(
        choices=[
            (BOOL, "Boolean"),
            (INT, "Integer"),
            (FLOAT, "Float"),
            (ONTOLOGY_TERM, "Ontology Term"),
        ],
    )
    value = models.TextField()

    class Meta:
        abstract = True

    def to_dict(self):
        rendered = {
            "name": self.name.to_dict(),
            "unit": None if self.unit is None else self.unit.to_dict(),
            "probability": "unknown" if self.probability is None else self.probability,
            "source": self.source.source_name,
            "methods": self.source.methods_url,
            "value": self.get_value(),
        }

        if self.value_type == ONTOLOGY_TERM:
            rendered["value"] = rendered["value"].to_dict()

        return rendered

    def set_value(self, value):
        """This method sets the attribute value and assigns the correct
        value_type. NOTE: we assume that all provided strings are ontology terms."""

        if type(value) == str:
            self.value_type = ONTOLOGY_TERM

            # make sure that we know how to deal with this ontology term, and error out if we don't
            OntologyTerm.get_or_create_from_api(value)

        elif type(value) == bool:
            self.value_type = BOOL
        elif type(value) == int:
            self.value_type = INT
        elif type(value) == float:
            self.value_type = FLOAT
        else:
            raise ValueError("Invalid metadata value type '{}'".format(type(value)))

        self.value = str(value)

    def get_value(self):
        """This method returns the value of this attribute using `value_type`
        to convert to the same type"""

        if self.value_type == ONTOLOGY_TERM:
            return OntologyTerm.get_or_create_from_api(self.value)
        elif self.value_type == BOOL:
            return bool(self.value)
        elif self.value_type == INT:
            return int(self.value)
        elif self.value_type == FLOAT:
            return float(self.value)
        else:
            raise ValueError("Invalid value_type '{}'".format(self.value_type))


class SampleAttribute(AbstractAttribute):
    sample = models.ForeignKey("Sample", on_delete=models.CASCADE, related_name="attributes")


class ExperimentAttribute(AbstractAttribute):
    experiment = models.ForeignKey(
        "Experiment", on_delete=models.CASCADE, related_name="attributes"
    )
