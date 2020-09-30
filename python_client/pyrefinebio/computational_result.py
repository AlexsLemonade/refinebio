import pyrefinebio.common.annotation as prb_annotation
import pyrefinebio.common.organism_index as prb_organism_index
import pyrefinebio.computed_file as prb_computed_file
import pyrefinebio.processor as prb_processor
from pyrefinebio.http import get_by_endpoint
from pyrefinebio.util import generator_from_pagination


class ComputationalResult:
    def __init__(
        self,
        id=None,
        commands=None,
        processor=None,
        is_ccdl=None,
        annotations=[],
        files=[],
        organism_index=None,
        time_start=None,
        time_end=None,
        created_at=None,
        last_modified=None,
    ):
        self.id = id
        self.commands = commands
        self.processor = prb_processor.Processor(**(processor or {}))
        self.is_ccdl = is_ccdl
        self.annotations = [prb_annotation.Annotation(**annotation) for annotation in annotations]
        self.files = [prb_computed_file.ComputedFile(**file) for file in files]
        self.organism_index = prb_organism_index.OrganismIndex(**(organism_index or {}))
        self.time_start = time_start
        self.time_end = time_end
        self.created_at = created_at
        self.last_modified = last_modified

    @classmethod
    def get(cls, id):
        result = get_by_endpoint("computational_results/" + str(id))
        return ComputationalResult(**result)

    @classmethod
    def search(cls, **kwargs):
        result = get_by_endpoint("computational_results", params=kwargs)
        return generator_from_pagination(result, cls)
