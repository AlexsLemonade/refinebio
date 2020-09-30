import pyrefinebio.common.annotation as prb_annotation
import pyrefinebio.sample as prb_sample
from pyrefinebio.http import get_by_endpoint
from pyrefinebio.util import generator_from_pagination


class Experiment:
    def __init__(
        self,
        id=None,
        title=None,
        description=None,
        annotations=[],
        samples=[],
        protocol_description=None,
        accession_code=None,
        alternate_accession_code=None,
        source_database=None,
        source_url=None,
        has_publication=None,
        publication_title=None,
        publication_doi=None,
        publication_authors=None,
        pubmed_id=None,
        source_first_published=None,
        source_last_modified=None,
        submitter_institution=None,
        last_modified=None,
        created_at=None,
        organism_names=None,
        sample_metadata=None,
        num_total_samples=None,
        num_processed_samples=None,
        num_downloadable_samples=None,
    ):
        self.id = id
        self.title = title
        self.description = description
        self.annotations = [prb_annotation.Annotation(**annotation) for annotation in annotations]
        self.samples = [prb_sample.Sample(**sample) for sample in samples]
        self.protocol_description = protocol_description
        self.accession_code = accession_code
        self.alternate_accession_code = alternate_accession_code
        self.source_database = source_database
        self.source_url = source_url
        self.has_publication = has_publication
        self.publication_title = publication_title
        self.publication_doi = publication_doi
        self.publication_authors = publication_authors
        self.pubmed_id = pubmed_id
        self.source_first_published = source_first_published
        self.source_last_modified = source_last_modified
        self.submitter_institution = submitter_institution
        self.last_modified = last_modified
        self.created_at = created_at
        self.organism_names = organism_names
        self.sample_metadata = sample_metadata
        self.num_total_samples = num_total_samples
        self.num_processed_samples = num_processed_samples
        self.num_downloadable_samples = num_downloadable_samples

    @classmethod
    def get(cls, accession_code):
        response = get_by_endpoint("experiments/" + accession_code)
        return Experiment(**response)

    @classmethod
    def search(cls, **kwargs):
        response = get_by_endpoint("search")
        return generator_from_pagination(response, cls)
