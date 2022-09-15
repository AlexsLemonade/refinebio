"""Abstract base class for accession gathering automation agents."""

from abc import ABC, abstractmethod
from datetime import datetime
from http.client import RemoteDisconnected

from requests.exceptions import ConnectionError, ConnectTimeout
from urllib3.exceptions import ProtocolError

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models.experiment import Experiment
from data_refinery_common.models.gathered_accession import GatheredAccession

logger = get_and_configure_logger(__name__)


class AccessionAgentBase(ABC):
    "Accession agent abstract base class."

    previous_accessions = set()
    retry_params = {
        "retry_on_exception": lambda e: isinstance(
            e, (ConnectionError, ConnectTimeout, ProtocolError, RemoteDisconnected)
        ),
        "stop_max_attempt_number": 5,
        "wait_exponential_multiplier": 1000,  # Seconds.
        "wait_exponential_max": 16000,  # Seconds.
    }

    def __init__(self, options) -> None:
        """Populates args and values for major variables."""
        self.options = options
        self.count = options["count"]
        self.keyword = options["keyword"]
        self.organism = options["organism"]
        self.since = options["since"]
        self.until = options["until"] or datetime.now().strftime("%Y-%m-%d")

        self.ids = self.get_ids()
        self.populate_previous_accessions()

    @abstractmethod
    def build_query(self):
        """Returns query/query dict depending on the accession data source."""
        pass

    @abstractmethod
    def collect_data(self):
        """Generates resulting entry collection."""
        pass

    @abstractmethod
    def fetch_data(self):
        """Fetches data from an external or local data source."""
        pass

    @abstractmethod
    def get_ids(self):
        """Gets IDs for query filtering depending on the accession technology."""
        pass

    def populate_previous_accessions(self) -> None:
        """Populates previous accession set from a provided excluded ids file."""
        if not self.options["exclude_previous"] or self.previous_accessions:
            return

        # Gathered accessions.
        self.previous_accessions.update(
            (
                entry["accession_code"]
                for entry in GatheredAccession.objects.values("accession_code")
            )
        )

        # Surveyed accessions.
        experiments = Experiment.objects.values("accession_code", "alternate_accession_code")
        self.previous_accessions.update(
            (experiment["accession_code"] for experiment in experiments)
        )
        self.previous_accessions.update(
            (experiment["alternate_accession_code"] for experiment in experiments)
        )
