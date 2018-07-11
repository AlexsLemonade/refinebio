import subprocess
import os
from data_refinery_common.models import OriginalFile

external_data_directory = os.curdir + "/data_store/external-data/"

if __name__ == "__main__":
    for filename in os.listdir(external_data_directory):
        if filename.endswith(".CEL"):
            file = OriginalFile()
            file.source_filename = filename
            file.filename = filename
            file.absolute_file_path = os.path.abspath(external_data_directory + filename)
            file.save()

            continue
        elif filename.endswith(".fastq.gz"):

            continue
        else:
            continue
