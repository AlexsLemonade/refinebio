"""A script to calculate how many RAM-hours a job took.

A RAM-hour shall be defined as using 1GB of random access memory for 1 hour.

This is a useful metric for two reasons:
  * Some jobs are fast but take a lot of RAM, others might be slow but
    only take a bit of RAM. A job that takes 10 minutes using 6GB of
    RAM and a job that takes an hour using only 1GB of RAM cost us the
    same.
  * Based on the amount of RAM an instance has and its cost, we can
    determine the cost of a ram-hour and use that for estimating the
    costs of jobs.
"""

import argparse
import datetime

import pyrefinebio
import pytz


def parse_args():
    description = """This script calculates how many RAM-hours a sample has consumed.

    A RAM-hour shall be defined as using 1GB of random access memory for 1 hour.

    All downloader and processor jobs that were run will be
    considered, regardless of whether they were successful or not. The
    goal is to determine how many resources were spent to download and
    process a sample, regardless of success."""
    parser = argparse.ArgumentParser(description=description)

    accession_help_text = "Specify the accession_code of the sample to calculate RAM-hours for."
    parser.add_argument("accession_code", help=accession_help_text)

    start_date_help_text = """Do not include jobs started previous to start_date in the calculations.
    The format is expected to be YYYY-MM-DD."""
    parser.add_argument(
        "--start-date",
        default="2021-06-23",
        type=lambda s: datetime.datetime.strptime(s, "%Y-%m-%d"),
        help=start_date_help_text,
    )

    sample_help_text = """If supplied only calculate RAM-hours for a sample rather than an experiment.
    If this is supplied then the accession code provided should be a sample accession code."""
    parser.add_argument("--sample", help=sample_help_text, action="store_true")

    return parser.parse_args()


def calculate_ram_hours_for_job(job):
    if not (job.start_time and job.end_time and job.ram_amount):
        return 0

    run_time = job.end_time - job.start_time
    hours = run_time.total_seconds() / 60 / 60

    # RAM is stored in MB so convert to GB
    gigabytes_of_RAM = job.ram_amount / 1024

    return hours * gigabytes_of_RAM


def calculate_ram_hours_for_sample(accession_code, start_date):
    utc_start_date = pytz.utc.localize(start_date)
    sample = pyrefinebio.Sample.get(accession_code)

    if len(sample.original_files) == 0:
        print(f"Found no original files for sample {accession_code}")
        return 0

    downloader_ram_hours = 0
    processor_ram_hours = 0
    seen_downloader_job_ids = set()
    seen_processor_job_ids = set()
    for original_file_id in sample.original_files:
        original_file = pyrefinebio.OriginalFile.get(original_file_id)
        for downloader_job in original_file.downloader_jobs:
            if (
                downloader_job.id not in seen_downloader_job_ids
                and downloader_job.start_time
                and downloader_job.start_time > utc_start_date
            ):
                downloader_ram_hours += calculate_ram_hours_for_job(downloader_job)
                seen_downloader_job_ids.add(downloader_job.id)
        for processor_job in original_file.processor_jobs:
            if (
                processor_job.id not in seen_processor_job_ids
                and processor_job.start_time
                and processor_job.start_time > utc_start_date
            ):
                processor_ram_hours += calculate_ram_hours_for_job(processor_job)
                seen_processor_job_ids.add(processor_job.id)

    total_ram_hours = downloader_ram_hours + processor_ram_hours
    return [accession_code, downloader_ram_hours, processor_ram_hours, total_ram_hours]


def calculate_ram_hours_for_experiment(accession_code, start_date):
    sample_data = []
    experiment = pyrefinebio.Experiment.get(accession_code)
    for sample in experiment.samples:
        sample_data.append(calculate_ram_hours_for_sample(sample.accession_code, start_date))

    return sample_data


if __name__ == "__main__":
    args = parse_args()

    print(
        ", ".join(
            ["accession_code", "downloader_ram_hours", "processor_ram_hours", "total_ram_hours"]
        )
    )
    if args.sample:
        sample_data = calculate_ram_hours_for_sample(args.accession_code, args.start_date)
        print(", ".join(list(map(str, sample_data))))
    else:
        all_sample_data = calculate_ram_hours_for_experiment(args.accession_code, args.start_date)
        for sample_data in all_sample_data:
            print(", ".join(list(map(str, sample_data))))

    # print(f"Sample {args.accession_code} took {ram_hours} RAM-hours to process.")
