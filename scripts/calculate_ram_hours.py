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

    return parser.parse_args()


def calculate_ram_hours_for_job(job):
    print("Hi")
    if not (job.start_time and job.end_time and job.ram_amount):
        print(job.start_time)
        print(job.end_time)
        print(job.ram_amount)
        return 0

    run_time = job.end_time - job.start_time
    hours = run_time.total_seconds() / 60 / 60

    # RAM is stored in MB so convert to GB
    gigabytes_of_RAM = job.ram_amount / 1024

    return hours * gigabytes_of_RAM


def calculate_ram_hours(accession_code, start_date):
    utc_start_date = pytz.utc.localize(start_date)
    sample = pyrefinebio.Sample.get(accession_code)

    if len(sample.original_files) == 0:
        print(f"Found no original files for sample {accession_code}")
        return 0

    ram_hours = 0
    for original_file_id in sample.original_files:
        original_file = pyrefinebio.OriginalFile.get(original_file_id)

        for downloader_job in original_file.downloader_jobs:
            print(downloader_job.start_time)
            if downloader_job.start_time and downloader_job.start_time > utc_start_date:
                ram_hours += calculate_ram_hours_for_job(downloader_job)

        for processor_job in original_file.processor_jobs:
            print(processor_job.start_time)
            if processor_job.start_time and processor_job.start_time > utc_start_date:
                ram_hours += calculate_ram_hours_for_job(processor_job)

    return ram_hours
    # print(accession_code)
    # print(start_date)


if __name__ == "__main__":
    args = parse_args()

    ram_hours = calculate_ram_hours(args.accession_code, args.start_date)

    print(f"Sample {args.accession_code} took {ram_hours} RAM-hours to process.")
