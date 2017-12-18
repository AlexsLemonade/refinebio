"""Provides an interface to send messages to the Message Queue."""

from __future__ import absolute_import, unicode_literals
import nomad


def send_job(job_name: str, job_id: int):
    nomad_client = nomad.Nomad("database", timeout=5)
    nomad_client.job.dispatch_job(job_name, meta={"JOB_ID": str(job_id)})
