import datetime

from django.conf import settings
from django.core.management.base import BaseCommand
from django.utils import timezone

import requests

from data_refinery_common.models import Dataset, DatasetAnnotation


class Command(BaseCommand):
    help = "Post downloads summary to slack"

    def add_arguments(self, parser):
        parser.add_argument(
            "--days",
            type=int,
            default=7,  # default to a week
            help=("Number of days in the past for which to build the stats"),
        )

    def handle(self, *args, **options):
        """
        """
        days = options["days"]
        start_time = timezone.now() - datetime.timedelta(days=-days)

        annotation_queryset = DatasetAnnotation.objects.filter(
            created_at__gt=start_time
        ).prefetch_related("dataset")
        annotations = [
            annotation
            for annotation in annotation_queryset
            if annotation.data["start"] and should_display_email(annotation.dataset.email_address)
        ]

        unique_users = list(set(annotation.dataset.email_address for annotation in annotations))
        unique_ips = list(set(annotation.data["ip"] for annotation in annotations))

        if unique_users:
            fallback_text = "In the last {0} days, {1} users downloaded datasets from {2} locations.".format(
                days, len(unique_users), len(unique_ips)
            )
        else:
            fallback_text = "There were no downloads in the last {0} days.".format(days)

        new_users = ""
        returning_users = ""
        for email in unique_users:
            user_annotations = annotation_queryset.filter(dataset__email_address=email)
            total_downloads = user_annotations.count()
            unique_locations = list(set(annotation.data["ip"] for annotation in user_annotations))
            locations = ", ".join(get_ip_location(ip) for ip in unique_locations)
            is_new_user = DatasetAnnotation.objects.filter(
                created_at__lt=start_time, dataset__email_address=email
            )
            text = "{0} | {1} downloads from {2}\n".format(email, total_downloads, locations)
            if is_new_user:
                new_users += text
            else:
                returning_users += text

        blocks = [
            {
                "type": "section",
                "text": {"type": "plain_text", "emoji": True, "text": fallback_text},
            }
        ]
        if new_users:
            blocks.append(
                {
                    "type": "section",
                    "text": {"type": "mrkdwn", "text": "*New users* \n" + new_users,},
                }
            )
        if returning_users:
            blocks.append(
                {
                    "type": "section",
                    "text": {"type": "mrkdwn", "text": "*Returning users* \n" + returning_users,},
                }
            )

        # Post to slack
        requests.post(
            "https://hooks.slack.com/services/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",  # settings.ENGAGEMENTBOT_WEBHOOK,
            json={
                "username": "EngagementBot",
                "icon_emoji": ":halal:",
                "channel": "#robots",  # "#ccdl-general",
                "text": fallback_text,
                "blocks": blocks,
            },
            headers={"Content-Type": "application/json"},
            timeout=10,
        )


def should_display_email(email: str) -> bool:
    """ Returns true if the given email is not associated with the CCDL suers """
    return (
        email is not None
        and email.find("cansav09") != 0
        and email.find("arielsvn") != 0
        and email.find("jaclyn.n.taroni") != 0
        and email.find("kurt.wheeler") != 0
        and email.find("greenescientist") != 0
        and email.find("@alexslemonade.org") == -1
        and email.find("miserlou") != 0
        and email.find("michael.zietz@gmail.com") != 0
        and email.find("d.prasad") != 0
        and email.find("daniel.himmelstein@gmail.com") != 0
        and email.find("dv.prasad991@gmail.com") != 0
    )


def get_ip_location(remote_ip):
    try:
        city = requests.get("https://ipapi.co/" + remote_ip + "/json/", timeout=10).json()["city"]
    except Exception:
        city = "COULD_NOT_DETERMINE"
    return city
