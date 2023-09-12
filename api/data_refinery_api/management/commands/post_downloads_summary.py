import datetime
from collections import Counter

from django.conf import settings
from django.core.management.base import BaseCommand
from django.template.defaultfilters import pluralize
from django.utils import timezone

import requests

from data_refinery_common.models import Dataset, DatasetAnnotation


class Command(BaseCommand):
    help = "Post downloads summary to Slack"

    def add_arguments(self, parser):
        parser.add_argument(
            "--channel",
            type=str,
            default="ccdl-general",
            help=("Optional parameter to choose the channel where the message will be posted."),
        )
        parser.add_argument(
            "--days",
            type=int,
            default=7,  # Default to a week.
            help=("Number of days in the past for which to build the stats."),
        )
        parser.add_argument(
            "--top-countries",
            type=int,
            default=5,
            help=("Number of countries to show in the per country downloads summary."),
        )

    def handle(self, *args, **options):
        post_downloads_summary(options["days"], options["channel"], options["top_countries"])


def format_user_data(header, data):
    """
    Formats user email, downloads count, location information sorted
    by downloads count.
    """
    # Allowed overhead for 2 column sorting: downloads count, email.
    lines = sorted(data, key=lambda u: u[0].lower())
    lines = [
        f"{email.lower()} | {downloads} download{pluralize(downloads)} from {location}"
        for email, downloads, location in sorted(lines, key=lambda u: u[1], reverse=True)
    ]
    lines.insert(0, header)

    return "\n".join(lines)


def get_user_location(ip_address):
    """Gets user location information based on their IP address."""
    try:
        data = requests.get(f"https://ipapi.co/{ip_address}/json/", timeout=10).json()
        # The list of available fields https://ipapi.co/api/#complete-location
        return ", ".join((data["city"], data["country_name"]))
    except (requests.exceptions.RequestException, KeyError, ValueError):
        return ip_address


def post_downloads_summary(days, channel, top_countries=5):
    """Posts downloads summary to Slack channel."""
    start_time = timezone.now() - datetime.timedelta(days=days)
    datasets = Dataset.processed_filtered_objects.filter(
        created_at__gt=start_time
    ).prefetch_related("datasetannotation_set")
    annotations = DatasetAnnotation.objects.filter(dataset__in=datasets)
    users_emails = set(dataset.email_address for dataset in datasets)

    locations = set()
    locations_cache = {}
    for annotation in annotations:
        if "location" not in annotation.data:
            ip_address = annotation.data["ip"]
            if ip_address not in locations_cache:
                locations_cache[ip_address] = get_user_location(ip_address)

            # Save the locations permanently, since IP addresses can cycle over time.
            annotation.data["location"] = locations_cache[ip_address]
            annotation.save()
        locations.add(annotation.data["location"])

    downloads_per_country = Counter()
    downloads_total = 0
    new_users = []
    returning_users = []
    for user_email in users_emails:
        user_annotations = annotations.filter(dataset__email_address=user_email)
        user_downloads = user_annotations.count()
        if user_downloads == 0:
            continue

        downloads_total += user_downloads
        user_locations = set()
        for user_annotation in user_annotations:
            user_locations.add(user_annotation.data["location"])
            try:
                country = user_annotation.data["location"].split(", ")[1]
                downloads_per_country.update({country: 1})
            except (IndexError, TypeError):
                pass

        user_locations = "; ".join(sorted(user_locations))
        user_data = (user_email, user_downloads, user_locations)

        is_returning_user = Dataset.processed_filtered_objects.filter(
            created_at__lt=start_time, email_address=user_email
        ).exists()
        if is_returning_user:
            returning_users.append(user_data)
        else:
            new_users.append(user_data)

    if downloads_total > 0:
        locations_count = len(locations)
        users_count = len(new_users) + len(returning_users)
        fallback_text = (
            f"In the last {days} day{pluralize(days)}, {users_count} "
            f"user{pluralize(users_count)} downloaded {downloads_total} "
            f"dataset{pluralize(downloads_total)} from {locations_count} "
            f"location{pluralize(locations_count)}."
        )
    else:
        fallback_text = f"There were no downloads in the last {days} day{pluralize(days)}."

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
                "text": {
                    "type": "mrkdwn",
                    "text": format_user_data("*New users*", new_users),
                },
            }
        )

    if returning_users:
        blocks.append(
            {
                "type": "section",
                "text": {
                    "type": "mrkdwn",
                    "text": format_user_data("*Returning users*", returning_users),
                },
            }
        )

    if top_countries and downloads_per_country:
        countries_count = downloads_per_country.most_common(top_countries)
        top_countries = min(top_countries, len(countries_count))
        lines = [f"*Top {top_countries} countr{pluralize(top_countries, 'y,ies')}*"]
        # Allowed overhead for 2 column sorting: downloads count, country.
        countries_count = sorted(countries_count, key=lambda cc: cc[0])
        countries_count = sorted(countries_count, key=lambda cc: cc[1], reverse=True)
        for country, count in countries_count:
            lines.append(f"{country}: {count} download{pluralize(count)}")

        blocks.append(
            {
                "type": "section",
                "text": {"type": "mrkdwn", "text": "\n".join(lines)},
            }
        )

    # Post to Slack.
    requests.post(
        settings.ENGAGEMENTBOT_WEBHOOK,
        json={
            "username": "EngagementBot",
            "icon_emoji": ":halal:",
            "channel": f"#{channel}",
            "text": fallback_text,
            "blocks": blocks,
        },
        headers={"Content-Type": "application/json"},
        timeout=10,
    )
