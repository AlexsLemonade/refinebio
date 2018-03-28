"""
WSGI config for derpo project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
"""

import os

from django.core.wsgi import get_wsgi_application

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "data_refinery_api.settings")

application = get_wsgi_application()
