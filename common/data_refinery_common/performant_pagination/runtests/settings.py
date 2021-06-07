#
#
#

from __future__ import absolute_import, print_function, unicode_literals

DEBUG = True
TEMPLATE_DEBUG = DEBUG
DEBUG_PROPAGATE_EXCEPTIONS = True

# https://dev.to/rubyflewtoo/upgrading-to-django-3-2-and-fixing-defaultautofield-warnings-518n
DEFAULT_AUTO_FIELD = "django.db.models.AutoField"

DATABASES = {"default": {"ENGINE": "django.db.backends.sqlite3", "NAME": "sqlite.db",}}

SECRET_KEY = "p!_%jd*snxvam@gs2(+w*mzg+x231s37rz#v*=f&kz07nn81j5"

INSTALLED_APPS = (
    "performant_pagination",
    "performant_pagination.tests",
)
