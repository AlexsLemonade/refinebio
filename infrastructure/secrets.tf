
resource "aws_secretsmanager_secret" "django_secret_key" {
  name_prefix = "data-refinery-${var.user}-${var.stage}-django-secret-key"
}

resource "aws_secretsmanager_secret_version" "django_secret_key" {
  secret_id = aws_secretsmanager_secret.django_secret_key.id
  secret_string = var.django_secret_key
}

resource "aws_secretsmanager_secret" "database_password" {
  name_prefix = "data-refinery-${var.user}-${var.stage}-database-password"
}

resource "aws_secretsmanager_secret_version" "database_password" {
  secret_id = aws_secretsmanager_secret.database_password.id
  secret_string = var.database_password
}

resource "aws_secretsmanager_secret" "sentry_dsn" {
  name_prefix = "data-refinery-${var.user}-${var.stage}-raven-dsn"
}

resource "aws_secretsmanager_secret_version" "sentry_dsn" {
  secret_id = aws_secretsmanager_secret.sentry_dsn.id
  secret_string = var.sentry_dsn
}
