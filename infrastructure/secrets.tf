
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

resource "aws_secretsmanager_secret" "raven_dsn" {
  name_prefix = "data-refinery-${var.user}-${var.stage}-raven-dsn"
}

resource "aws_secretsmanager_secret_version" "raven_dsn" {
  secret_id = aws_secretsmanager_secret.raven_dsn.id
  secret_string = var.raven_dsn
}

resource "aws_secretsmanager_secret" "raven_dsn_api" {
  name_prefix = "data-refinery-${var.user}-${var.stage}-raven-dsn-api"
}

resource "aws_secretsmanager_secret_version" "raven_dsn_api" {
  secret_id = aws_secretsmanager_secret.raven_dsn_api.id
  secret_string = var.raven_dsn_api
}
