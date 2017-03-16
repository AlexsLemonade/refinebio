from downloaders.tasks import save_test_batch

# This is just a temporary way to test calling a task of the Celery worker.

save_test_batch.delay()
