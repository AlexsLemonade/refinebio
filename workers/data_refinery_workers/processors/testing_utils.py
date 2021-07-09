class ProcessorJobTestCaseMixin:
    def assertSucceeded(self, pj):
        pj.refresh_from_db()

        if pj.success:
            return

        if pj.failure_reason is not None:
            msg = f"Processor job failed with reason '{pj.failure_reason}'"
        else:
            msg = f"Processor job failed without a given reason"

        raise self.failureException(msg)

    def assertFailed(self, pj):
        pj.refresh_from_db()

        if not pj.success:
            return

        raise self.failureException("Expected processor job to fail, but it succeeded")
