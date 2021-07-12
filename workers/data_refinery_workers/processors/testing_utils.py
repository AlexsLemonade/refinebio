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

    def assertFailed(self, pj, expected_reason=""):
        pj.refresh_from_db()

        if not pj.success and expected_reason in pj.failure_reason:
            return

        if pj.success:
            raise self.failureException("Expected processor job to fail, but it succeeded")

        raise self.failureException(
            "Processor job failed, but for an incorrect reason.\n"
            + f"We expected the reason to contain '{expected_reason}' but the actual reason was '{pj.failure_reason}'"
        )
