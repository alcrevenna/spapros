import time

import pytest

pytest.importorskip("fastapi")

from spapros.server.jobs import JobQueue  # noqa: E402
from spapros.server.jobs import JobStore  # noqa: E402

FAKES = "tests.server.fake_pipelines"


def wait_for(store, job_id, statuses, timeout=60):
    """Poll until the job reaches one of ``statuses`` and return it."""
    statuses = {statuses} if isinstance(statuses, str) else set(statuses)
    deadline = time.time() + timeout
    while time.time() < deadline:
        job = store.get(job_id)
        if job["status"] in statuses:
            return job
        time.sleep(0.05)
    raise AssertionError(f"job {job_id} stayed {store.get(job_id)['status']}, expected {statuses}")


@pytest.fixture()
def store(tmp_path):
    return JobStore(tmp_path)


@pytest.fixture()
def make_queue(store):
    queues = []

    def make(fake, workers=1):
        q = JobQueue(store, pipeline=f"{FAKES}:{fake}", workers=workers, poll_interval=0.05)
        queues.append(q)
        return q

    yield make
    for q in queues:
        q.stop()
