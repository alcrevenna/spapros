import json

from spapros.server.jobs import JobQueue
from tests.server.conftest import wait_for


def submit(store, queue, n=3):
    job = store.create({"n": n}, input_filename="data.h5ad")
    queue.submit(job["id"])
    return job["id"]


def test_job_succeeds_and_stores_result(store, make_queue):
    queue = make_queue("succeed")
    queue.start()
    job_id = submit(store, queue, n=7)

    job = wait_for(store, job_id, "succeeded")
    assert job["attempts"] == 1
    assert job["stage"] == "selecting"
    assert job["started_at"] and job["finished_at"]
    assert store.result(job_id)["n"] == 7
    assert (store.path(job_id) / "results" / "probeset.csv").exists()
    assert "stage: selecting" in (store.path(job_id) / "log.txt").read_text()


def test_job_failure_reports_error(store, make_queue):
    queue = make_queue("fail")
    queue.start()
    job_id = submit(store, queue)

    job = wait_for(store, job_id, "failed")
    assert job["error"] == "ValueError: boom"
    assert "Traceback" in (store.path(job_id) / "error.txt").read_text()
    assert store.result(job_id) is None


def test_jobs_run_one_at_a_time_in_order(store, make_queue, tmp_path):
    queue = make_queue("wait_for_release")
    queue.start()
    first, second = submit(store, queue), submit(store, queue)

    wait_for(store, first, "running")
    assert store.get(second)["status"] == "queued"
    (store.root / "release").touch()
    wait_for(store, first, "succeeded")
    wait_for(store, second, "succeeded")


def test_cancel_running_and_queued_jobs(store, make_queue):
    queue = make_queue("wait_for_release")
    queue.start()
    running, queued = submit(store, queue), submit(store, queue)
    wait_for(store, running, "running")

    assert queue.cancel(queued)["status"] == "cancelled"
    queue.cancel(running)
    assert wait_for(store, running, "cancelled")["finished_at"]
    # The cancelled queued job is skipped, not run.
    assert store.get(queued)["attempts"] == 0


def test_retry_resumes_from_checkpoint(store, make_queue):
    queue = make_queue("fail_once")
    queue.start()
    job_id = submit(store, queue)
    assert wait_for(store, job_id, "failed")["error"] == "RuntimeError: interrupted"

    queue.retry(job_id)
    job = wait_for(store, job_id, "succeeded")
    assert job["attempts"] == 2
    assert job["error"] is None
    assert store.result(job_id)["resumed_from"] == "half done"
    assert not (store.path(job_id) / "error.txt").exists()


def test_retry_rejects_unfinished_job(store, make_queue):
    queue = make_queue("succeed")
    job_id = submit(store, queue)  # queue not started, so the job stays queued
    try:
        queue.retry(job_id)
    except ValueError as e:
        assert "queued" in str(e)
    else:
        raise AssertionError("retry of a queued job should fail")


def test_restart_requeues_unfinished_jobs(store, make_queue):
    # A server that stops while a job runs leaves it "running" ...
    queue = make_queue("wait_for_release")
    queue.start()
    job_id = submit(store, queue)
    wait_for(store, job_id, "running")
    queue.stop()
    assert store.get(job_id)["status"] == "running"

    # ... and the next server picks it up again.
    (store.root / "release").touch()
    restarted = JobQueue(store, pipeline="tests.server.fake_pipelines:wait_for_release", poll_interval=0.05)
    restarted.start()
    try:
        job = wait_for(store, job_id, "succeeded")
    finally:
        restarted.stop()
    assert job["attempts"] == 2


def test_job_json_is_valid_after_updates(store):
    job = store.create({"n": 1}, input_filename="x.h5ad")
    store.update(job["id"], status="running", attempts=1)
    data = json.loads((store.path(job["id"]) / "job.json").read_text())
    assert data["status"] == "running" and data["attempts"] == 1
    assert store.get("../etc") is None
