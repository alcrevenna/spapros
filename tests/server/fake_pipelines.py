"""Fast stand-ins for spapros.server.pipeline.run_pipeline, loaded by the job worker process."""

import time
from pathlib import Path


def succeed(job_dir, options, report_stage):
    report_stage("selecting")
    results = Path(job_dir) / "results"
    results.mkdir(exist_ok=True)
    (results / "probeset.csv").write_text("gene,selection\nA,True\nB,True\n")
    return {"genes": ["A", "B"], "n_genes": 2, "summary": {}, "files": ["probeset.csv"], "n": options["n"]}


def fail(job_dir, options, report_stage):
    report_stage("selecting")
    raise ValueError("boom")


def wait_for_release(job_dir, options, report_stage):
    """Runs until a ``release`` file appears next to the job directories, so tests control when it ends."""
    report_stage("waiting")
    release = Path(job_dir).parent / "release"
    for _ in range(600):
        if release.exists():
            return succeed(job_dir, options, report_stage)
        time.sleep(0.05)
    raise TimeoutError("never released")


def fail_once(job_dir, options, report_stage):
    """Leaves a checkpoint and fails on the first attempt; finishes from the checkpoint on the next one."""
    checkpoint = Path(job_dir) / "selection" / "checkpoint"
    if not checkpoint.exists():
        checkpoint.parent.mkdir(exist_ok=True)
        checkpoint.write_text("half done")
        raise RuntimeError("interrupted")
    result = succeed(job_dir, options, report_stage)
    result["resumed_from"] = checkpoint.read_text()
    return result
