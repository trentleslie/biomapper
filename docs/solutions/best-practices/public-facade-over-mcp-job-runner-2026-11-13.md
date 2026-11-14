---
title: "Opening an MCP-only job runner as public library API: a facade, caller-side batching for progress, and a restart test that separates saved from in-flight"
date: 2026-11-13
category: best-practices/
module: biomapper/jobs
problem_type: best_practice
component: background_job
severity: medium
applies_when:
  - "A second consumer (a web backend) needs a resumable job runner that so far only an MCP server can reach"
  - "Adding per-batch progress to a function whose mapper is injected by the caller"
  - "Writing a test that a job killed mid-run does not repeat finished work"
  - "A job outlives its process and its result must be readable afterwards"
related_components:
  - testing_framework
  - development_workflow
tags:
  - background-jobs
  - harmonize
  - public-api
  - progress-callback
  - resume
  - anyio
  - release
  - mcp
---

# Opening an MCP-only job runner as public library API

Shipped in biomapper 1.9.0 (fork PR trentleslie/biomapper#49, merge commit `793e7e6`).

## Context

biomapper-ui needed to run `harmonize_cohorts()` and `harmonize_many()` as long background jobs
from a FastAPI backend. Four things stood in the way, all confirmed against `main` at `81ccb83`:

- The resumable wrapper (gated mapper, per-job entity cache, staging and single placement,
  review items) lived in `src/biomapper/mcp/harmonize_jobs.py` and `src/biomapper/mcp/jobs.py`,
  so only the MCP server could reach it.
- `on_message` carried only the "run directory" and "wrote results" lines. Per-entity progress
  went to stderr.
- The report was an in-memory object plus TSV tables. Nothing held the whole report as data.
- `harmonize_cohorts()` re-mapped both cohorts on every call.

The work was to close those gaps without a second implementation of anything, and to leave
`harmonize()`'s linking rules untouched.

## Guidance

### 1. Put a facade over the existing runner. Do not write a second one.

`biomapper.jobs.HarmonizeJobs` (`src/biomapper/jobs.py`) builds the same `ServerConfig`,
`RequestGate`, `BioMapperClient` and `JobManager` the MCP server's lifespan builds, and calls the
same `submit` / `status` / `resume` / `cancel`. The runner itself gained only a progress hook and
one field. Two details make the facade safe:

**Shutdown order must mirror the server's lifespan.** A harmonize job runs a synchronous function
on a worker thread. The order on exit is: set the stop flag, cancel the worker, wait for the task
group (which waits for the thread), release the job locks, close the client. With an
`AsyncExitStack`, callbacks unwind in reverse, so they are registered backwards:

```python
client = await stack.enter_async_context(BioMapperClient(...))
manager = JobManager(cfg, gate=RequestGate(), client=client)
stack.callback(manager.close)                      # runs 4th: locks, after the thread is gone
group = await stack.enter_async_context(anyio.create_task_group())   # exits 3rd: waits
stack.callback(group.cancel_scope.cancel)          # runs 2nd
stack.callback(manager.begin_shutdown)             # runs 1st: stop flag for the thread
group.start_soon(manager.run_worker)
await manager.auto_resume([])                      # re-queue unfinished jobs under the roots
```

Releasing the locks before the thread has exited would let a job keep a saved-mapping store's
lock after its owner is gone. The pre-flight doc records the same order for the server.

**Keep the heavy import out.** Only `biomapper/mcp/server.py` imports the MCP SDK. The facade
imports the runner modules, not the server, and a subprocess test asserts that
`import biomapper.jobs` leaves no `mcp` module in `sys.modules`.

Roots are passed once at construction as the config's `env_roots`, and every call passes
`client_roots=[]`. That keeps the runner's path checks (inputs, `output_dir`, `store_dir` must be
inside a root, never through a link) in force for a web backend.

### 2. For progress with an injected mapper, the caller has to batch.

Both functions take `mapper=`, any callable that maps a list of requests. Batching happened
inside the default mapper, where an injected mapper or a subclass that overrides `__call__` never
sees a hook. So progress can only be reported from the calling side, by sending the mapper one
batch at a time:

```python
step = total if on_batch is None else max(batch_size, 1)
for start in range(0, total, step):
    got = list(mapper(requests[start : start + step], ...))
    results.extend(got)
    if on_batch is not None:
        on_batch(len(results), total)
```

Do this only when a callback is given. Without one the mapper still gets every request in one
call, so existing callers and tests that count mapper calls are unaffected. The default mapper
splits into the same batches itself, so what reaches the API is identical either way; a test
asserts `loud.sent == quiet.sent` while the call counts differ.

Two more points:

- **Say what the counts count.** `harmonize_cohorts()` maps one cohort at a time, so an update
  names the cohort. `harmonize_many()` sends each distinct entity once for all cohorts, so
  per-cohort progress does not exist there: `cohort` is `None`, `cohorts` lists every label, and
  the counts are distinct entities per arm. Reused saved mappings count as done from the first
  update, which is why a resumed run starts at, say, 5 of 10.
- **A progress callback must not cost the run.** Mapping is the expensive part. A callback that
  raises is reported once with `warnings.warn`, is not called again, and the run continues.

In a job the function runs on a worker thread, so the runner marshals each update to the event
loop with `anyio.from_thread.run_sync`. The listener then runs on the caller's loop and needs no
locking. The latest update is also kept on the job record, so a caller that polls status sees it.

### 3. A restart test must separate "saved" from "in flight".

The first version of the restart test waited until five names had been sent, exited, reopened,
and asserted every entity was sent exactly once. It failed with five extra sends. That was not a
bug: the batch in flight when the process died had been sent, its answer never arrived, and it
was correctly sent again.

The existing MCP shutdown test did not catch this because it compared `sorted(set(sent_names))`,
which cannot see a repeat.

The guarantee worth testing is narrower and stronger: no entity whose answer was saved is sent
again.

```python
saved = (await jobs.status(job_id)).entities_mapped   # wait until this is non-zero, then exit
...
for name in sent[:saved]:           # batches go out in order, one at a time
    assert sent.count(name) == 1
assert set(sent) == unique
assert len(sent) <= len(unique) + chunk_size          # only the in-flight batch may repeat
```

State the same thing in the docs: "at most the one batch in flight at the kill is sent twice".
"Without re-mapping" on its own promises more than a checkpointing runner can deliver.

### 4. A job that survives a restart needs its result on disk as data.

After a restart there is no report object to hand back. Every report class got a `to_dict()`
that is JSON-safe with no absolute paths, and `write()` saves it as `report.json` next to the
tables, under a `schema_version`. Keys may be added within a version; a rename or removal raises
it. Tests pin the key sets.

Adding an output file has a consequence in the runner: `REPORT_PATTERNS` in
`harmonize_jobs.py` lists the names an existing `output_dir` may not already hold, and it must
include the new file. This is the "no-overwrite for every artifact" rule from the file-jobs
review-gap doc applied to a new artifact.

Where a table is streamed to disk to save memory (the cross-node row pairs of a saved
`harmonize_many` run), the JSON holds `null` and names the file. Reading it back into memory to
fill the JSON would undo the reason it was streamed.

### 5. When behaviour is deliberate, report it. Do not add a switch.

Name linking is derived from the category (small molecules only), because gene and protein
symbols are case-sensitive. The request offered a choice between an override and reporting.
Reporting was chosen: a `name_linking` block (`enabled`, `category`, `reason`,
`overridable: false`) in the settings, the summary and `report.json`. The UI shows why, and no
linking rule changed.

### 6. Release notes for this repo

- **A version bump breaks `tests/test_version.py` until the install metadata is refreshed.** It
  compares the installed distribution's version with `pyproject.toml`. Run `poetry install` after
  editing the version, before running the suite.
- **PyPI's unversioned JSON lags a fresh upload.** Minutes after publishing,
  `https://pypi.org/pypi/biomapper/json` still said 1.8.0 while
  `https://pypi.org/pypi/biomapper/1.9.0/json` and `pip install biomapper==1.9.0` were correct.
  Verify with the versioned URL and a pinned install in a clean environment.
- **The version bump rode in the feature PR** as its own last commit, which cost one Greptile
  review and one merge instead of two. The fork merges with a merge commit, so the commit
  survives in history.
- **Order used:** fork PR merged, build and publish from fork `main`, tag on the fork, then the
  `sync/fork-main-<sha>` PR to the org, then the tag on the org. Earlier releases published after
  the org sync. The artifact is the same commit either way; publishing first only avoids waiting
  on a merge the agent is not allowed to make.
- **A migration note is logic, not prose.** Greptile's one finding on the PR was in
  `docs/migration-1.9.0.md`: it claimed `medium` could only shrink and the mean could only fall.
  When rows can drop a tier or leave the scored set, middle tiers and averages can move either
  way. Only per-row statements ("equal, lower or blank") and the end tiers are one-directional.

## Why This Matters

- A second wrapper would have drifted from the MCP one within a release. The runner had already
  absorbed two rounds of lifecycle fixes (cancel races, placement recovery, store privacy) that a
  rewrite would have had to rediscover.
- A facade that releases locks in the wrong order leaves a store locked by a dead owner, and the
  next run is refused.
- A resume test that asserts exactly-once either fails for a correct runner or, written against
  a set, passes for a broken one.
- Progress added by changing the default call pattern would have broken every replay test that
  counts mapper calls, for a feature most callers do not use.
- A downstream UI that reparses TSV tables, or reimplements the job wrapper, owns a copy of logic
  that belongs to the package.

## When to Apply

- Before exposing any internal runner, service or manager to a second consumer: look for the
  lifespan that already wires it, and copy its order exactly.
- Whenever a hook has to observe work done inside a caller-supplied callable.
- Whenever a test claims "no repeated work after a crash": find what was durable at the moment
  of the kill and assert on that set.
- Whenever a long-running job's result is needed after the process that produced it is gone.

## Examples

Using the result from a FastAPI lifespan:

```python
from biomapper.jobs import CohortSpec, HarmonizeJobs

async with HarmonizeJobs(home="/var/lib/app/biomapper", roots=["/var/lib/app/runs"],
                         on_progress=lambda job_id, p: publish(job_id, p.to_dict())) as jobs:
    job = await jobs.submit_many(
        [CohortSpec(label="a", path=".../a.tsv"), CohortSpec(label="b", path=".../b.tsv"),
         CohortSpec(label="c", path=".../c.tsv")],
        entity="metabolites", output_dir=".../report", store_dir=".../store",
    )
    done = await jobs.wait(job.job_id)
    report = jobs.load_report(job.job_id)
```

A live check on 2026-11-13 against the public deployment (three invented cohorts, 10 distinct
names, `chunk_size` 5) killed the service after the first saved batch and reopened it. Progress
resumed at 5 of 10 and the job completed. This was one small run, not a benchmark.

Tests that carry these lessons:

- `tests/mcp_server/test_public_jobs.py`: restart for two-cohort and many-cohort jobs, progress
  on the event loop, a failing listener, the no-SDK-import check.
- `tests/test_harmonize_progress.py`: the exact update sequence per arm and cohort, and the call
  pattern with and without a callback.
- `tests/test_harmonize_plain_data.py`: pinned key sets for every `to_dict()` and `report.json`.

## Related

- `docs/solutions/best-practices/mcp-harmonize-jobs-preflight-2026-11-11.md`: the decisions
  behind the MCP harmonize job, including the shutdown order this facade mirrors. Its note that
  jobs run with `progress=False` and only stdout messages are routed is still true, and is now
  joined by `on_progress`.
- `docs/solutions/best-practices/mcp-server-file-jobs-review-gap-2026-11-03.md`: the job
  lifecycle failures found in review, including "no-overwrite for every artifact".
- `docs/solutions/best-practices/mcp-registry-pypi-launch-contract-2026-11-04.md`: the other
  release-process doc for this package.
- `docs/migration-1.9.0.md`: the `confidence_score` change released alongside.
- Phenome-Health/biomapper#16: a three-or-more cohort example with a recorded replay, left out of
  1.9.0 because no third cohort's licence could be confirmed.
