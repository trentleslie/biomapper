"""Request pacing (offline; time.sleep is captured, never actually slept)."""

from __future__ import annotations

import pytest

from biomapper.benchmarks.pacing import PUBCHEM_MIN_INTERVAL_S, Pacer


@pytest.fixture
def sleeps(monkeypatch):
    recorded: list[float] = []
    monkeypatch.setattr("time.sleep", lambda s: recorded.append(s))
    return recorded


def test_the_first_call_does_not_sleep(sleeps):
    Pacer(0.25).wait()
    assert sleeps == []


def test_every_subsequent_call_is_paced(sleeps):
    # The defect this replaces slept once per ROW while issuing two REQUESTS, so a fallback lookup
    # went out unspaced and could draw a lookup_failed of its own.
    pacer = Pacer(0.25)
    for _ in range(3):
        pacer.wait()
    assert len(sleeps) == 2
    assert all(0 < s <= 0.25 for s in sleeps)


def test_a_disabled_pacer_never_sleeps(sleeps):
    pacer = Pacer(0)
    for _ in range(5):
        pacer.wait()
    assert sleeps == []


def test_a_slow_caller_is_not_penalised(monkeypatch, sleeps):
    # If more than the interval already elapsed, there is nothing to wait for.
    clock = iter([0.0, 0.0, 10.0, 10.0])
    monkeypatch.setattr("time.monotonic", lambda: next(clock))
    pacer = Pacer(0.25)
    pacer.wait()
    pacer.wait()
    assert sleeps == []


def test_the_default_respects_pubchems_published_ceiling():
    # PUG-REST publishes 5 requests per second.
    assert PUBCHEM_MIN_INTERVAL_S >= 1 / 5
    assert Pacer().min_interval_s == PUBCHEM_MIN_INTERVAL_S
