"""CLI for the external benchmark suite.

``--endpoint`` defaults to **production**, because a benchmark that measures a dev checkout is
not measuring the service a reader can call. ``dev`` is available for future testing.

The API key is read from ``BIOMAPPER_API_KEY`` and never taken on the command line: argv is
visible to every process on the host and lands in shell history. A deployment with no keys
configured is open, and the suite runs against it unauthenticated without needing a placeholder.

There is no ``--no-save``. The expensive part of a run is live API traffic, and a flag that
discards it is not an acceptable failure mode; ``--out`` overrides *where* results land, never
*whether* they do.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from pathlib import Path

from biomapper.benchmarks.config import SUITE_DATASETS, SUITE_SKIPPED
from biomapper.benchmarks.provenance import DEFAULT_KESTREL_URL
from biomapper.benchmarks.suite import ENDPOINTS, run_suite


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m biomapper.benchmarks",
        description="Run the BioMapper external benchmark suite against a deployment (API only).",
    )
    parser.add_argument(
        "--endpoint",
        default="production",
        help=(
            f"API endpoint: one of {sorted(ENDPOINTS)} or a full http(s) URL. "
            f"Default: production."
        ),
    )
    parser.add_argument(
        "--kestrel-url",
        default=DEFAULT_KESTREL_URL,
        help=(
            "Kestrel base URL, read for run provenance (/health) and for node-name lookups in the "
            "structure oracle's name-fallback path. The public host is keyless and is never sent a "
            f"key. Default: {DEFAULT_KESTREL_URL}"
        ),
    )
    parser.add_argument(
        "--out", default=None, help="Override the output dir (default: timestamped)."
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=20,
        help="Entities per /map/batch request. Default: 20.",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Debug logging.")

    sub = parser.add_subparsers(dest="command", required=True)

    all_parser = sub.add_parser("all", help=f"Run all {len(SUITE_DATASETS)} arms.")
    all_parser.add_argument(
        "--only",
        nargs="+",
        default=None,
        choices=SUITE_DATASETS,
        help="Restrict to these arms (still writes one suite manifest).",
    )

    arm_parser = sub.add_parser("arm", help="Run a single arm.")
    arm_parser.add_argument("name", choices=SUITE_DATASETS, help="Arm to run.")

    sub.add_parser("list", help="List the arms and the deliberate skips, then exit.")

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    if args.command == "list":
        print(f"Suite arms ({len(SUITE_DATASETS)}):")
        for key in SUITE_DATASETS:
            print(f"  {key}")
        print(f"\nDeliberate skips ({len(SUITE_SKIPPED)}), recorded in every manifest:")
        for key, reason in SUITE_SKIPPED.items():
            print(f"  {key}: {reason}")
        return 0

    datasets = [args.name] if args.command == "arm" else args.only

    # Read from the environment only. See the module docstring on why not from argv.
    api_key = os.getenv("BIOMAPPER_API_KEY")

    outcome = run_suite(
        out_dir=Path(args.out) if args.out else None,
        datasets=datasets,
        endpoint=args.endpoint,
        api_key=api_key,
        kestrel_url=args.kestrel_url,
        batch_size=args.batch_size,
    )
    manifest = outcome["manifest"]
    print(json.dumps({k: v for k, v in manifest.items() if k != "datasets"}, indent=2, default=str))
    print(f"\nResults saved to: {outcome['out_dir']}")
    for entry in manifest["datasets"]:
        note = entry.get("reason") or entry.get("error") or ""
        print(f"  {entry['status']:8s} {entry['dataset']:24s} {note}")

    # A failed arm is a non-zero exit so a scheduled run is not reported as green. A SKIP is not a
    # failure: it is a recorded, deliberate outcome with a reason attached.
    return 1 if manifest["n_failed"] else 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
