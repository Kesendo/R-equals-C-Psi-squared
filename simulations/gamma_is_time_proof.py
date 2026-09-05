#!/usr/bin/env python3
"""Withdrawal-safe tombstone for the former gamma/time argument."""

from pathlib import Path


OUTPUT_PATH = Path(__file__).parent / "results" / "gamma_is_time_proof.txt"


def current_report():
    return (
        "WITHDRAWN: this script previously promoted finite two-qubit trajectory\n"
        "comparisons into an unsupported ontology. Those calculations do not\n"
        "establish a necessary or sufficient condition for experience.\n"
        "\n"
        "Current scoped account: docs/GAMMA_TIME_DISTINCTION.md\n"
        "Current computational reading: gamma sets a dissipative scale in the\n"
        "tested Lindblad models; coherent motion and other dimensionless ratios\n"
        "prevent a universal collapse onto tau=gamma*t.\n"
    )


def main():
    report = current_report()
    OUTPUT_PATH.write_text(report, encoding="utf-8")
    print(report, end="")


if __name__ == "__main__":
    main()
