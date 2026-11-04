"""Check that every KlebTrace script called from a Nextflow module gets all its required arguments.

`-stub-run` never runs the real commands, so a module that forgets an argument
only fails at the end of a real run (issue #1: LOAD_DATABASE without --fasta).
This test reads each module's command and compares it with the script's own
required arguments from `--help`.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULES = sorted((ROOT / "nextflow" / "modules").glob("*.nf"))


def required_flags(usage: str) -> set[str]:
    """Flags outside square brackets in an argparse usage line are required."""
    depth, flags = 0, set()
    for token in re.findall(r"\[|\]|--[\w-]+", usage):
        if token == "[":
            depth += 1
        elif token == "]":
            depth -= 1
        elif depth == 0:
            flags.add(token)
    return flags


def script_calls(module_text: str) -> dict[str, str]:
    """Return {script name: command text} for each KlebTrace script in the module's script: blocks."""
    calls = {}
    for block in re.findall(r'script:\s*"""(.*?)"""', module_text, re.S):
        for name in re.findall(r"\b(\w+\.py)\b", block):
            if (ROOT / "scripts" / name).exists():   # skip the tools' own scripts, e.g. plasmidfinder.py
                calls[name] = block
    return calls


def usage_of(script: str) -> str:
    out = subprocess.run([sys.executable, str(ROOT / "scripts" / script), "--help"],
                         capture_output=True, text=True, check=True).stdout
    return out.split("\n\n")[0]


def missing(command: str, usage: str) -> set[str]:
    used = set(re.findall(r"--[\w-]+", command))
    return required_flags(usage) - used


CALLS = [(m.name, s, c) for m in MODULES for s, c in script_calls(m.read_text()).items()]


def test_modules_call_scripts():
    assert {s for _, s, _ in CALLS} >= {"merge_annotations.py", "load_database.py", "qc_checks.py"}


@pytest.mark.parametrize("module, script, command", CALLS, ids=[f"{m}:{s}" for m, s, _ in CALLS])
def test_module_passes_required_arguments(module, script, command):
    assert not missing(command, usage_of(script)), f"{module} calls {script} without {missing(command, usage_of(script))}"


def test_issue_1_command_would_fail():
    # The LOAD_DATABASE command as it was before the fix
    broken = ("load_database.py --schema $schema --samples $samples --portal-calls $portal_calls "
              "--gff3 $gff3s --amrfinderplus $amrfinder_reports --kleborate $kleborate_reports --out klebtrace.sqlite")
    assert missing(broken, usage_of("load_database.py")) == {"--fasta"}
