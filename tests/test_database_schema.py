"""Tests for database/schema.sql and the database creation step."""

from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import load_database as ld  # noqa: E402

SCHEMA = ROOT / "database" / "schema.sql"


@pytest.fixture
def con(tmp_path):
    con = ld.create(tmp_path / "k.sqlite", SCHEMA)
    con.execute("INSERT INTO genomes (assembly, biosample, selection_group) VALUES ('G.1', 'S1', 'susceptible')")
    con.execute("INSERT INTO sequences VALUES ('G.1', 'CP1.1', 1000)")
    return con


def feature(con, gff_id="f1", seqid="CP1.1", start=1, end=100, strand="+"):
    con.execute(
        "INSERT INTO features (assembly, seqid, source, type, start, end, strand, gff_id) "
        "VALUES ('G.1', ?, 'Bakta', 'CDS', ?, ?, ?, ?)", (seqid, start, end, strand, gff_id))


def test_schema_creates_all_tables(con):
    names = {r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type IN ('table', 'view')")}
    assert {"genomes", "sequences", "phenotypes", "features", "resistance_calls",
            "tool_versions", "carbapenem_summary"} <= names


@pytest.mark.parametrize("kwargs", [
    {"seqid": "contig_1"},          # contig name not in the FASTA
    {"start": 0},                   # GFF3 counts from 1
    {"start": 50, "end": 10},       # start after end
    {"strand": ""},                 # empty strand, as in the Phase 0 ISEScan GFF
])
def test_bad_features_are_rejected(con, kwargs):
    with pytest.raises(sqlite3.IntegrityError):
        feature(con, **kwargs)


def test_duplicate_feature_id_is_rejected(con):
    feature(con)
    with pytest.raises(sqlite3.IntegrityError):
        feature(con)
