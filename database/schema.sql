-- KlebTrace database schema (SQLite).
--
-- Every table is keyed by assembly accession, so annotations, resistance calls
-- and laboratory results for one genome can be joined in one query. Standard SQL
-- only (apart from the SQLite pragma), so the same schema can move to
-- PostgreSQL or MySQL.

PRAGMA foreign_keys = ON;

-- One row per study genome: why it was chosen and the assembly checks from
-- download_genomes.py.
CREATE TABLE genomes (
    assembly              TEXT PRIMARY KEY,        -- e.g. GCA_000764615.1
    biosample             TEXT NOT NULL,
    country               TEXT,
    collection_year       INTEGER,
    region                TEXT,
    selection_group       TEXT NOT NULL,
    carbapenemase_family  TEXT,
    in_ensembl            INTEGER CHECK (in_ensembl IN (0, 1)),
    n_sequences           INTEGER,
    total_bp              INTEGER,
    n50                   INTEGER,
    gc_percent            REAL,
    md5                   TEXT
);

-- Sequences (contigs, chromosomes, plasmids) in each FASTA. Features must be
-- on one of these, which is how mismatched contig names are caught.
CREATE TABLE sequences (
    assembly  TEXT NOT NULL REFERENCES genomes (assembly),
    seqid     TEXT NOT NULL,                       -- accession, e.g. CP008933.1
    length    INTEGER NOT NULL CHECK (length > 0),
    PRIMARY KEY (assembly, seqid)
);

-- Laboratory meropenem result from the AMR Portal.
CREATE TABLE phenotypes (
    assembly          TEXT NOT NULL REFERENCES genomes (assembly),
    antibiotic        TEXT NOT NULL,
    phenotype         TEXT NOT NULL CHECK (phenotype IN ('resistant', 'susceptible')),
    mic_mg_l          TEXT,                        -- with sign, e.g. '>32'
    strict_resistant  INTEGER CHECK (strict_resistant IN (0, 1)),
    PRIMARY KEY (assembly, antibiotic)
);

-- Features from the merged GFF3 (merge_annotations.py), one row per GFF3 line.
CREATE TABLE features (
    feature_id  INTEGER PRIMARY KEY,
    assembly    TEXT NOT NULL,
    seqid       TEXT NOT NULL,
    source      TEXT NOT NULL,                     -- tool, e.g. Bakta, AMRFinderPlus
    type        TEXT NOT NULL,                     -- Sequence Ontology term
    start       INTEGER NOT NULL CHECK (start >= 1),
    end         INTEGER NOT NULL,
    strand      TEXT NOT NULL CHECK (strand IN ('+', '-', '.', '?')),
    gff_id      TEXT NOT NULL,                     -- the ID attribute
    name        TEXT,
    product     TEXT,
    attributes  TEXT,                              -- the full column 9
    CHECK (end >= start),
    UNIQUE (assembly, gff_id),
    FOREIGN KEY (assembly, seqid) REFERENCES sequences (assembly, seqid)
);

-- Resistance calls from every source, in one shape so they can be compared.
-- The portal's per-antibiotic repeat rows are collapsed before loading.
CREATE TABLE resistance_calls (
    call_id          INTEGER PRIMARY KEY,
    assembly         TEXT NOT NULL REFERENCES genomes (assembly),
    source           TEXT NOT NULL CHECK (source IN ('portal', 'amrfinderplus', 'kleborate', 'split_gene_screen')),
    gene             TEXT NOT NULL,                -- e.g. blaKPC-2, ompK36_D135DGD
    element_subtype  TEXT,                         -- AMR, POINT, POINT_DISRUPT
    subclass         TEXT,
    is_carbapenemase INTEGER NOT NULL CHECK (is_carbapenemase IN (0, 1)),
    method           TEXT,                         -- e.g. ALLELEX, PARTIAL_CONTIG_ENDX, BLASTP
    is_partial       INTEGER NOT NULL DEFAULT 0 CHECK (is_partial IN (0, 1)),
    seqid            TEXT,
    start            INTEGER,
    end              INTEGER,
    coverage_pct     REAL,
    identity_pct     REAL,
    closest_reference TEXT
);

-- Tool and database versions behind each source (from galaxy/tool_versions.tsv
-- and the pipeline), so every call can be traced to a version.
CREATE TABLE tool_versions (
    source            TEXT PRIMARY KEY,
    tool              TEXT NOT NULL,
    tool_version      TEXT NOT NULL,
    database_version  TEXT,
    run_where         TEXT                         -- e.g. usegalaxy.org, Nextflow
);

CREATE INDEX features_by_location ON features (assembly, seqid, start, end);
CREATE INDEX calls_by_gene ON resistance_calls (assembly, gene);

-- Carbapenem mechanisms per genome and source, next to the laboratory result.
CREATE VIEW carbapenem_summary AS
SELECT g.assembly,
       g.selection_group,
       p.phenotype,
       p.mic_mg_l,
       c.source,
       group_concat(CASE WHEN c.is_carbapenemase = 1 THEN c.gene END, ';') AS carbapenemases,
       group_concat(CASE WHEN c.element_subtype LIKE 'POINT%' AND c.subclass LIKE '%CARBAPENEM%'
                         THEN c.gene END, ';') AS porin_changes
FROM genomes g
LEFT JOIN phenotypes p ON p.assembly = g.assembly AND p.antibiotic = 'meropenem'
LEFT JOIN resistance_calls c ON c.assembly = g.assembly
GROUP BY g.assembly, c.source;
