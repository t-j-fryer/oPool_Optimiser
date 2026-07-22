"""Reusable end-to-end oPool cloning workflow.

This module powers both the terminal CLI and the streamlined notebook. It runs:

1. amino-acid reverse translation and DNA Chisel optimization (or accepts an
   existing optimized CSV),
2. fast fragment/pool assignment with native or synonymous 4-nt overhangs,
3. primer, Type IIS, vector-overhang, and optional stuffer assembly.
"""

from __future__ import annotations

import bisect
import math
import time
from dataclasses import dataclass, field
from itertools import combinations, product
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from Bio.Data import CodonTable
from Bio.Seq import Seq
from dnachisel import (
    AvoidPattern,
    CodonOptimize,
    DnaOptimizationProblem,
    EnforceGCContent,
    EnforceTranslation,
)
from python_codon_tables import available_codon_tables_names


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data"
DNA_BASES = set("ACGT")
BUILTIN_CODON_SPECIES = tuple(sorted({
    "_".join(table_name.split("_")[:-1])
    for table_name in available_codon_tables_names
}))
BUILTIN_CODON_TABLES = frozenset(available_codon_tables_names)


def _path(value: str | Path | None) -> Path | None:
    return None if value is None else Path(value).expanduser()


def _normalize_dna(value: str, label: str, exact_length: int | None = None) -> str:
    sequence = str(value).strip().upper()
    if not sequence or set(sequence) - DNA_BASES:
        raise ValueError(f"{label} must contain only A/C/G/T and cannot be empty.")
    if exact_length is not None and len(sequence) != exact_length:
        raise ValueError(f"{label} must be exactly {exact_length} nucleotide(s).")
    return sequence


def reverse_complement(sequence: str) -> str:
    return str(Seq(sequence).reverse_complement())


@dataclass
class WorkflowConfig:
    """User-facing workflow settings with conservative defaults."""

    input_path: str | Path
    input_kind: str = "auto"  # auto, aa, optimized
    output_dir: str | Path | None = None
    run_name: str | None = None
    overhangs_path: str | Path | None = None
    primers_path: str | Path | None = None

    opool_length: int = 250
    genes_per_subpool: int | None = None
    vector_oh1: str = "GCTT"
    vector_oh2: str = "AGTG"
    short_pool_max_size: int | None = 1000
    max_fragments: int = 32
    strip_nterm_met: bool = True

    typeiis_site: str = "GGTCTC"
    typeiis_n: str = "A"
    additional_avoid_sequences: tuple[str, ...] = (
        "AAAAA",
        "GGGGG",
        "CCCCC",
        "TTTTT",
    )
    gc_min: float = 0.30
    gc_max: float = 0.70
    gc_window: int = 50
    codon_species: str = "e_coli"
    codon_method: str = "match_codon_usage"

    primer_mode: str = "combinatorial"  # combinatorial, unique_pairs
    primer_start_at: str | None = None
    add_stuffer: bool = True
    stuffer_gc_min: float = 0.40
    stuffer_gc_max: float = 0.60
    stuffer_max_homopolymer: int = 4
    stuffer_max_tries: int = 10_000
    random_seed: int | None = 123
    skip_primer_assignment: bool = False
    exclude_unsafe_primers: bool = True

    show_progress: bool = True
    force: bool = False

    def __post_init__(self) -> None:
        self.input_path = _path(self.input_path)
        self.output_dir = _path(self.output_dir)
        self.overhangs_path = _path(self.overhangs_path)
        self.primers_path = _path(self.primers_path)
        self.input_kind = str(self.input_kind).strip().lower()
        self.primer_mode = str(self.primer_mode).strip().lower()
        self.codon_species = str(self.codon_species).strip().lower()
        self.typeiis_site = _normalize_dna(self.typeiis_site, "typeiis_site")
        self.typeiis_n = _normalize_dna(self.typeiis_n, "typeiis_n", exact_length=1)
        self.vector_oh1 = _normalize_dna(self.vector_oh1, "vector_oh1", exact_length=4)
        self.vector_oh2 = _normalize_dna(self.vector_oh2, "vector_oh2", exact_length=4)
        self.additional_avoid_sequences = tuple(
            _normalize_dna(sequence, "additional_avoid_sequences entry")
            for sequence in self.additional_avoid_sequences
        )
        self._validate()

    @property
    def forbidden_typeiis_sites(self) -> set[str]:
        return {self.typeiis_site, reverse_complement(self.typeiis_site)}

    @property
    def avoid_sequences(self) -> list[str]:
        return list(dict.fromkeys([
            *sorted(self.forbidden_typeiis_sites),
            *self.additional_avoid_sequences,
        ]))

    def _validate(self) -> None:
        if self.input_kind not in {"auto", "aa", "optimized"}:
            raise ValueError("input_kind must be 'auto', 'aa', or 'optimized'.")
        if self.primer_mode not in {"combinatorial", "unique_pairs"}:
            raise ValueError("primer_mode must be 'combinatorial' or 'unique_pairs'.")
        if not self.codon_species:
            raise ValueError("codon_species cannot be empty.")
        if (
            not self.codon_species.isdigit()
            and self.codon_species not in BUILTIN_CODON_SPECIES
            and self.codon_species not in BUILTIN_CODON_TABLES
        ):
            choices = ", ".join(BUILTIN_CODON_SPECIES)
            raise ValueError(
                f"Unknown codon_species {self.codon_species!r}. "
                f"Use a built-in keyword ({choices}) or an NCBI taxonomy ID."
            )
        if not isinstance(self.opool_length, int) or self.opool_length <= 62:
            raise ValueError("opool_length must be an integer greater than 62.")
        for label, value in (
            ("genes_per_subpool", self.genes_per_subpool),
            ("short_pool_max_size", self.short_pool_max_size),
        ):
            if value is not None and (not isinstance(value, int) or value < 1):
                raise ValueError(f"{label} must be None or a positive integer.")
        if not isinstance(self.max_fragments, int) or self.max_fragments < 2:
            raise ValueError("max_fragments must be an integer of at least 2.")
        if not 0 <= self.gc_min <= self.gc_max <= 1:
            raise ValueError("Require 0 <= gc_min <= gc_max <= 1.")
        if not 0 <= self.stuffer_gc_min <= self.stuffer_gc_max <= 1:
            raise ValueError("Require 0 <= stuffer_gc_min <= stuffer_gc_max <= 1.")
        if self.gc_window < 1:
            raise ValueError("gc_window must be positive.")
        if len(self.typeiis_site) < 5:
            raise ValueError("typeiis_site must be at least five nucleotides long.")


@dataclass(frozen=True)
class RunPaths:
    run_name: str
    run_dir: Path
    optimized: Path
    assigned: Path
    strip_log: Path
    unassigned: Path
    overhangs_used: Path
    full_info: Path
    references_fasta: Path
    fragments: Path
    all_pairs: Path
    unused_pairs: Path
    unused_primers: Path

    @classmethod
    def from_config(cls, config: WorkflowConfig) -> "RunPaths":
        input_path = config.input_path.resolve()
        run_dir = (
            config.output_dir.resolve()
            if config.output_dir is not None
            else PROJECT_ROOT / "outputs"
            if input_path.parent.resolve() == DATA_DIR.resolve()
            else input_path.parent
        )
        if config.run_name and config.run_name.strip():
            run_name = config.run_name.strip()
        elif input_path.parent.resolve() == DATA_DIR.resolve():
            run_name = input_path.stem
        else:
            run_name = input_path.parent.name
        return cls(
            run_name=run_name,
            run_dir=run_dir,
            optimized=run_dir / f"{run_name}_Optimised.csv",
            assigned=run_dir / f"{run_name}_Assigned.csv",
            strip_log=run_dir / f"{run_name}_stripped_ATG_log.csv",
            unassigned=run_dir / f"{run_name}_unassigned.csv",
            overhangs_used=run_dir / f"{run_name}_overhangs_used.csv",
            full_info=run_dir / f"{run_name}_FULL_INFO.csv",
            references_fasta=run_dir / f"{run_name}_references.fasta",
            fragments=run_dir / f"{run_name}_oPool_Order_Fragments.csv",
            all_pairs=run_dir / f"{run_name}_orthogonal_oligos_pairs_all.csv",
            unused_pairs=run_dir / f"{run_name}_orthogonal_oligos_pairs_unused.csv",
            unused_primers=run_dir / f"{run_name}_orthogonal_oligos_unused.csv",
        )

    def as_dict(self) -> dict[str, Path | str]:
        return {
            field_name: getattr(self, field_name)
            for field_name in self.__dataclass_fields__
        }


@dataclass(frozen=True)
class GeneRecord:
    name: str
    sequence: str


@dataclass
class WorkflowResult:
    paths: RunPaths
    optimized_genes: int
    assigned_genes: int
    unassigned_genes: int
    blocks: int
    runtime_seconds: float
    output_files: list[Path] = field(default_factory=list)


def _resolve_inventory(path: Path | None, default_name: str) -> Path:
    if path is None:
        result = DATA_DIR / default_name
    elif path.is_absolute():
        result = path
    elif len(path.parts) == 1:
        result = DATA_DIR / path
    else:
        result = PROJECT_ROOT / path
    result = result.resolve()
    if not result.is_file():
        raise FileNotFoundError(f"Inventory file not found: {result}")
    return result


def _detect_input_kind(path: Path, requested: str) -> str:
    if requested != "auto":
        return requested
    header = pd.read_csv(path, nrows=0)
    columns = {str(column).strip().lower() for column in header.columns}
    if "dna_seq_optimized" in columns:
        return "optimized"
    if {"name", "aa_seq"}.issubset(columns):
        return "aa"
    return "aa"


def _read_aa_input(path: Path) -> pd.DataFrame:
    header = pd.read_csv(path, nrows=0)
    normalized = {str(column).strip().lower(): column for column in header.columns}
    if {"name", "aa_seq"}.issubset(normalized):
        df = pd.read_csv(path, dtype={normalized["name"]: str})
        df = df.rename(columns={normalized["name"]: "name", normalized["aa_seq"]: "aa_seq"})
        df = df[["name", "aa_seq"]]
    else:
        df = pd.read_csv(path, header=None, names=["name", "aa_seq"], dtype={"name": str})
    df["name"] = df["name"].astype(str).str.strip()
    df["aa_seq"] = df["aa_seq"].astype(str).str.replace(r"\s+", "", regex=True).str.upper()
    if (df["name"] == "").any() or (df["aa_seq"] == "").any():
        raise ValueError("Input contains an empty sequence name or amino-acid sequence.")
    if df["name"].duplicated().any():
        duplicates = df.loc[df["name"].duplicated(), "name"].tolist()[:5]
        raise ValueError(f"Sequence names must be unique. Duplicates include: {duplicates}")
    valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY*")
    bad = [name for name, sequence in zip(df["name"], df["aa_seq"]) if set(sequence) - valid_amino_acids]
    if bad:
        raise ValueError(f"Invalid amino-acid character(s) in sequences including: {bad[:5]}")
    return df


def _read_optimized_input(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype={"name": str})
    required = {"name", "dna_seq_optimized"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Optimized input is missing columns: {sorted(missing)}")
    if "aa_seq" not in df.columns:
        df["aa_seq"] = ""
    df["name"] = df["name"].astype(str).str.strip()
    df["dna_seq_optimized"] = (
        df["dna_seq_optimized"].astype(str).str.replace(r"\s+", "", regex=True).str.upper()
    )
    for name, sequence in zip(df["name"], df["dna_seq_optimized"]):
        if not sequence or set(sequence) - DNA_BASES or len(sequence) % 3:
            raise ValueError(f"Optimized DNA for {name!r} must be A/C/G/T and divisible by 3.")
    if df["name"].duplicated().any():
        raise ValueError("Optimized input contains duplicate sequence names.")
    return df[["name", "aa_seq", "dna_seq_optimized"]]


def reverse_translate(protein_sequence: str) -> str:
    table = CodonTable.unambiguous_dna_by_id[1]
    aa_to_codon: dict[str, str] = {}
    for codon, aa in table.forward_table.items():
        aa_to_codon.setdefault(aa, codon.upper())
    aa_to_codon["*"] = table.stop_codons[0].upper()
    return "".join(aa_to_codon[aa] for aa in protein_sequence.upper())


def optimize_dna(dna_sequence: str, config: WorkflowConfig) -> str:
    sequence_length = len(dna_sequence)
    problem = DnaOptimizationProblem(
        sequence=dna_sequence,
        logger=None,
        constraints=[
            *(AvoidPattern(pattern) for pattern in config.avoid_sequences),
            EnforceGCContent(
                mini=config.gc_min,
                maxi=config.gc_max,
                window=config.gc_window,
            ),
            EnforceTranslation(location=(0, sequence_length)),
        ],
        objectives=[
            CodonOptimize(
                species=config.codon_species,
                method=config.codon_method,
                location=(0, sequence_length),
            )
        ],
    )
    problem.resolve_constraints()
    problem.optimize()
    return str(problem.sequence).upper()


def prepare_optimized_input(config: WorkflowConfig, paths: RunPaths) -> pd.DataFrame:
    input_path = config.input_path.resolve()
    if not input_path.is_file():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    input_kind = _detect_input_kind(input_path, config.input_kind)
    if input_kind == "optimized":
        optimized = _read_optimized_input(input_path)
        if input_path != paths.optimized.resolve():
            optimized.to_csv(paths.optimized, index=False)
        if config.show_progress:
            print(f"[1/3] Reusing {len(optimized)} optimized DNA sequences from {input_path}")
        return optimized

    aa_df = _read_aa_input(input_path)
    if config.show_progress:
        print(f"[1/3] Optimizing {len(aa_df)} amino-acid sequences")
    optimized_sequences = []
    for index, amino_acid_sequence in enumerate(aa_df["aa_seq"], start=1):
        optimized_sequences.append(optimize_dna(reverse_translate(amino_acid_sequence), config))
        if config.show_progress:
            print(f"  optimized {index}/{len(aa_df)}", end="\r")
    if config.show_progress:
        print()
    result = aa_df.copy()
    result["dna_seq_optimized"] = optimized_sequences
    result.to_csv(paths.optimized, index=False)
    return result


class PoolAssigner:
    """Fast pool assignment using indexed cut options and dynamic path search."""

    def __init__(
        self,
        config: WorkflowConfig,
        paths: RunPaths,
        optimized_df: pd.DataFrame,
        overhangs_path: Path,
    ) -> None:
        self.config = config
        self.paths = paths
        self.optimized_df = optimized_df
        self.overhangs_path = overhangs_path
        self.overhangs = self._load_overhangs(overhangs_path)
        self.overhang_rank = {overhang: index for index, overhang in enumerate(self.overhangs)}
        self.overhang_bits = {overhang: 1 << index for index, overhang in enumerate(self.overhangs)}
        self.codon_table = CodonTable.unambiguous_dna_by_name["Standard"]
        self.synonyms: dict[str, list[str]] = {}
        for codon, aa in self.codon_table.forward_table.items():
            self.synonyms.setdefault(aa, []).append(codon)
        for aa in self.synonyms:
            self.synonyms[aa] = sorted(self.synonyms[aa])
        forbidden_internal = {config.vector_oh1, config.vector_oh2}
        self.internal_overhangs = [
            overhang for overhang in self.overhangs
            if overhang not in forbidden_internal
        ]
        if config.vector_oh2 in self.overhangs:
            raise ValueError("The overhang inventory contains vector_oh2; remove it.")
        if config.vector_oh1 in self.overhangs:
            raise ValueError("The overhang inventory contains vector_oh1; remove it.")

    @staticmethod
    def _load_overhangs(path: Path) -> list[str]:
        df = pd.read_csv(path, header=None)
        raw = str(df.iloc[0, 0])
        overhangs = [value.strip().upper() for value in raw.split(",") if value.strip()]
        if not overhangs:
            raise ValueError(f"No overhangs found in {path}")
        if len(overhangs) != len(set(overhangs)):
            raise ValueError("The overhang inventory contains duplicates.")
        if any(len(value) != 4 or set(value) - DNA_BASES for value in overhangs):
            raise ValueError("Every internal overhang must be four A/C/G/T bases.")
        return overhangs

    def _load_records(self) -> tuple[list[GeneRecord], list[dict[str, Any]]]:
        records: list[GeneRecord] = []
        strip_log: list[dict[str, Any]] = []
        for _, row in self.optimized_df.iterrows():
            name = str(row["name"]).strip()
            sequence = str(row["dna_seq_optimized"]).strip().upper()
            original_length = len(sequence)
            if self.config.strip_nterm_met and sequence.startswith("ATG"):
                sequence = sequence[3:]
                strip_log.append({
                    "Sequence Name": name,
                    "Action": "Removed leading ATG",
                    "Original_Length": original_length,
                    "New_Length": len(sequence),
                    "Original_Start": "ATG",
                    "New_Start": sequence[:12],
                })
            records.append(GeneRecord(name=name, sequence=sequence))
        return records, strip_log

    def _contains_forbidden_site(self, sequence: str) -> bool:
        return any(site in sequence for site in self.config.forbidden_typeiis_sites)

    @staticmethod
    def _mutated_slice(
        sequence: str,
        start: int,
        end: int,
        codon_start: int,
        new_codon: str,
    ) -> str:
        chars = list(sequence[start:end])
        for position in range(max(start, codon_start), min(end, codon_start + 3)):
            chars[position - start] = new_codon[position - codon_start]
        return "".join(chars)

    def _mutation_is_safe(self, sequence: str, codon_start: int, new_codon: str) -> bool:
        flank = max(len(site) for site in self.config.forbidden_typeiis_sites) - 1
        start = max(0, codon_start - flank)
        end = min(len(sequence), codon_start + 3 + flank)
        local = self._mutated_slice(sequence, start, end, codon_start, new_codon)
        return not any(site in local for site in self.config.forbidden_typeiis_sites)

    def _build_cut_index(self, sequence: str) -> tuple[dict[int, list[tuple]], list[int], str | None]:
        if self._contains_forbidden_site(sequence):
            return {}, [], "Input sequence already contains a forbidden Type IIS site"
        allowed = set(self.internal_overhangs)
        by_cut: dict[int, list[tuple]] = {}
        for cut in range(4, len(sequence) + 1):
            options: list[tuple] = []
            native = sequence[cut - 4:cut]
            if native in allowed:
                options.append((native, "Native", None, None))
            codon_starts = {
                position - (position % 3)
                for position in range(cut - 4, cut)
            }
            for codon_start in codon_starts:
                old_codon = sequence[codon_start:codon_start + 3]
                aa = self.codon_table.forward_table.get(old_codon)
                if aa is None:
                    continue
                for new_codon in self.synonyms.get(aa, []):
                    if new_codon == old_codon:
                        continue
                    mutated_window = self._mutated_slice(
                        sequence, cut - 4, cut, codon_start, new_codon
                    )
                    if mutated_window in allowed and self._mutation_is_safe(
                        sequence, codon_start, new_codon
                    ):
                        options.append((mutated_window, "Synonymous", codon_start, new_codon))
            if options:
                unique: list[tuple] = []
                seen = set()
                for option in options:
                    key = (option[0], option[2], option[3])
                    if key not in seen:
                        seen.add(key)
                        unique.append(option)
                by_cut[cut] = unique
        return by_cut, sorted(by_cut), None

    def _precompute(self, records: list[GeneRecord]) -> dict[str, dict[str, Any]]:
        started = time.perf_counter()
        info: dict[str, dict[str, Any]] = {}
        for index, record in enumerate(records, start=1):
            by_cut, cuts, error = self._build_cut_index(record.sequence)
            info[record.name] = {
                "seq": record.sequence,
                "length": len(record.sequence),
                "by_cut": by_cut,
                "cuts": cuts,
                "option_count": sum(len(options) for options in by_cut.values()),
                "error": error,
            }
            if self.config.show_progress:
                print(
                    f"  indexed {index}/{len(records)}: {record.name} "
                    f"({len(record.sequence)} nt, {len(cuts)} cut positions)"
                )
        if self.config.show_progress:
            print(f"  cut indexing: {time.perf_counter() - started:.3f} seconds")
        return info

    def _frag_max(self, index: int, count: int) -> int:
        if index == count:
            return self.config.opool_length - 62
        return self.config.opool_length - 58

    def _max_bases(self, count: int) -> int:
        return (
            sum(self._frag_max(index, count) for index in range(1, count + 1))
            - 4 * (count - 1)
        )

    def _remaining_capacity(self, first_index: int, count: int) -> int:
        indices = list(range(first_index, count + 1))
        if not indices:
            return 0
        return (
            sum(self._frag_max(index, count) for index in indices)
            - 4 * (len(indices) - 1)
        )

    @staticmethod
    def _sequence_with_patches(sequence: str, patches: tuple[tuple[int, str], ...]) -> str:
        chars = list(sequence)
        for codon_start, new_codon in patches:
            chars[codon_start:codon_start + 3] = new_codon
        return "".join(chars)

    @staticmethod
    def _window_with_patches(
        sequence: str,
        start: int,
        end: int,
        patches: tuple[tuple[int, str], ...],
    ) -> str:
        chars = list(sequence[start:end])
        for codon_start, new_codon in patches:
            for position in range(max(start, codon_start), min(end, codon_start + 3)):
                chars[position - start] = new_codon[position - codon_start]
        return "".join(chars)

    def _find_path(
        self,
        info: dict[str, Any],
        free_overhangs: set[str],
        fragment_count: int,
    ) -> tuple[dict[str, Any] | None, int]:
        sequence = info["seq"]
        length = info["length"]
        by_cut = info["by_cut"]
        all_cuts = info["cuts"]
        free_mask = sum(self.overhang_bits[overhang] for overhang in free_overhangs)
        failed_states: set[tuple] = set()
        visited_nodes = 0

        def search(
            fragment_index: int,
            fragment_start: int,
            used_mask: int,
            patches: tuple[tuple[int, str], ...],
        ) -> tuple[list[int], list[str], tuple[tuple[int, str], ...]] | None:
            nonlocal visited_nodes
            visited_nodes += 1
            if fragment_index == fragment_count:
                if length - fragment_start > self._frag_max(fragment_index, fragment_count):
                    return None
                modified = self._sequence_with_patches(sequence, patches)
                return None if self._contains_forbidden_site(modified) else ([], [], patches)

            state = (fragment_index, fragment_start, used_mask, patches)
            if state in failed_states:
                return None
            remaining_capacity = self._remaining_capacity(fragment_index + 1, fragment_count)
            minimum_cut = max(
                4 if fragment_index == 1 else fragment_start + 5,
                length + 4 - remaining_capacity,
            )
            maximum_cut = min(
                length,
                fragment_start + self._frag_max(fragment_index, fragment_count),
            )
            left = bisect.bisect_left(all_cuts, minimum_cut)
            right = bisect.bisect_right(all_cuts, maximum_cut)
            fragments_remaining = fragment_count - fragment_index + 1
            target_length = math.ceil(
                (length - fragment_start + 4 * (fragments_remaining - 1))
                / fragments_remaining
            )
            target_cut = fragment_start + target_length
            feasible_cuts = sorted(
                all_cuts[left:right],
                key=lambda cut: (abs(cut - target_cut), -cut),
            )
            current_patches = dict(patches)
            for cut in feasible_cuts:
                options = sorted(
                    by_cut[cut],
                    key=lambda option: (
                        option[1] != "Native",
                        self.overhang_rank[option[0]],
                    ),
                )
                for overhang, _kind, codon_start, new_codon in options:
                    bit = self.overhang_bits[overhang]
                    if not (free_mask & bit) or (used_mask & bit):
                        continue
                    next_patches = dict(current_patches)
                    if codon_start is not None:
                        if codon_start in next_patches and next_patches[codon_start] != new_codon:
                            continue
                        next_patches[codon_start] = new_codon
                    next_tuple = tuple(sorted(next_patches.items()))
                    if self._window_with_patches(sequence, cut - 4, cut, next_tuple) != overhang:
                        continue
                    result = search(
                        fragment_index + 1,
                        cut - 4,
                        used_mask | bit,
                        next_tuple,
                    )
                    if result is not None:
                        later_cuts, later_overhangs, final_patches = result
                        return [cut, *later_cuts], [overhang, *later_overhangs], final_patches
            failed_states.add(state)
            return None

        result = search(1, 0, 0, tuple())
        if result is None:
            return None, visited_nodes
        cuts, overhangs, patches = result
        return {
            "cuts": cuts,
            "overhangs": overhangs,
            "sequence": self._sequence_with_patches(sequence, patches),
        }, visited_nodes

    def _long_pool_capacity(self, fragment_count: int) -> int:
        capacities = [len(self.internal_overhangs) // (fragment_count - 1)]
        if self.config.genes_per_subpool is not None:
            capacities.append(self.config.genes_per_subpool)
        return min(capacities)

    @staticmethod
    def _fragments_from_cuts(sequence: str, cuts: list[int]) -> list[str]:
        fragments: list[str] = []
        start = 0
        for cut in cuts:
            fragments.append(sequence[start:cut])
            start = cut - 4
        fragments.append(sequence[start:])
        return fragments

    def _build_pools_for_count(
        self,
        genes: list[GeneRecord],
        gene_info: dict[str, dict[str, Any]],
        fragment_count: int,
        first_block: int,
    ) -> tuple[list[dict[str, Any]], list[GeneRecord], int, int]:
        rows: list[dict[str, Any]] = []
        remaining = sorted(
            genes,
            key=lambda record: (
                gene_info[record.name]["option_count"],
                -gene_info[record.name]["length"],
                record.name,
            ),
        )
        block = first_block
        capacity = self._long_pool_capacity(fragment_count)
        total_nodes = 0
        if capacity <= 0:
            return rows, remaining, block, total_nodes

        while remaining:
            selected: list[tuple] = []
            free_internal = set(self.internal_overhangs)
            index = 0
            while index < len(remaining) and len(selected) < capacity:
                record = remaining[index]
                info = gene_info[record.name]
                if info["error"] or not info["cuts"]:
                    index += 1
                    continue
                if info["length"] > self._max_bases(fragment_count):
                    index += 1
                    continue
                result, nodes = self._find_path(info, free_internal, fragment_count)
                total_nodes += nodes
                if result is None:
                    index += 1
                    continue
                selected.append((record, result))
                free_internal.difference_update(result["overhangs"])
                remaining.pop(index)

            if not selected:
                break
            for record, result in selected:
                fragments = self._fragments_from_cuts(result["sequence"], result["cuts"])
                row: dict[str, Any] = {
                    "Block": block,
                    "Length Distribution": f"Long-{fragment_count}part",
                    "Sequence Name": record.name,
                    "VectorOH1": self.config.vector_oh1,
                    "VectorOH1_Source": "Fixed",
                    "VectorOH2": self.config.vector_oh2,
                    "Full Sequence": result["sequence"],
                }
                for number, overhang in enumerate(result["overhangs"], start=1):
                    row[f"Overhang{number}"] = overhang
                for number, fragment in enumerate(fragments, start=1):
                    row[f"DNA Fragment {number}"] = fragment
                rows.append(row)
            if self.config.show_progress:
                print(
                    f"  block {block}: {len(selected)} gene(s), "
                    f"{fragment_count} fragments each"
                )
            block += 1
        return rows, remaining, block, total_nodes

    def _short_capacity(self, number_of_genes: int) -> int:
        capacities = [max(1, number_of_genes)]
        if self.config.short_pool_max_size is not None:
            capacities.append(self.config.short_pool_max_size)
        if self.config.genes_per_subpool is not None:
            capacities.append(self.config.genes_per_subpool)
        return min(capacities)

    def _assign_short(self, records: list[GeneRecord], first_block: int = 1) -> list[dict[str, Any]]:
        if not records:
            return []
        rows: list[dict[str, Any]] = []
        capacity = self._short_capacity(len(records))
        for start in range(0, len(records), capacity):
            block = first_block + start // capacity
            for record in records[start:start + capacity]:
                rows.append({
                    "Block": block,
                    "Length Distribution": "Short",
                    "Sequence Name": record.name,
                    "VectorOH1": self.config.vector_oh1,
                    "VectorOH1_Source": "Fixed",
                    "VectorOH2": self.config.vector_oh2,
                    "Overhang1": "N/A",
                    "Overhang2": "N/A",
                    "DNA Fragment 1": record.sequence,
                    "DNA Fragment 2": "",
                    "DNA Fragment 3": "",
                    "Full Sequence": record.sequence,
                })
        return rows

    @staticmethod
    def _nonempty(value: Any) -> bool:
        if value is None or (isinstance(value, float) and np.isnan(value)):
            return False
        return str(value).strip().upper() not in {"", "N/A", "NAN"}

    def _validate(self, df: pd.DataFrame, source_records: list[GeneRecord]) -> None:
        if df.empty:
            return
        source = {record.name: record.sequence for record in source_records}
        fragment_columns = sorted(
            [column for column in df if column.startswith("DNA Fragment ")],
            key=lambda column: int(column.replace("DNA Fragment ", "")),
        )
        overhang_columns = sorted(
            [column for column in df if column.startswith("Overhang")],
            key=lambda column: int(column.replace("Overhang", "")),
        )
        for _, row in df.iterrows():
            name = str(row["Sequence Name"])
            full = str(row["Full Sequence"])
            fragments = [str(row[column]).strip().upper() for column in fragment_columns if self._nonempty(row.get(column))]
            overhangs = [str(row[column]).strip().upper() for column in overhang_columns if self._nonempty(row.get(column))]
            if len(overhangs) != max(0, len(fragments) - 1):
                raise AssertionError(f"{name}: fragment/overhang count mismatch")
            reconstructed = fragments[0] + "".join(fragment[4:] for fragment in fragments[1:])
            if reconstructed != full:
                raise AssertionError(f"{name}: fragments do not reconstruct Full Sequence")
            for index, overhang in enumerate(overhangs):
                if fragments[index][-4:] != overhang or fragments[index + 1][:4] != overhang:
                    raise AssertionError(f"{name}: invalid overlap {index + 1}")
            for index, fragment in enumerate(fragments, start=1):
                maximum = self._frag_max(index, len(fragments))
                if len(fragment) > maximum:
                    raise AssertionError(f"{name}: fragment {index} exceeds {maximum} nt")
            if self._contains_forbidden_site(full):
                raise AssertionError(f"{name}: forbidden Type IIS site in Full Sequence")
            if str(Seq(full).translate()) != str(Seq(source[name]).translate()):
                raise AssertionError(f"{name}: a synonymous edit changed translation")
        for block, block_df in df.groupby("Block"):
            internal = [
                str(value).strip().upper()
                for column in overhang_columns
                for value in block_df[column]
                if self._nonempty(value)
            ]
            if len(internal) != len(set(internal)):
                raise AssertionError(f"Block {block}: an internal overhang is reused")
            if self.config.genes_per_subpool is not None and len(block_df) > self.config.genes_per_subpool:
                raise AssertionError(f"Block {block}: exceeds genes_per_subpool")

    def _write_used_overhangs(self, df: pd.DataFrame) -> None:
        columns = [column for column in df if column.startswith("Overhang")]
        used = {
            str(value).strip().upper()
            for column in columns
            for value in df[column]
            if self._nonempty(value)
        }
        ordered = [overhang for overhang in self.overhangs if overhang in used]
        pd.DataFrame([[",".join(ordered)]]).to_csv(
            self.paths.overhangs_used, index=False, header=False
        )

    def run(self) -> tuple[pd.DataFrame, list[str], int]:
        records, strip_log = self._load_records()
        pd.DataFrame(strip_log).to_csv(self.paths.strip_log, index=False)
        short_limit = self.config.opool_length - 62
        short = [record for record in records if len(record.sequence) <= short_limit]
        pending = [record for record in records if len(record.sequence) > short_limit]
        rows = self._assign_short(short)
        block = max((row["Block"] for row in rows), default=0) + 1
        gene_info = self._precompute(pending)
        total_nodes = 0
        for fragment_count in range(2, self.config.max_fragments + 1):
            if not pending or fragment_count - 1 > len(self.internal_overhangs):
                break
            eligible = [
                record for record in pending
                if len(record.sequence) <= self._max_bases(fragment_count)
            ]
            if not eligible:
                continue
            new_rows, _failed, block, nodes = self._build_pools_for_count(
                eligible, gene_info, fragment_count, block
            )
            total_nodes += nodes
            rows.extend(new_rows)
            assigned = {str(row["Sequence Name"]) for row in new_rows}
            pending = [record for record in pending if record.name not in assigned]
        unassigned = [record.name for record in pending]
        df = pd.DataFrame(rows)
        if not df.empty:
            df.sort_values(["Block", "Sequence Name"], inplace=True)
            overhang_columns = sorted(
                [column for column in df if column.startswith("Overhang")],
                key=lambda column: int(column.replace("Overhang", "")),
            )
            fragment_columns = sorted(
                [column for column in df if column.startswith("DNA Fragment ")],
                key=lambda column: int(column.replace("DNA Fragment ", "")),
            )
            order = [
                "Block", "Length Distribution", "Sequence Name",
                "VectorOH1", "VectorOH1_Source", "VectorOH2",
                *overhang_columns, *fragment_columns, "Full Sequence",
            ]
            df = df[[column for column in order if column in df]]
        self._validate(df, records)
        df.to_csv(self.paths.assigned, index=False)
        pd.DataFrame({"Sequence Name": unassigned}).to_csv(self.paths.unassigned, index=False)
        self._write_used_overhangs(df)
        if self.config.show_progress:
            print(
                f"[2/3] Pool assignment: {len(df)} assigned, {len(unassigned)} unassigned, "
                f"{df['Block'].nunique() if not df.empty else 0} blocks, {total_nodes:,} search states"
            )
        return df, unassigned, total_nodes


class PrimerAssembler:
    """Assign orthogonal primer pairs and build order-ready oligos."""

    def __init__(self, config: WorkflowConfig, paths: RunPaths, primers_path: Path) -> None:
        self.config = config
        self.paths = paths
        self.primers_path = primers_path
        self.rng = np.random.default_rng(config.random_seed)

    @staticmethod
    def _nonempty(value: Any) -> bool:
        if value is None or (isinstance(value, float) and np.isnan(value)):
            return False
        return str(value).strip().upper() not in {"", "N/A", "NAN"}

    def _load_primers(self) -> list[tuple[str, str]]:
        df = pd.read_csv(self.primers_path, header=None)
        df = df.replace(r"^\s*$", np.nan, regex=True).dropna(how="all").reset_index(drop=True)
        if df.shape[1] < 2:
            raise ValueError("Primer CSV requires at least two columns: name and sequence.")
        primers = [
            (str(name).strip(), _normalize_dna(sequence, f"primer {name}"))
            for name, sequence in zip(df.iloc[:, 0], df.iloc[:, -1])
            if str(name).strip() and str(sequence).strip()
        ]
        if not primers:
            raise ValueError(f"No primers found in {self.primers_path}")
        if len({name for name, _ in primers}) != len(primers):
            raise ValueError("Primer names must be unique.")
        return primers

    def _primer_is_safe(self, sequence: str) -> bool:
        return not any(site in sequence for site in self.config.forbidden_typeiis_sites)

    @staticmethod
    def _names_equal(left: Any, right: Any) -> bool:
        try:
            return int(str(left).strip()) == int(str(right).strip())
        except ValueError:
            return str(left).strip() == str(right).strip()

    def _all_combinatorial_pairs(self, primers: list[tuple[str, str]]) -> pd.DataFrame:
        if self.config.exclude_unsafe_primers:
            excluded = [name for name, sequence in primers if not self._primer_is_safe(sequence)]
            primers = [(name, sequence) for name, sequence in primers if self._primer_is_safe(sequence)]
            if excluded and self.config.show_progress:
                print(
                    f"  excluded {len(excluded)} primer(s) containing a Type IIS site: "
                    f"{', '.join(excluded[:8])}"
                )
        rows = []
        pair_id = 1
        for forward_index in range(len(primers)):
            forward_name, forward_sequence = primers[forward_index]
            for reverse_index in range(forward_index + 1, len(primers)):
                reverse_name, reverse_sequence = primers[reverse_index]
                rows.append({
                    "Pair_ID": pair_id,
                    "Fwd_Name": forward_name,
                    "Fwd_Seq": forward_sequence,
                    "Rev_Name": reverse_name,
                    "Rev_Seq": reverse_sequence,
                })
                pair_id += 1
        result = pd.DataFrame(rows)
        result.to_csv(self.paths.all_pairs, index=False)
        return result

    def _pairs_for_blocks(
        self,
        primers: list[tuple[str, str]],
        blocks: list[int],
    ) -> tuple[dict[int, tuple[str, str, str, str]], pd.DataFrame | None]:
        if self.config.primer_mode == "combinatorial":
            pairs = self._all_combinatorial_pairs(primers)
            if len(pairs) < len(blocks):
                raise ValueError(
                    f"Need {len(blocks)} safe primer pairs but only {len(pairs)} are available."
                )
            mapping = {
                block: (
                    str(pairs.iloc[index]["Fwd_Name"]),
                    str(pairs.iloc[index]["Fwd_Seq"]),
                    str(pairs.iloc[index]["Rev_Name"]),
                    str(pairs.iloc[index]["Rev_Seq"]),
                )
                for index, block in enumerate(blocks)
            }
            unused = pairs.iloc[len(blocks):].reset_index(drop=True)
            unused.to_csv(self.paths.unused_pairs, index=False)
            return mapping, unused

        if self.config.primer_start_at is None or not str(self.config.primer_start_at).strip():
            start_index = 0
        else:
            start_index = next(
                (
                    index for index, (name, _sequence) in enumerate(primers)
                    if self._names_equal(name, self.config.primer_start_at)
                ),
                -1,
            )
            if start_index < 0:
                raise ValueError(f"Primer {self.config.primer_start_at!r} was not found.")
        mapping = {}
        for offset, block in enumerate(blocks):
            forward_index = start_index + 2 * offset
            reverse_index = forward_index + 1
            if reverse_index >= len(primers):
                raise ValueError(f"Not enough sequential primer pairs for block {block}.")
            forward_name, forward_sequence = primers[forward_index]
            reverse_name, reverse_sequence = primers[reverse_index]
            if self.config.exclude_unsafe_primers and (
                not self._primer_is_safe(forward_sequence)
                or not self._primer_is_safe(reverse_sequence)
            ):
                raise ValueError(
                    f"Sequential primer pair {forward_name}/{reverse_name} contains a Type IIS site. "
                    "Choose a different primer_start_at or disable exclude_unsafe_primers."
                )
            mapping[block] = (
                forward_name,
                forward_sequence,
                reverse_name,
                reverse_sequence,
            )
        return mapping, None

    def _vec1_prefix(self, row: pd.Series) -> str:
        vec1 = (
            str(row.get("VectorOH1", "")).strip().upper()
            if self._nonempty(row.get("VectorOH1"))
            else self.config.vector_oh1
        )
        return vec1

    def _gc_fraction(self, sequence: str) -> float:
        return sum(base in "GC" for base in sequence) / len(sequence) if sequence else 0.0

    @staticmethod
    def _has_long_homopolymer(sequence: str, maximum: int) -> bool:
        run = 1
        for index in range(1, len(sequence)):
            run = run + 1 if sequence[index] == sequence[index - 1] else 1
            if run > maximum:
                return True
        return False

    @staticmethod
    def _junction_creates_site(left: str, right: str, sites: Iterable[str]) -> bool:
        sites = list(sites)
        flank = max(len(site) for site in sites) - 1
        left_tail = left[-flank:] if flank else ""
        right_head = right[:flank] if flank else ""
        combined = left_tail + right_head
        boundary = len(left_tail)
        return any(
            start < boundary < start + len(site)
            for site in sites
            for start in range(len(combined) - len(site) + 1)
            if combined.startswith(site, start)
        )

    @staticmethod
    def _same_order_subsequences(sequence: str, length: int) -> set[str]:
        return {
            "".join(sequence[index] for index in indices)
            for indices in combinations(range(len(sequence)), length)
        }

    def _stuffer_prefixes(self) -> set[str]:
        return set().union(*(
            self._same_order_subsequences(site, 5)
            for site in self.config.forbidden_typeiis_sites
        ))

    def _stuffer_valid(
        self,
        sequence: str,
        left: str,
        right: str,
        enforce_gc: bool,
    ) -> bool:
        if set(sequence) - DNA_BASES:
            return False
        if enforce_gc and not self.config.stuffer_gc_min <= self._gc_fraction(sequence) <= self.config.stuffer_gc_max:
            return False
        if self._has_long_homopolymer(sequence, self.config.stuffer_max_homopolymer):
            return False
        if any(site in sequence for site in self.config.forbidden_typeiis_sites):
            return False
        if len(sequence) >= 5 and sequence[:5] in self._stuffer_prefixes():
            return False
        if self._junction_creates_site(left, sequence, self.config.forbidden_typeiis_sites):
            return False
        if self._junction_creates_site(sequence, right, self.config.forbidden_typeiis_sites):
            return False
        return True

    def _make_stuffer(self, length: int, left: str, right: str) -> str:
        if length < 0:
            raise ValueError("The fixed oligo elements exceed opool_length.")
        if length == 0:
            if self._stuffer_valid("", left, right, enforce_gc=False):
                return ""
            raise ValueError("A zero-length stuffer creates a forbidden Type IIS junction.")
        if length < 5:
            candidates = ["".join(chars) for chars in product("ACGT", repeat=length)]
            self.rng.shuffle(candidates)
            distance = lambda candidate: max(
                self.config.stuffer_gc_min - self._gc_fraction(candidate),
                0.0,
                self._gc_fraction(candidate) - self.config.stuffer_gc_max,
            )
            best = min(map(distance, candidates))
            for candidate in candidates:
                if distance(candidate) == best and self._stuffer_valid(
                    candidate, left, right, enforce_gc=False
                ):
                    return candidate
            raise ValueError(f"No valid {length}-nt stuffer exists.")
        for _ in range(self.config.stuffer_max_tries):
            candidate = "".join(self.rng.choice(list("ACGT"), size=length))
            if self._stuffer_valid(candidate, left, right, enforce_gc=True):
                return candidate
        raise RuntimeError(f"Could not generate a valid {length}-nt stuffer.")

    @staticmethod
    def _count_overlapping(sequence: str, motif: str) -> int:
        return sum(sequence.startswith(motif, index) for index in range(len(sequence) - len(motif) + 1))

    def _build_oligo(
        self,
        insert: str,
        forward_primer: str,
        reverse_primer: str,
    ) -> str:
        cut_prefix = self.config.typeiis_site + self.config.typeiis_n
        cut_suffix = reverse_complement(cut_prefix)
        reverse_primer_suffix = reverse_complement(reverse_primer)
        core = cut_prefix + insert + cut_suffix
        fixed_length = len(forward_primer) + len(core) + len(reverse_primer_suffix)
        if self.config.add_stuffer:
            stuffer = self._make_stuffer(
                self.config.opool_length - fixed_length,
                forward_primer,
                cut_prefix,
            )
        else:
            stuffer = ""
        oligo = forward_primer + stuffer + core + reverse_primer_suffix
        if self.config.add_stuffer and len(oligo) != self.config.opool_length:
            raise AssertionError("Oligo length does not match opool_length.")
        expected = {
            site: self._count_overlapping(cut_prefix + cut_suffix, site)
            for site in self.config.forbidden_typeiis_sites
        }
        for site, expected_count in expected.items():
            actual = self._count_overlapping(oligo, site)
            if actual != expected_count:
                raise ValueError(
                    f"Final oligo contains {actual} copies of {site}; expected {expected_count}. "
                    "Check primers, vector-overhang junctions, and inserts."
                )
        return oligo

    def run(self, assigned_df: pd.DataFrame) -> pd.DataFrame:
        if assigned_df.empty:
            raise ValueError("No assigned genes are available for primer assignment.")
        primers = self._load_primers()
        blocks = sorted(pd.to_numeric(assigned_df["Block"]).astype(int).unique())
        pairs, _unused_pairs = self._pairs_for_blocks(primers, blocks)
        fragment_columns = sorted(
            [column for column in assigned_df if column.startswith("DNA Fragment ")],
            key=lambda column: int(column.replace("DNA Fragment ", "")),
        )
        output_blocks = []
        for block in blocks:
            forward_name, forward_sequence, reverse_name, reverse_sequence = pairs[block]
            block_df = assigned_df[assigned_df["Block"] == block].copy()
            for column in fragment_columns:
                block_df[column] = block_df[column].astype(object)
            for index, row in block_df.iterrows():
                nonempty = [column for column in fragment_columns if self._nonempty(row.get(column))]
                for fragment_index, column in enumerate(fragment_columns):
                    if column not in nonempty:
                        block_df.at[index, column] = ""
                        continue
                    insert = str(row[column]).strip().upper()
                    if column == nonempty[0]:
                        insert = self._vec1_prefix(row) + insert
                    if column == nonempty[-1]:
                        vec2 = (
                            str(row.get("VectorOH2", "")).strip().upper()
                            if self._nonempty(row.get("VectorOH2"))
                            else self.config.vector_oh2
                        )
                        insert += vec2
                    block_df.at[index, column] = self._build_oligo(
                        insert, forward_sequence, reverse_sequence
                    )
                block_df.at[index, "Primer_Forward_Name"] = forward_name
                block_df.at[index, "Primer_Forward_Seq"] = forward_sequence
                block_df.at[index, "Primer_Reverse_Name"] = reverse_name
                block_df.at[index, "Primer_Reverse_Seq"] = reverse_sequence
            output_blocks.append(block_df)
        output = pd.concat(output_blocks, ignore_index=True).sort_values("Block")
        output.to_csv(self.paths.full_info, index=False)

        used_names = set(output["Primer_Forward_Name"].astype(str)) | set(output["Primer_Reverse_Name"].astype(str))
        unused_primers = [(name, sequence) for name, sequence in primers if name not in used_names]
        pd.DataFrame(unused_primers).to_csv(self.paths.unused_primers, index=False, header=False)

        with self.paths.references_fasta.open("w") as handle:
            for _, row in output.iterrows():
                handle.write(f">Block_{int(row['Block'])}_{row['Sequence Name']}\n{row['Full Sequence']}\n")
        fragment_rows = [
            {
                "Gene_Fragment": f"{row['Sequence Name']}_Fragment{number}",
                "Sequence": str(row[column]),
            }
            for _, row in output.iterrows()
            for number, column in enumerate(fragment_columns, start=1)
            if self._nonempty(row.get(column))
        ]
        pd.DataFrame(fragment_rows).to_csv(self.paths.fragments, index=False)
        if self.config.show_progress:
            print(
                f"[3/3] Primer assignment: {len(blocks)} pair(s), "
                f"{len(fragment_rows)} order fragments"
            )
        return output


def _planned_outputs(config: WorkflowConfig, paths: RunPaths, input_kind: str) -> list[Path]:
    outputs = [
        paths.assigned,
        paths.strip_log,
        paths.unassigned,
        paths.overhangs_used,
    ]
    if input_kind == "aa" or config.input_path.resolve() != paths.optimized.resolve():
        outputs.insert(0, paths.optimized)
    if not config.skip_primer_assignment:
        outputs.extend([
            paths.full_info,
            paths.references_fasta,
            paths.fragments,
            paths.unused_primers,
        ])
        if config.primer_mode == "combinatorial":
            outputs.extend([paths.all_pairs, paths.unused_pairs])
    return outputs


def _check_overwrite(config: WorkflowConfig, paths: RunPaths, input_kind: str) -> None:
    existing = [path for path in _planned_outputs(config, paths, input_kind) if path.exists()]
    if existing and not config.force:
        formatted = "\n".join(f"  - {path}" for path in existing)
        raise FileExistsError(
            "Refusing to overwrite existing output files:\n"
            f"{formatted}\n"
            "Choose a different --run-name/--output-dir or pass --force."
        )


def describe_config(
    config: WorkflowConfig,
    paths: RunPaths,
    overhangs_path: Path,
    primers_path: Path,
) -> str:
    genes_per_pool = "automatic" if config.genes_per_subpool is None else str(config.genes_per_subpool)
    lines = [
        "oPool workflow configuration",
        f"  input:              {config.input_path.resolve()}",
        f"  input kind:         {config.input_kind}",
        f"  output directory:   {paths.run_dir}",
        f"  run name:           {paths.run_name}",
        f"  overhang inventory: {overhangs_path}",
        f"  primer inventory:   {primers_path}",
        f"  oligo length:       {config.opool_length}",
        f"  genes per sub-pool: {genes_per_pool}",
        f"  vector overhangs:   {config.vector_oh1} / {config.vector_oh2}",
        f"  Type IIS element:   {config.typeiis_site} + {config.typeiis_n}",
        f"  primer assignment:  {'skipped' if config.skip_primer_assignment else config.primer_mode}",
    ]
    return "\n".join(lines)


def run_workflow(config: WorkflowConfig) -> WorkflowResult:
    """Run the complete workflow and return paths plus a compact summary."""
    started = time.perf_counter()
    paths = RunPaths.from_config(config)
    paths.run_dir.mkdir(parents=True, exist_ok=True)
    overhangs_path = _resolve_inventory(config.overhangs_path, "overhangs.csv")
    primers_path = _resolve_inventory(config.primers_path, "orthogonal_oligos.csv")
    input_kind = _detect_input_kind(config.input_path.resolve(), config.input_kind)
    _check_overwrite(config, paths, input_kind)
    if config.show_progress:
        print(describe_config(config, paths, overhangs_path, primers_path))

    optimized_df = prepare_optimized_input(config, paths)
    pool_assigner = PoolAssigner(config, paths, optimized_df, overhangs_path)
    assigned_df, unassigned, _search_nodes = pool_assigner.run()
    if not config.skip_primer_assignment and not assigned_df.empty:
        PrimerAssembler(config, paths, primers_path).run(assigned_df)

    output_files = [
        path for path in _planned_outputs(config, paths, input_kind)
        if path.exists()
    ]
    result = WorkflowResult(
        paths=paths,
        optimized_genes=len(optimized_df),
        assigned_genes=len(assigned_df),
        unassigned_genes=len(unassigned),
        blocks=int(assigned_df["Block"].nunique()) if not assigned_df.empty else 0,
        runtime_seconds=time.perf_counter() - started,
        output_files=output_files,
    )
    if config.show_progress:
        print(
            f"[DONE] {result.assigned_genes}/{result.optimized_genes} genes assigned "
            f"in {result.blocks} blocks ({result.runtime_seconds:.3f} seconds)"
        )
        print("Outputs:")
        for path in result.output_files:
            print(f"  {path}")
    return result
