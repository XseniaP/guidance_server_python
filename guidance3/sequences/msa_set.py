"""A user-supplied set of alternative MSAs of the same sequences: read, validate and rank with the DL model.

Accepted input: a folder of FASTA alignments, or a .zip / .tar / .tar.gz / .tgz archive of them.
Every MSA must contain the same sequence names (each once), be aligned (all rows the same length) and
contain at least 4 sequences. Files are recognised by extension: .fasta .fas .fa .faa .fna .aln .afa .mfa

The DL model (Serok et al.) predicts each MSA's distance from the true alignment, so LOWER predicted_score = better.
"""
import json
import os
import shutil
import subprocess
import tarfile
import tempfile
import zipfile

MSA_EXTENSIONS = (".fasta", ".fas", ".fa", ".faa", ".fna", ".aln", ".afa", ".mfa")
MAX_MSAS = 2000
MAX_ARCHIVE_BYTES = 2 * 1000 ** 3  # uncompressed
TRUE_MSA_FILENAME = "MSA_TRUE.fas"  # placeholder reference required by features_for_msas (not used for ranking)


class MSASetError(ValueError):
    """Problem with the user's MSA set (shown to the user)."""


# ------------------------------------------------------------------------------------------------ reading

def _is_msa_file(name):
    base = os.path.basename(name)
    return not base.startswith(".") and not base.startswith("__MACOSX") and base.lower().endswith(MSA_EXTENSIONS)


def _safe_member(name):
    norm = os.path.normpath(name)
    return not (os.path.isabs(norm) or norm.startswith("..") or "__MACOSX" in norm.split(os.sep))


def collect_msa_files(src, dest_dir):
    """Copy/extract the MSAs of `src` (folder or archive) into dest_dir (flat). Returns sorted file paths."""
    os.makedirs(dest_dir, exist_ok=True)
    out = []

    def add(name, data_or_path, from_file):
        base = os.path.basename(name)
        if base in {os.path.basename(p) for p in out}:
            raise MSASetError(f"Two MSA files have the same name: {base}")
        dst = os.path.join(dest_dir, base)
        if from_file:
            shutil.copyfile(data_or_path, dst)
        else:
            with open(dst, "wb") as fh:
                fh.write(data_or_path)
        out.append(dst)

    if os.path.isdir(src):
        for root, _dirs, files in os.walk(src):
            for f in sorted(files):
                if _is_msa_file(f):
                    add(f, os.path.join(root, f), True)
    elif zipfile.is_zipfile(src):
        with zipfile.ZipFile(src) as z:
            infos = [i for i in z.infolist() if not i.is_dir() and _safe_member(i.filename) and _is_msa_file(i.filename)]
            if sum(i.file_size for i in infos) > MAX_ARCHIVE_BYTES:
                raise MSASetError("The MSA archive is too large.")
            for i in infos:
                add(i.filename, z.read(i), False)
    elif tarfile.is_tarfile(src):
        with tarfile.open(src) as t:
            members = [m for m in t.getmembers() if m.isfile() and _safe_member(m.name) and _is_msa_file(m.name)]
            if sum(m.size for m in members) > MAX_ARCHIVE_BYTES:
                raise MSASetError("The MSA archive is too large.")
            for m in members:
                add(m.name, t.extractfile(m).read(), False)
    else:
        raise MSASetError("The MSA set must be a folder or a .zip / .tar.gz archive of FASTA alignments.")

    if len(out) > MAX_MSAS:
        raise MSASetError(f"Too many MSAs ({len(out)}); at most {MAX_MSAS} are allowed.")
    return sorted(out)


def read_fasta(path):
    """[(name, sequence)] of a FASTA file (name = full header line after '>')."""
    records, name, chunks = [], None, []
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name, chunks = line[1:].strip(), []
            elif name is None:
                raise MSASetError(f"{os.path.basename(path)} is not in FASTA format.")
            else:
                chunks.append(line.replace(" ", ""))
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def validate_msa_set(files):
    """Raise MSASetError unless the files are >= 2 aligned FASTA MSAs of the same >= 4 sequences."""
    if len(files) < 2:
        raise MSASetError(f"The MSA set must contain at least 2 alignments (found {len(files)}). "
                          f"Accepted extensions: {' '.join(MSA_EXTENSIONS)}")
    ref_names, ref_file, ref_seqs = None, None, None
    for f in files:
        base = os.path.basename(f)
        recs = read_fasta(f)
        if len(recs) < 4:
            raise MSASetError(f"{base}: at least 4 sequences are required (found {len(recs)}).")
        names = [n for n, _s in recs]
        if len(set(names)) != len(names):
            raise MSASetError(f"{base}: sequence names must be unique.")
        lengths = {len(s) for _n, s in recs}
        if len(lengths) != 1:
            raise MSASetError(f"{base}: is not aligned (sequences have different lengths).")
        ungapped = {n: s.replace("-", "").replace(".", "").upper() for n, s in recs}
        if ref_names is None:
            ref_names, ref_file, ref_seqs = set(names), base, ungapped
            continue
        if set(names) != ref_names:
            diff = sorted(set(names) ^ ref_names)[:5]
            raise MSASetError(f"{base} and {ref_file} do not contain the same sequence names "
                              f"(e.g. {', '.join(diff)}).")
        bad = [n for n in names if ungapped[n] != ref_seqs[n]]
        if bad:
            raise MSASetError(f"{base} and {ref_file} align different sequences under the same name "
                              f"(e.g. {bad[0]}). All MSAs must be alignments of the same sequences.")


# ------------------------------------------------------------------------------------------------ ranking

def _copy_stripped(src, dst):
    """features_for_msas crashes on blank lines."""
    with open(src) as fin, open(dst, "w") as fout:
        for line in fin:
            if line.strip():
                fout.write(line)


def rank_msas(files, is_nucleotide, out_dir, verbose=False, timeout=3600):
    """Rank MSAs with the pretrained DL model. Returns a pandas DataFrame (code, predicted_score), best first
    (lowest predicted distance from the true MSA). Writes <out_dir>/features.csv. Raises RuntimeError on failure."""
    import pandas as pd
    from guidance3.constants import (
        FEATURES_EXTRACTION_PROG, FEATURES_EXTRACTION_MATRIX_DIR, DL_MODEL_PYTHON, DL_MODEL_PREDICT_SCRIPT,
        DL_MODEL_PATH, DL_MODEL_SCALER_PATH, DL_MODEL_NUC_PATH, DL_MODEL_NUC_SCALER_PATH,
    )
    from guidance3.dl_model.scripts.predict_best_msa_main import AA_MODELS, NUC_MODELS

    os.makedirs(out_dir, exist_ok=True)
    work = tempfile.mkdtemp(prefix="guidance3_rank_", dir=out_dir)
    try:
        all_msas = os.path.join(work, "features_input", "all_msas")
        os.makedirs(all_msas)
        for f in files:
            _copy_stripped(f, os.path.join(all_msas, os.path.basename(f)))
        _copy_stripped(files[0], os.path.join(all_msas, TRUE_MSA_FILENAME))  # structural placeholder only

        cfg = {
            "models_list": NUC_MODELS if is_nucleotide else AA_MODELS,
            "sop_calc_type": 1,
            "additional_weights": ["HENIKOFF_WG", "HENIKOFF_WOG", "CLUSTAL_MID_ROOT", "CLUSTAL_DIFFERENTIAL_SUM"],
            "k_values": [8, 16, 32] if is_nucleotide else [5, 10, 20],
            "stats_output": ["ALL"] if is_nucleotide else ["ALL_NO_SUBS_MATRIX"],
            "input_files_dir_path": os.path.join(work, "features_input"),
            "output_file_dir_path": work,
            "is_unified_file": True,
            "matrix_dir_path": FEATURES_EXTRACTION_MATRIX_DIR,
        }
        cfg_path = os.path.join(work, "config_features.json")
        with open(cfg_path, "w") as fh:
            json.dump(cfg, fh, indent=4)
        res = subprocess.run([FEATURES_EXTRACTION_PROG, cfg_path], capture_output=True, text=True, timeout=timeout)
        if verbose and res.stdout:
            print(res.stdout)
        features_file = os.path.join(work, "unified_stats_all_msas.csv")
        if res.returncode != 0 or not os.path.exists(features_file):
            raise RuntimeError(f"Feature extraction failed (rc={res.returncode}): {res.stderr[-2000:]}")

        df = pd.read_csv(features_file)
        df.insert(1, "code1", "msa_set")
        df.dropna(axis=1, how="all", inplace=True)
        df = df[df[df.columns[0]] != "all_msas"]
        df = df.drop(columns=[c for c in ("ssp_from_true", "dseq_from_true", "dpos_from_true") if c in df.columns])
        df.to_csv(features_file, index=False)
        df[df["code"] != TRUE_MSA_FILENAME].to_csv(os.path.join(out_dir, "features.csv"), index=False)

        pred_dir = os.path.join(work, "pred")
        cmd = DL_MODEL_PYTHON + [
            DL_MODEL_PREDICT_SCRIPT, "--features-file", features_file,
            "--model-path", DL_MODEL_NUC_PATH if is_nucleotide else DL_MODEL_PATH,
            "--scaler-path", DL_MODEL_NUC_SCALER_PATH if is_nucleotide else DL_MODEL_SCALER_PATH,
            "--scaler-type-features", "rank", "--scaler-type-labels", "rank",
            "--out-dir", pred_dir, "--no-metrics",
        ]
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        env = {**os.environ, "PYTHONPATH": os.path.dirname(root) + os.pathsep + os.environ.get("PYTHONPATH", "")}
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, env=env)
        if verbose and res.stdout:
            print(res.stdout)
        pred_file = os.path.join(pred_dir, "prediction_pretrained_0_mode1_dseq_from_true.csv")
        if res.returncode != 0 or not os.path.exists(pred_file):
            raise RuntimeError(f"DL model failed (rc={res.returncode}): {res.stderr[-2000:]}")
        preds = pd.read_csv(pred_file)
        if not {"code", "predicted_score"} <= set(preds.columns):
            raise RuntimeError("DL predictions are missing the code / predicted_score columns")
        names = {os.path.basename(f) for f in files}
        preds = preds[preds["code"].isin(names)][["code", "predicted_score"]]
        if preds.empty:
            raise RuntimeError("The DL model returned no predictions for the MSA set")
        return preds.sort_values("predicted_score", ascending=True).reset_index(drop=True)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def select_best_from_msa_set(src, work_dir, seq_type, verbose=False):
    """Collect, validate and rank a user's MSA set. seq_type: 'aa' | 'nuc' | 'codon' (codon MSAs are ranked on
    their translation). Writes <work_dir>/MSA_set/ (the MSAs) and <work_dir>/MSA_set_ranking.csv.
    Returns the path of the best MSA (original, untranslated file)."""
    set_dir = os.path.join(work_dir, "MSA_set")
    files = collect_msa_files(src, set_dir)
    validate_msa_set(files)
    rank_files = files
    if seq_type == "codon":
        from guidance3.sequences.fasta import _translate_codon_msa_inplace
        tr_dir = os.path.join(work_dir, "MSA_set_translated")
        os.makedirs(tr_dir, exist_ok=True)
        rank_files = []
        for f in files:
            dst = os.path.join(tr_dir, os.path.basename(f))
            shutil.copyfile(f, dst)
            _translate_codon_msa_inplace(dst)
            rank_files.append(dst)
    print(f"[msa_set] ranking {len(files)} MSAs with the DL model...")
    ranking = rank_msas(rank_files, seq_type == "nuc", os.path.join(work_dir, "MSA_set_ranking"), verbose=verbose)
    ranking.insert(0, "rank", range(1, len(ranking) + 1))
    ranking.to_csv(os.path.join(work_dir, "MSA_set_ranking.csv"), index=False)
    if seq_type == "codon":
        shutil.rmtree(os.path.join(work_dir, "MSA_set_translated"), ignore_errors=True)
    best = os.path.join(set_dir, ranking.loc[0, "code"])
    print(f"[msa_set] best MSA: {ranking.loc[0, 'code']} (predicted score {ranking.loc[0, 'predicted_score']:.6f})")
    return best
