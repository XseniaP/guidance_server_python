import os
import tarfile
import zipfile

import pytest

from guidance3.sequences.msa_set import MSASetError, collect_msa_files, validate_msa_set

A = ">s1\nAC-GT\n>s2\nACG-T\n>s3\nA-CGT\n>s4\nACGT-\n"
B = ">s1\nACGT-\n>s2\nACG-T\n>s3\nAC-GT\n>s4\n-ACGT\n"


def _write(d, name, text):
    p = os.path.join(d, name)
    with open(p, "w") as fh:
        fh.write(text)
    return p


def test_folder_ignores_non_msa_files(tmp_path):
    src = tmp_path / "src"
    src.mkdir()
    _write(src, "a.fasta", A)
    _write(src, "b.fas", B)
    _write(src, "notes.txt", "x")
    files = collect_msa_files(str(src), str(tmp_path / "out"))
    assert [os.path.basename(f) for f in files] == ["a.fasta", "b.fas"]
    validate_msa_set(files)


def test_zip_and_targz(tmp_path):
    src = tmp_path / "src"
    src.mkdir()
    a, b = _write(src, "a.fasta", A), _write(src, "b.fasta", B)
    z = tmp_path / "set.zip"
    with zipfile.ZipFile(z, "w") as zf:
        zf.write(a, "sub/a.fasta")
        zf.write(b, "b.fasta")
        zf.writestr("__MACOSX/._a.fasta", "junk")
    assert len(collect_msa_files(str(z), str(tmp_path / "z"))) == 2
    t = tmp_path / "set.tar.gz"
    with tarfile.open(t, "w:gz") as tf:
        tf.add(a, "a.fasta")
        tf.add(b, "b.fasta")
    assert len(collect_msa_files(str(t), str(tmp_path / "t"))) == 2


def test_zip_path_traversal_is_skipped(tmp_path):
    z = tmp_path / "bad.zip"
    with zipfile.ZipFile(z, "w") as zf:
        zf.writestr("../evil.fasta", A)
        zf.writestr("ok.fasta", A)
    files = collect_msa_files(str(z), str(tmp_path / "out"))
    assert [os.path.basename(f) for f in files] == ["ok.fasta"]
    assert not (tmp_path / "evil.fasta").exists()


@pytest.mark.parametrize("second, message", [
    (">s1\nACGT\n>s2\nACG-T\n>s3\nA-CGT\n>s4\nACGT-\n", "not aligned"),
    (">s1\nAC-GT\n>s2\nACG-T\n>s3\nA-CGT\n>sX\nACGT-\n", "same sequence names"),
    (">s1\nAC-GA\n>s2\nACG-T\n>s3\nA-CGT\n>s4\nACGT-\n", "different sequences"),
    (">s1\nAC-GT\n>s2\nACG-T\n>s3\nA-CGT\n", "at least 4"),
])
def test_invalid_sets(tmp_path, second, message):
    files = [_write(tmp_path, "a.fasta", A), _write(tmp_path, "b.fasta", second)]
    with pytest.raises(MSASetError, match=message):
        validate_msa_set(files)


def test_single_msa_is_rejected(tmp_path):
    with pytest.raises(MSASetError, match="at least 2"):
        validate_msa_set([_write(tmp_path, "a.fasta", A)])
