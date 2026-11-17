"""gs://-aware file I/O helpers shared by the EUR GCTB reference-LD scripts.

These scripts were originally written for local paths (they ran on randi). Inside the
AoU perimeter there's no "randi" to egress results to (see multi_ancestry/aou_ld/README.md),
so they instead read/write gs:// paths directly. Local paths still behave exactly as
before -- every function here is a no-op passthrough to the plain stdlib call for
anything not starting with "gs://".
"""
import gzip
import io
import os


def is_gcs(path):
    return str(path).startswith("gs://")


def _fs():
    import gcsfs
    return gcsfs.GCSFileSystem()


def exists(path):
    if is_gcs(path):
        return _fs().exists(path)
    return os.path.exists(path)


def makedirs(path):
    if not is_gcs(path):
        os.makedirs(path, exist_ok=True)


def open_binary(path, mode="rb"):
    """Binary file handle for a gs:// or local path."""
    if is_gcs(path):
        return _fs().open(path, mode)
    return open(path, mode)


def open_maybe_gz(path, mode="rt", **kwargs):
    """Text (or binary) handle for a local or gs:// path, transparently gzipped if the
    path ends in .gz. For local paths this is byte-for-byte what the plain
    gzip.open(path)/open(path) call it replaces would have done. Extra kwargs (e.g.
    encoding=, errors=) are passed through to gzip.open/open, same as before."""
    if not is_gcs(path):
        return gzip.open(path, mode, **kwargs) if path.endswith(".gz") else open(path, mode, **kwargs)
    binmode = "rb" if "r" in mode else "wb"
    fh = _fs().open(path, binmode)
    if path.endswith(".gz"):
        return gzip.open(fh, mode, **kwargs)
    if "b" in mode:
        return fh
    # Text mode on a non-gzipped gs:// path: gcsfs hands back a binary handle, so wrap
    # it the same way a plain (non-.gz) local text-mode open() would.
    return io.TextIOWrapper(fh, **kwargs)


def read_binary(path):
    """All bytes of a local or gs:// file (for np.frombuffer-style binary reads)."""
    with open_binary(path, "rb") as f:
        return f.read()
