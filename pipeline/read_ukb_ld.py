import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gcs_io import read_binary  # noqa: E402


def read_gctb_ldm_block(block_id, bin_dir, ldm_info_file, snp_info_file):
    """
    Read one block from a GCTB blockwise full LD matrix.

    Parameters
    ----------
    block_id      : int, 1-based block number (matches Block column in ldm.info)
    bin_dir       : directory containing block{N}.ldm.bin files
    ldm_info_file : path to ldm.info  (block registry)
    snp_info_file : path to snp.info  (SNP registry)
    """
    # --- Block registry ---
    ldm_info = pd.read_csv(ldm_info_file, sep=r'\s+')
    block_row = ldm_info[ldm_info['Block'] == block_id].squeeze()
    n_snps = int(block_row['NumSnps'])

    # --- SNP metadata for this block ---
    # Chunked + Block-filtered read: snp.info is ~13M rows / ~1.4 GB, too big to
    # materialize whole on memory-constrained Dataproc masters (e.g. n2-standard-2,
    # ~2 GB free after YARN/HDFS/Hive/Spark). We only need the ~n_snps rows for this
    # one block, so filter per-chunk and keep only the small matched subset.
    chunks = [c[c['Block'] == block_id] for c in pd.read_csv(snp_info_file, sep=r'\s+', chunksize=1_000_000)]
    snp_info = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
    assert len(snp_info) == n_snps, \
        f"Block {block_id}: ldm.info says {n_snps} SNPs, snp.info has {len(snp_info)}"

    # --- Binary --- (gs:// or local; read_binary is a plain np.fromfile-equivalent
    # for local paths and a gcsfs read for gs://)
    bin_file = f"{bin_dir}/block{block_id}.ldm.bin"
    raw = np.frombuffer(read_binary(bin_file), dtype=np.float32)

    n_upper = n_snps * (n_snps + 1) // 2
    n_full = n_snps * n_snps

    if len(raw) == n_upper:  # triangular storage (v2.5.5+): GCTB stores lower triangle
        ld = np.zeros((n_snps, n_snps), dtype=np.float32)
        ld[np.tril_indices(n_snps)] = raw
        ld += ld.T
        np.fill_diagonal(ld, np.diag(ld) / 2)  # diagonal was doubled by ld.T

    elif len(raw) == n_full:  # older: full symmetric
        ld = raw.reshape(n_snps, n_snps)

    else:
        raise ValueError(
            f"Block {block_id}: {len(raw)} floats for {n_snps} SNPs "
            f"(expected {n_upper} upper-tri or {n_full} full)"
        )

    return ld, snp_info


def extract_region_ld(chrom, start, end, bin_dir, ldm_info_file, snp_info_file):
    """
    Extract the LD (r²) matrix for SNPs in a genomic region.
    If the region spans multiple blocks, all overlapping blocks are loaded and
    the matrices are stitched (off-diagonal blocks are zeros — cross-block LD
    is not stored in GCTB blockwise LDMs).

    Returns
    -------
    r2  : np.ndarray (n, n) float32  — r² values
    snps: pd.DataFrame               — SNP metadata, sorted by PhysPos
    """
    snp_info_all = pd.read_csv(snp_info_file, sep=r'\s+')
    region_snps  = snp_info_all[
        (snp_info_all['Chrom'] == chrom) &
        (snp_info_all['PhysPos'] >= start) &
        (snp_info_all['PhysPos'] <= end)
    ].copy().reset_index(drop=True)

    if region_snps.empty:
        raise ValueError(f"No SNPs found in chr{chrom}:{start}-{end}")

    blocks_needed = sorted(region_snps['Block'].unique())
    print(f"Region chr{chrom}:{start}-{end} — {len(region_snps)} SNPs across block(s): {blocks_needed}")

    # Load each required block, keep only SNPs in region
    block_matrices = {}   # block_id -> (ld_sub, snps_sub)
    for bid in blocks_needed:
        ld_block, snps_block = read_gctb_ldm_block(bid, bin_dir, ldm_info_file, snp_info_file)
        mask      = snps_block['PhysPos'].between(start, end)
        idx       = np.where(mask)[0]
        block_matrices[bid] = (ld_block[np.ix_(idx, idx)], snps_block[mask].reset_index(drop=True))

    if len(blocks_needed) == 1:
        r2_full = block_matrices[blocks_needed[0]][0] ** 2
        snps    = block_matrices[blocks_needed[0]][1]
    else:
        # Stitch: block-diagonal only (cross-block LD not available)
        snps   = pd.concat([block_matrices[b][1] for b in blocks_needed], ignore_index=True)
        n      = len(snps)
        r2_full = np.zeros((n, n), dtype=np.float32)
        offset = 0
        for bid in blocks_needed:
            sub = block_matrices[bid][0] ** 2
            k   = sub.shape[0]
            r2_full[offset:offset+k, offset:offset+k] = sub
            offset += k
        np.fill_diagonal(r2_full, 1.0)

    return r2_full, snps


def plot_haploview(r2, snps, label_col='ID', max_labels=50,
                   cmap='Reds', title=None, out_file=None):
    """
    Haploview-style triangle LD plot using imshow on a rotated array.
    Scales to thousands of SNPs — single image draw call regardless of n.

    The triangle opens downward: diagonal (r²=1) at the top,
    most-distant pairs at the bottom apex.

    Parameters
    ----------
    r2         : (n, n) r² matrix
    snps       : DataFrame with SNP metadata (must align row/col order with r2)
    label_col  : column in snps to use for x-axis tick labels
    max_labels : thin labels if there are more SNPs than this
    cmap       : matplotlib colormap
    title      : optional plot title
    out_file   : if given, save figure to this path instead of showing
    """
    n = r2.shape[0]

    # Rotate upper triangle into a (n × 2n-1) image.
    # Cell (i, j) with j>=i maps to: row = j-i, col = i+j
    rotated = np.full((n, 2 * n - 1), np.nan, dtype=np.float32)
    for d in range(n):
        diag = np.diagonal(r2, offset=d)         # length n-d
        cols = np.arange(len(diag)) * 2 + d
        rotated[d, cols] = diag

    fig_w = max(10, n * 0.05)
    fig_h = max(5,  n * 0.025)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    im = ax.imshow(rotated, aspect='auto', cmap=cmap, vmin=0, vmax=1,
                   interpolation='nearest', origin='upper')

    ax.set_yticks([])

    # X-axis ticks: SNP i sits at column 2*i in the rotated image
    step        = max(1, n // max_labels)
    tick_idx    = list(range(0, n, step))
    tick_cols   = [2 * i for i in tick_idx]
    tick_labels = snps[label_col].iloc[tick_idx].tolist()
    ax.set_xticks(tick_cols)
    ax.set_xticklabels(tick_labels, rotation=90, fontsize=7)
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')

    plt.colorbar(im, ax=ax, label='r²', fraction=0.02, pad=0.01)
    if title:
        ax.set_title(title, pad=30)

    plt.tight_layout()
    if out_file:
        plt.savefig(out_file, dpi=150, bbox_inches='tight')
        print(f"Saved: {out_file}")
    else:
        plt.show()
    plt.close(fig)


if __name__ == '__main__':
    # Example: one locus from the unpacked GCTB ukbEUR_13M_FullLDM (set GCTB_LDM to its directory).
    UKB = os.path.expanduser(os.environ.get("GCTB_LDM", "~/work/ukbEUR_13M_FullLDM"))
    r2, snps = extract_region_ld(
        chrom=2, start=57771591, end=58946078,
        bin_dir=UKB,
        ldm_info_file=f"{UKB}/ldm.info",
        snp_info_file=f"{UKB}/snp.info",
    )
    plot_haploview(r2, snps, label_col='ID', title="chr2:57771591-58946078", out_file="ld_triangle_chr2.png")


