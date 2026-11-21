from __future__ import annotations

import argparse
from collections import Counter
import io
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact


SAMPLE_TARGET_FASTA = """>target_seq_1
ATGGCCGGTTAAATGGCCGGTTAA
>target_seq_2
CCGGTTAACCGGTTAACCGGTTAA
>target_seq_3
ATGGCCGGTTAA
"""

SAMPLE_BG_FASTA = """>bg_seq_1
ATATATATATATATATATATATAT
>bg_seq_2
GCGCGCGCGCGCGCGCGCGCGCGC
>bg_seq_3
TATATATATATATATATATATATA
"""


def read_fasta(fasta_path: str | None, fallback_text: str) -> list[str]:
    """Reads sequences from a FASTA file (or built-in sample) and returns uppercase sequences."""
    raw_text = Path(fasta_path).read_text(encoding="utf-8") if fasta_path else fallback_text
    sequences: list[str] = []
    current_seq: list[str] = []

    for line in io.StringIO(raw_text):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_seq:
                sequences.append("".join(current_seq).upper())
                current_seq = []
        else:
            current_seq.append(line)

    if current_seq:
        sequences.append("".join(current_seq).upper())

    return sequences


def count_kmers(sequences: list[str], k: int) -> Counter[str]:
    """Counts all non-N k-mers across a list of sequences."""
    counts: Counter[str] = Counter()
    for seq in sequences:
        for i in range(len(seq) - k + 1):
            kmer = seq[i : i + k]
            # 'N' (bilinmeyen baz) içeren k-mer'leri dışarıda bırakmak temiz analiz sağlar
            if "N" not in kmer:
                counts[kmer] += 1
    return counts


def analyze_enrichment(
    target_counts: Counter[str],
    bg_counts: Counter[str],
    min_target_count: int = 2,
) -> pd.DataFrame:
    """Compares k-mer counts between target and background using Fisher's exact test."""
    target_total = sum(target_counts.values())
    bg_total = sum(bg_counts.values())
    all_kmers = set(target_counts.keys()) | set(bg_counts.keys())

    rows: list[dict[str, object]] = []
    for kmer in all_kmers:
        t_in = target_counts.get(kmer, 0)
        b_in = bg_counts.get(kmer, 0)

        # Çok nadir (örneğin sadece 1 kez görülen) gürültü k-mer'leri elemek için eşik
        if t_in < min_target_count:
            continue

        t_out = max(target_total - t_in, 0)
        b_out = max(bg_total - b_in, 0)

        # 2x2 kontenjans tablosu: hedef içi/dışı vs. arka plan içi/dışı
        _, p_val = fisher_exact([[t_in, t_out], [b_in, b_out]], alternative="greater")

        t_freq = t_in / target_total if target_total else 0.0
        b_freq = b_in / bg_total if bg_total else 0.0
        # Sıfıra bölünmeyi önlemek için küçük pseudocount ekliyoruz
        enrichment_ratio = (t_freq + 1e-6) / (b_freq + 1e-6)

        rows.append(
            {
                "kmer": kmer,
                "target_count": t_in,
                "bg_count": b_in,
                "target_freq": round(t_freq, 4),
                "bg_freq": round(b_freq, 4),
                "enrichment_ratio": round(enrichment_ratio, 2),
                "p_value": float(p_val),
            }
        )

    if not rows:
        return pd.DataFrame(
            columns=[
                "kmer",
                "target_count",
                "bg_count",
                "target_freq",
                "bg_freq",
                "enrichment_ratio",
                "p_value",
            ]
        )

    df = pd.DataFrame(rows)
    return df.sort_values(["p_value", "enrichment_ratio"], ascending=[True, False]).reset_index(
        drop=True
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="K-mer enrichment analysis between two FASTA files.")
    parser.add_argument("--target", help="Target FASTA path (uses built-in sample if omitted)")
    parser.add_argument("--background", help="Background FASTA path (uses built-in sample if omitted)")
    parser.add_argument("-k", type=int, default=4, help="K-mer length (default: 4)")
    parser.add_argument("--min-count", type=int, default=2, help="Minimum target count (default: 2)")
    parser.add_argument("--out", default="kmer_enrichment.csv", help="Output CSV path")
    args = parser.parse_args()

    target_seqs = read_fasta(args.target, SAMPLE_TARGET_FASTA)
    bg_seqs = read_fasta(args.background, SAMPLE_BG_FASTA)

    target_counts = count_kmers(target_seqs, args.k)
    bg_counts = count_kmers(bg_seqs, args.k)
    result = analyze_enrichment(target_counts, bg_counts, min_target_count=args.min_count)

    print(result.to_string(index=False))
    result.to_csv(args.out, index=False)
    print(f"\nSaved to {args.out}")


if __name__ == "__main__":
    main()
