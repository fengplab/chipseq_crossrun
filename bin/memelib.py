"""Minimal MEME motif format reader/writer (DNA) plus simple PWM descriptors."""
import math
import re

HEADER = """MEME version 5

ALPHABET= ACGT

strands: + -

Background letter frequencies
A {A:.3f} C {C:.3f} G {G:.3f} T {T:.3f}

"""


class Motif:
    def __init__(self, mid, alt, rows, nsites=None, evalue=None):
        self.id, self.alt, self.rows = mid, alt, rows
        self.nsites, self.evalue = nsites, evalue

    @property
    def width(self):
        return len(self.rows)

    def consensus(self):
        iupac = {"A": "A", "C": "C", "G": "G", "T": "T", "AC": "M", "AG": "R", "AT": "W",
                 "CG": "S", "CT": "Y", "GT": "K"}
        out = []
        for r in self.rows:
            order = sorted(range(4), key=lambda i: -r[i])
            if r[order[0]] >= 0.6:
                out.append("ACGT"[order[0]])
            elif r[order[0]] + r[order[1]] >= 0.75:
                out.append(iupac["".join(sorted("ACGT"[order[0]] + "ACGT"[order[1]]))])
            else:
                out.append("N")
        return "".join(out)

    def info_content(self):
        ic = 0.0
        for r in self.rows:
            ic += 2 + sum(p * math.log2(p) for p in r if p > 0)
        return ic

    def self_periodicity(self):
        """Best internal shift (1..w//2) by mean Pearson correlation between column vectors.
        High correlation at a small shift = the PWM is a tandem copy of a shorter unit."""
        best = (0, 0.0)
        w = self.width
        for k in range(1, max(1, w // 2) + 1):
            cors = []
            for i in range(w - k):
                a, b = self.rows[i], self.rows[i + k]
                ma, mb = sum(a) / 4, sum(b) / 4
                va = sum((x - ma) ** 2 for x in a)
                vb = sum((y - mb) ** 2 for y in b)
                if va > 0 and vb > 0:
                    cors.append(sum((x - ma) * (y - mb) for x, y in zip(a, b)) / math.sqrt(va * vb))
            if cors:
                m = sum(cors) / len(cors)
                if m > best[1]:
                    best = (k, m)
        return best

    def to_meme(self):
        lines = [f"MOTIF {self.id}" + (f" {self.alt}" if self.alt else ""), ""]
        ns = self.nsites if self.nsites else 20
        ev = self.evalue if self.evalue is not None else 0
        lines.append(f"letter-probability matrix: alength= 4 w= {self.width} nsites= {ns} E= {ev}")
        for r in self.rows:
            lines.append(" " + "  ".join(f"{p:.6f}" for p in r))
        lines.append("")
        return "\n".join(lines) + "\n"


def trimer_diversity(seq):
    s = re.sub(r"[^ACGT]", "", seq.upper())
    if len(s) < 3:
        return float("nan")
    k = {s[i:i + 3] for i in range(len(s) - 2)}
    return len(k) / (len(s) - 2)


def read_meme(path):
    """Return (background dict, list[Motif])."""
    bg = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}
    motifs = []
    with open(path) as fh:
        lines = fh.read().splitlines()
    i = 0
    while i < len(lines):
        ln = lines[i].strip()
        if ln.startswith("Background letter frequencies"):
            toks = []
            j = i + 1
            while j < len(lines) and lines[j].strip() and not lines[j].startswith("MOTIF"):
                toks += lines[j].split()
                j += 1
            for a, b in zip(toks[0::2], toks[1::2]):
                if a in bg:
                    try:
                        bg[a] = float(b)
                    except ValueError:
                        pass
            i = j
            continue
        if ln.startswith("MOTIF"):
            parts = ln.split()
            mid = parts[1]
            alt = parts[2] if len(parts) > 2 else ""
            j = i + 1
            while j < len(lines) and not lines[j].strip().startswith("letter-probability"):
                if lines[j].strip().startswith("MOTIF"):
                    break
                j += 1
            if j >= len(lines) or not lines[j].strip().startswith("letter-probability"):
                i = j
                continue
            hdr = lines[j]
            w = int(re.search(r"w=\s*(\d+)", hdr).group(1))
            ns = re.search(r"nsites=\s*([\d.eE+-]+)", hdr)
            ev = re.search(r"E=\s*([\d.eE+-]+)", hdr)
            rows = []
            k = j + 1
            while len(rows) < w and k < len(lines):
                vals = lines[k].split()
                if len(vals) >= 4:
                    try:
                        r = [float(v) for v in vals[:4]]
                        s = sum(r)
                        rows.append([v / s for v in r] if s > 0 else [0.25] * 4)
                    except ValueError:
                        pass
                k += 1
            motifs.append(Motif(mid, alt, rows,
                                float(ns.group(1)) if ns else None,
                                ev.group(1) if ev else None))
            i = k
            continue
        i += 1
    return bg, motifs


def write_meme(path, motifs, bg=None):
    bg = bg or {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}
    with open(path, "w") as fh:
        fh.write(HEADER.format(**bg))
        for m in motifs:
            fh.write(m.to_meme())
