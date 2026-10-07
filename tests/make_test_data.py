"""Build a small synthetic genome + fake pipeline outputs (10 samples, 3 runs, one planted satellite
time-course effect) for exercising the bundle. Usage: python3 tests/make_test_data.py <dir>"""
import random, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
from memelib import Motif, write_meme
random.seed(7)
R = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else "test_data")
os.makedirs(R, exist_ok=True)
L = 300000
chroms = ["chr1", "chr2", "chr3", "chr4"]
seq = {c: [random.choice("ACGT") for _ in range(L)] for c in chroms}
def put(c, pos, s):
    seq[c][pos:pos+len(s)] = list(s)
# satellites
for p in range(100000, 130000, 5): put("chr1", p, "ATTCC")
cenpb = "TTCGTTGGAAACGGGA"
for p in range(50000, 70000, 171):
    put("chr2", p, "".join(random.choice("ACGT") for _ in range(171))); put("chr2", p+20, cenpb)
tf_sites = []
for i in range(400):
    c = random.choice(chroms); p = random.randint(1000, L-1000)
    if c == "chr1" and 95000 < p < 135000: continue
    if c == "chr2" and 45000 < p < 75000: continue
    put(c, p, "CCGGAAGT"); tf_sites.append((c, p))
alus = []
for i in range(200):
    c = random.choice(chroms); p = random.randint(1000, L-1000)
    if (c == "chr1" and 95000 < p < 135000) or (c == "chr2" and 45000 < p < 75000): continue
    alus.append((c, p, p+300))
with open(f"{R}/genome.fa", "w") as fh, open(f"{R}/genome.fa.fai", "w") as fi:
    off = 0
    for c in chroms:
        h = f">{c}\n"; fh.write(h); off += len(h)
        s = "".join(seq[c])
        fi.write(f"{c}\t{L}\t{off}\t60\t61\n")
        for i in range(0, L, 60):
            fh.write(s[i:i+60] + "\n"); off += len(s[i:i+60]) + 1
with open(f"{R}/genome.sizes", "w") as fh:
    for c in chroms: fh.write(f"{c}\t{L}\n")
with open(f"{R}/rmsk.bed", "w") as fh:
    for c, a, b in sorted(alus): fh.write(f"{c}\t{a}\t{b}\tAluY\t0\t+\tSINE\tAlu\n")
    fh.write("chr1\t100000\t130000\t(ATTCC)n\t0\t+\tSimple_repeat\tSimple_repeat\n")
    fh.write("chr2\t50000\t70000\tALR/Alpha\t0\t+\tSatellite\tcentr\n")
with open(f"{R}/censat.bed", "w") as fh:
    fh.write("chr1\t95000\t100000\tct_1_1(p_arm)\n")
    fh.write("chr1\t100000\t130000\thsat3_1(A1)\n")
    fh.write("chr2\t50000\t70000\tactive_hor(S1C2H1L)\n")
    fh.write("chr2\t70000\t75000\tmon_2\n")
with open(f"{R}/include.bed", "w") as fh:
    for c in chroms: fh.write(f"{c}\t0\t{L}\n")
with open(f"{R}/genes.bed", "w") as fh:
    for i, (c, p) in enumerate(sorted(tf_sites)[::4]): fh.write(f"{c}\t{p+100}\t{p+5000}\tGENE{i}\t0\t+\n")
def pwm(cons, p=0.85):
    rows = []
    for ch in cons:
        r = [(1-p)/3]*4; r["ACGT".index(ch)] = p; rows.append(r)
    return rows
known = [Motif("MA0637.1", "CENPB", pwm(cenpb)), Motif("MA0076.2", "ELK4", pwm("CCGGAAGT")),
         Motif("MA0079.5", "SP1", pwm("GGGGCGGGG")), Motif("MA9999.1", "FAKE1", pwm("ACGTACGTAAC"))]
write_meme(f"{R}/known.meme", known)
samples = [  # id, run, series, virus, cell, hpi, construct, sat_frac, tf_frac
 ("HSV_THP1_WT_0h_REP1","run01","HSV_THP1_WT","HSV-1","THP-1",0,"WT",0.05,0.4),
 ("HSV_THP1_WT_6h_REP1","run01","HSV_THP1_WT","HSV-1","THP-1",6,"WT",0.15,0.4),
 ("HSV_THP1_WT_12h_REP1","run01","HSV_THP1_WT","HSV-1","THP-1",12,"WT",0.25,0.35),
 ("HSV_THP1_WT_24h_REP1","run01","HSV_THP1_WT","HSV-1","THP-1",24,"WT",0.35,0.3),
 ("HSV_THP1_dN_0h_REP1","run02","HSV_THP1_dN","HSV-1","THP-1",0,"dN",0.03,0.4),
 ("HSV_THP1_dN_12h_REP1","run02","HSV_THP1_dN","HSV-1","THP-1",12,"dN",0.04,0.4),
 ("HSV_THP1_dN_24h_REP1","run02","HSV_THP1_dN","HSV-1","THP-1",24,"dN",0.05,0.4),
 ("VACV_HEK_WT_24h_REP1","run03","","VACV","HEK293T",24,"WT",0.10,0.6),
 ("mock_HEK_WT_REP1","run03","","none","HEK293T",0,"WT",0.02,0.6),
 ("VACV_HEK_dN_24h_REP1","run03","","VACV","HEK293T",24,"dN",0.02,0.6),
]
man = open(f"{R}/manifest.tsv", "w")
man.write("sample_id\trun_id\tresults_dir\tseries_id\tvirus\tcell_type\thpi\tcgas_construct\tcontrol_id\tgenome_build\n")
hot = [("chr3", 150000), ("chr4", 42000)]
for sid, run, ser, vir, cell, hpi, con, sf, tfr in samples:
    res = f"{R}/runs/{run}/results"
    mac = f"{res}/bowtie2/merged_library/macs3/narrow_peak"
    os.makedirs(f"{mac}/qc", exist_ok=True); os.makedirs(f"{res}/bowtie2/merged_library/samtools_stats", exist_ok=True)
    man.write(f"{sid}\t{run}\t{res}\t{ser}\t{vir}\t{cell}\t{hpi}\t{con}\tINPUT_{run}\ttest\n")
    n = random.randint(200, 400); peaks = []
    for i in range(n):
        u = random.random(); w = random.randint(200, 500)
        if u < sf:
            if random.random() < 0.5:   # HSat3 (chr1:100000-130000): peaks near the array edges
                c, p = "chr1", random.choice([random.randint(100000, 100800), random.randint(129200, 129600)])
            else:                       # alpha-satellite HOR: spread uniformly
                c, p = "chr2", random.randint(50500, 69000)
        elif u < sf + tfr:
            c, p = random.choice(tf_sites); p = p - random.randint(20, w-40)
        elif u < sf + tfr + 0.1:
            c, p = random.choice(hot); p += random.randint(-50, 50)
        else:
            c = random.choice(chroms); p = random.randint(1000, L - 1000)
        p = max(0, p)
        peaks.append((c, p, p + w))
    peaks.sort()
    with open(f"{mac}/{sid}_peaks.narrowPeak", "w") as fh:
        for i, (c, a, b) in enumerate(peaks):
            fe = round(random.lognormvariate(1.4, 0.5), 3)
            fh.write(f"{c}\t{a}\t{b}\t{sid}_peak_{i+1}\t{min(1000, int(fe * 30))}\t.\t{fe}\t{fe * 3:.3f}\t{fe * 2:.3f}\t{(b-a)//2}\n")
    with open(f"{mac}/{sid}_peaks.annotatePeaks.txt", "w") as fh:
        fh.write("PeakID (cmd=annotatePeaks.pl)\tChr\tStart\tEnd\tStrand\tPeak Score\tFocus Ratio/Region Size\tAnnotation\tDetailed Annotation\tDistance to TSS\n")
        for i, (c, a, b) in enumerate(peaks):
            an = random.choice(["promoter-TSS (NM_1)", "intron (NM_2, intron 1 of 3)", "Intergenic", "exon (NM_3)", "TTS (NM_4)"])
            fh.write(f"{sid}_peak_{i+1}\t{c}\t{a}\t{b}\t+\t100\tNA\t{an}\tx\t100\n")
    open(f"{mac}/qc/{sid}_peaks.FRiP_mqc.tsv", "w").write(f"# id: 'frip'\n{sid}\t{random.uniform(0.01,0.06):.4f}\n")
    open(f"{res}/bowtie2/merged_library/samtools_stats/{sid}.mLb.clN.sorted.bam.flagstat", "w").write(
        f"{random.randint(10**7, 4*10**7)} + 0 in total\n{random.randint(10**7, 3*10**7)} + 0 mapped (98.00% : N/A)\n")
    xd = f"{mac}/xstreme/{sid}_xstreme"; os.makedirs(f"{xd}/sea_out", exist_ok=True)
    mots = []
    if sf >= 0.1: mots.append(("1-ATTCCATTCCATTCC", "STREME-1", "ATTCCATTCCATTCC", "STREME", 1e-50))
    if sf >= 0.15: mots.append(("TTCGTTGGAAACGGGA", "MEME-1", cenpb, "MEME", 1e-20))
    if tfr >= 0.3: mots.append(("2-CCGGAAGTN" if random.random() < 0.5 else "3-ACCGGAAGT", "STREME-2", "CCGGAAGT" if random.random() < 0.5 else "ACCGGAAGT", "STREME", 1e-15))
    mots.append((f"4-{''.join(random.choice('ACGT') for _ in range(8))}", "STREME-4", "".join(random.choice("ACGT") for _ in range(8)), "STREME", 1e-3))
    write_meme(f"{xd}/combined.meme", [Motif(i, a, pwm(c)) for i, a, c, s, e in mots])
    with open(f"{xd}/xstreme.tsv", "w") as fh:
        fh.write("SEED_MOTIF\tRANK\tSOURCE\tID\tALT_ID\tCONSENSUS\tWIDTH\tSITES\tSEA_PVALUE\tEVALUE\tEVALUE_ACC\tSIM_SOURCE\tSIM_MOTIF\tMOTIF_URL\n")
        for r, (i, a, c, s, e) in enumerate(mots):
            fh.write(f"1\t{r+1}\t{s}\t{i}\t{a}\t{c}\t{len(c)}\t{random.randint(20,200)}\t{e}\t{e*10}\tyes\t\t\t\n")
        fh.write("# XSTREME stub\n")
    with open(f"{xd}/sea_out/sea.tsv", "w") as fh:
        fh.write("RANK\tDB\tID\tALT_ID\tCONSENSUS\tTP\tTP%\tFP\tFP%\tENR_RATIO\tSCORE_THR\tPVALUE\tLOG_PVALUE\tEVALUE\tLOG_EVALUE\tQVALUE\tLOG_QVALUE\n")
        for r, (i, a, c, s, e) in enumerate(mots):
            fh.write(f"{r+1}\tprimary\t{i}\t{a}\t{c}\t50\t20\t5\t2\t10\t5\t{e}\t0\t{e}\t0\t{e}\t0\n")
    open(f"{xd}/pipeline_maxw.txt", "w").write("maxw_used\t30\n")
man.close()

# ---- paired-end BAMs (IP + one control per run) for motif-discovery controls and insert-size QC ----
import shutil, subprocess
if shutil.which("samtools"):
    rc = str.maketrans("ACGT", "TGCA")
    genome = {c: "".join(seq[c]) for c in chroms}
    def write_bam(path, frags, rl=50):
        sam = path[:-4] + ".sam"
        with open(sam, "w") as fh:
            fh.write("@HD\tVN:1.6\tSO:unsorted\n")
            for c in chroms: fh.write(f"@SQ\tSN:{c}\tLN:{L}\n")
            for i, frag in enumerate(frags):
                c, st, fl = frag[:3]; fwd = frag[3] if len(frag) > 3 else (random.random() < 0.5)
                if st < 0 or st + fl > L or fl < rl: continue
                left = genome[c][st:st + rl]; right = genome[c][st + fl - rl:st + fl]
                q = "I" * rl
                if fwd:   # read 1 on the forward strand at the fragment's left end
                    fh.write(f"f{i}\t99\t{c}\t{st+1}\t60\t{rl}M\t=\t{st+fl-rl+1}\t{fl}\t{left}\t{q}\n")
                    fh.write(f"f{i}\t147\t{c}\t{st+fl-rl+1}\t60\t{rl}M\t=\t{st+1}\t{-fl}\t{right}\t{q}\n")
                else:     # read 1 on the reverse strand at the fragment's right end
                    fh.write(f"f{i}\t83\t{c}\t{st+fl-rl+1}\t60\t{rl}M\t=\t{st+1}\t{-fl}\t{right}\t{q}\n")
                    fh.write(f"f{i}\t163\t{c}\t{st+1}\t60\t{rl}M\t=\t{st+fl-rl+1}\t{fl}\t{left}\t{q}\n")
        subprocess.run(["samtools", "sort", "-o", path, sam], check=True, stderr=subprocess.DEVNULL)
        subprocess.run(["samtools", "index", path], check=True)
        os.remove(sam)
    runs_done = set()
    for sid, run, ser, vir, cell, hpi, con, sf, tfr in samples:
        ml = f"{R}/runs/{run}/results/bowtie2/merged_library"
        frags = []
        for i in range(6000):
            u = random.random()
            if u < sf + 0.05:   # IP over HSat3 / alpha-sat (planted QC signals)
                if random.random() < 0.5:   # HSat3: short fragments, stacked on few starts near the peak positions, one strand
                    st = random.choice([100050, 100300, 100600, 129300, 129450]) + random.choice([0, 0, 0, 7, 13])
                    frags.append(("chr1", st, int(random.gauss(150, 25)), True)); continue
                c, st = "chr2", random.randint(50000, 69500)
                fl = int(random.gauss(200, 25))
            else:
                c, st, fl = random.choice(chroms), random.randint(0, L - 600), int(random.gauss(200, 30))
            frags.append((c, st, fl))
        write_bam(f"{ml}/{sid}.mLb.clN.sorted.bam", frags)
        if run not in runs_done:
            runs_done.add(run)
            cf = [(random.choice(chroms), random.randint(0, L - 600), int(random.gauss(200, 30))) for _ in range(9000)]
            write_bam(f"{ml}/INPUT_{run}_REP1.mLb.clN.sorted.bam", cf)
print("ok")
