"""
DNA DIPOLE CANARY
==================
Question: Does the info-layer dipole appear in DNA sequences?

DNA has two natural binary partitions per base — both are
"matched pairs that are different" in the dipole sense:
  R/Y axis:  Purine (A,G) vs Pyrimidine (T,C)   — size asymmetry
  S/W axis:  Strong (G,C)  vs Weak (A,T)        — bond-count asymmetry
                                                   (G-C = 3 H-bonds,
                                                    A-T = 2 H-bonds)

Mapping to Level 1 framework:
  H_R(p) = binary entropy of R-vs-Y in sliding window at position p
  H_S(p) = binary entropy of S-vs-W in sliding window at position p
  MI(p)  = mutual information between R/Y stream and S/W stream
           within the window

"Time" axis is now position along the sequence. We compute
dH_R/dp and run OD against the same kind of operator library
that found the Level 1 dipole.

Canary design (synthetic):
  Region A: uniform random DNA (null - no structure)
  Region B: codon-biased DNA (real coding has 3-periodic structure
            and amino-acid-driven base composition)
  Region C: tandem repeat (highly periodic)

Falsification-first:
  - Dipole absent in all three      -> wrong detector for sequences
  - Dipole in B and C but not A     -> dipole reads structure,
                                       methodology valid, move to chr22
  - Dipole in A (random)            -> methodology artifact, redesign
"""
import numpy as np
from numpy.linalg import lstsq
import time
import warnings
warnings.filterwarnings('ignore')

np.random.seed(7)
T_START = time.time()

# ============================================================
# SYNTHETIC DNA GENERATORS
# ============================================================
BASES = np.array(['A', 'C', 'G', 'T'])
B2I = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

def gen_random(n):
    """Uniform random DNA, no structure."""
    return np.random.choice(BASES, size=n)

def gen_codon_biased(n):
    """
    DNA with realistic codon-level bias. We pick codons from a weighted
    distribution that mimics human coding-region statistics: some
    amino acids are more common, GC3 (third codon position) skews high,
    stop codons rare.
    """
    n = (n // 3) * 3
    # 64 codons with rough human-coding weights (not exact, but realistic)
    codons = []
    weights = []
    for b1 in BASES:
        for b2 in BASES:
            for b3 in BASES:
                c = b1 + b2 + b3
                if c in ('TAA', 'TAG', 'TGA'):
                    w = 0.001  # stops, rare
                else:
                    # GC3 bias, slight purine preference at pos 1
                    w = 1.0
                    if b3 in ('G', 'C'):
                        w *= 1.6
                    if b1 in ('A', 'G'):
                        w *= 1.15
                    if b2 == 'A':
                        w *= 1.2
                codons.append(c)
                weights.append(w)
    weights = np.array(weights) / np.sum(weights)
    idx = np.random.choice(len(codons), size=n // 3, p=weights)
    seq = ''.join(codons[i] for i in idx)
    return np.array(list(seq))

def gen_tandem_repeat(n, motif='CAGCTG'):
    """Tandem repeat — extreme periodicity."""
    m = len(motif)
    reps = (n // m) + 1
    s = (motif * reps)[:n]
    return np.array(list(s))

# ============================================================
# SLIDING-WINDOW INFO OBSERVABLES
# ============================================================
def binary_entropy(p):
    """Binary entropy in nats; p in [0,1]."""
    if p <= 0 or p >= 1:
        return 0.0
    return -p * np.log(p) - (1 - p) * np.log(1 - p)

def sliding_info(seq, window=120, stride=4):
    """
    Returns H_R(p), H_S(p), MI(p) along the sequence.

    H_R: binary entropy of Purine(A,G) vs Pyrimidine(T,C)
    H_S: binary entropy of Strong(G,C) vs Weak(A,T)
    MI:  mutual information between R/Y and S/W binary streams
         INSIDE the window (4-cell joint distribution over RS,RW,YS,YW)
    """
    n = len(seq)
    is_R = np.isin(seq, ['A', 'G']).astype(np.float64)
    is_S = np.isin(seq, ['G', 'C']).astype(np.float64)
    starts = np.arange(0, n - window + 1, stride)
    H_R = np.zeros(len(starts))
    H_S = np.zeros(len(starts))
    MI = np.zeros(len(starts))
    pos = np.zeros(len(starts))

    for k, s in enumerate(starts):
        e = s + window
        rw = is_R[s:e]
        sw = is_S[s:e]
        p_R = rw.mean()
        p_S = sw.mean()
        H_R[k] = binary_entropy(p_R)
        H_S[k] = binary_entropy(p_S)

        # joint p(R=r, S=s) over 4 cells
        # G = R&S, A = R&W, C = Y&S, T = Y&W
        p_GR = (rw * sw).mean()                          # R=1,S=1  -> G
        p_AR = (rw * (1 - sw)).mean()                    # R=1,S=0  -> A
        p_CY = ((1 - rw) * sw).mean()                    # R=0,S=1  -> C
        p_TY = ((1 - rw) * (1 - sw)).mean()              # R=0,S=0  -> T
        joint = np.array([p_GR, p_AR, p_CY, p_TY])
        # MI = sum p(r,s) log [p(r,s) / (p(r) p(s))]
        marg_R = np.array([p_R, p_R, 1 - p_R, 1 - p_R])
        marg_S = np.array([p_S, 1 - p_S, p_S, 1 - p_S])
        mi = 0.0
        for q, mr, ms in zip(joint, marg_R, marg_S):
            if q > 1e-12 and mr > 1e-12 and ms > 1e-12:
                mi += q * np.log(q / (mr * ms))
        MI[k] = max(0.0, mi)
        pos[k] = s + window / 2.0
    return pos, H_R, H_S, MI

# ============================================================
# OD ON SEQUENCE-POSITION DERIVATIVES
# ============================================================
def od_dna(pos, H_R, H_S, MI, label):
    """Run OD: dH_R/dp = library @ coefficients."""
    dp = pos[1] - pos[0]
    dHR = np.gradient(H_R, dp)
    # operator library (mirrors the Level 1 library that found the dipole)
    lib = np.column_stack([
        H_R, H_S, H_R**2, H_S**2, H_R * H_S, MI, H_R * MI, H_S * MI,
        np.ones(len(H_R))
    ])
    labels = ['H_R', 'H_S', 'H_R^2', 'H_S^2', 'H_R*H_S', 'MI', 'H_R*MI', 'H_S*MI', 'const']
    c, *_ = lstsq(lib, dHR, rcond=None)
    pred = lib @ c
    ss_res = np.sum((dHR - pred) ** 2)
    ss_tot = np.sum((dHR - dHR.mean()) ** 2) + 1e-15
    r2 = 1 - ss_res / ss_tot

    # also try dH_S/dp
    dHS = np.gradient(H_S, dp)
    c_S, *_ = lstsq(lib, dHS, rcond=None)
    ss_res_S = np.sum((dHS - lib @ c_S) ** 2)
    ss_tot_S = np.sum((dHS - dHS.mean()) ** 2) + 1e-15
    r2_S = 1 - ss_res_S / ss_tot_S

    # and dMI/dp (this is where Level 1 dipole was strongest)
    dMI = np.gradient(MI, dp)
    c_MI, *_ = lstsq(lib, dMI, rcond=None)
    ss_res_MI = np.sum((dMI - lib @ c_MI) ** 2)
    ss_tot_MI = np.sum((dMI - dMI.mean()) ** 2) + 1e-15
    r2_MI = 1 - ss_res_MI / ss_tot_MI

    return {
        'dH_R/dp': {'r2': r2, 'coeffs': dict(zip(labels, c))},
        'dH_S/dp': {'r2': r2_S, 'coeffs': dict(zip(labels, c_S))},
        'dMI/dp': {'r2': r2_MI, 'coeffs': dict(zip(labels, c_MI))},
    }

def print_result(label, res):
    print(f"\n  {label}")
    print(f"  {'-' * len(label)}")
    for target, r in res.items():
        print(f"    {target:<10s}  R^2 = {r['r2']:+.3f}")
        # top 3 operators by |coefficient|
        items = [(k, v) for k, v in r['coeffs'].items() if k != 'const']
        items.sort(key=lambda x: -abs(x[1]))
        for op, c in items[:5]:
            print(f"      {op:<10s} {c:+.4e}")

# ============================================================
# RUN
# ============================================================
print("=" * 65)
print("  DNA DIPOLE CANARY")
print("=" * 65)
print("""
  Detector test: does the Level 1 dipole library find structure
  in DNA sequences? Synthetic canary across 3 region types.
""")

N = 60000  # 60 kb per region
WINDOW = 120
STRIDE = 4

regions = [
    ("RANDOM (null)",          gen_random(N)),
    ("CODON-BIASED (coding-like)", gen_codon_biased(N)),
    ("TANDEM REPEAT (structured)", gen_tandem_repeat(N, motif='CAGCTGAAGC')),
]

all_results = {}
for label, seq in regions:
    t0 = time.time()
    pos, H_R, H_S, MI = sliding_info(seq, window=WINDOW, stride=STRIDE)
    res = od_dna(pos, H_R, H_S, MI, label)
    elapsed = time.time() - t0
    print_result(f"{label}    [n_windows={len(pos)}, {elapsed:.1f}s]", res)

    # observable summary
    print(f"    H_R: mean={H_R.mean():.3f} std={H_R.std():.4f}")
    print(f"    H_S: mean={H_S.mean():.3f} std={H_S.std():.4f}")
    print(f"    MI:  mean={MI.mean():.4f} std={MI.std():.4f}")
    all_results[label] = res

# ============================================================
# DIPOLE TEST: does H_a^2 oppose H_a*H_b like in Level 1?
# ============================================================
print(f"\n{'=' * 65}")
print(f"  DIPOLE OPPOSITION TEST")
print(f"{'=' * 65}")
print(f"\n  Level 1 finding: H_a^2 and H_a*H_b are sign-opposite.")
print(f"  Does this hold in DNA?\n")
print(f"  Region                              dH_R/dp       dH_S/dp       dMI/dp")
print(f"  ----------------------------------  ------------  ------------  ------------")
for label, res in all_results.items():
    line = f"  {label:<34s}"
    for target in ['dH_R/dp', 'dH_S/dp', 'dMI/dp']:
        c = res[target]['coeffs']
        # check H_R^2 vs H_R*H_S signs (or H_S^2 vs H_R*H_S for dH_S)
        if target == 'dH_R/dp':
            self_term = c['H_R^2']
            cross = c['H_R*H_S']
        elif target == 'dH_S/dp':
            self_term = c['H_S^2']
            cross = c['H_R*H_S']
        else:
            # for dMI/dp, dipole-style would be H_R^2 vs H_R*H_S
            self_term = c['H_R^2']
            cross = c['H_R*H_S']
        opposing = (self_term * cross) < 0
        ratio = abs(self_term) / (abs(cross) + 1e-12)
        glyph = 'OPPOSE' if opposing else 'same'
        line += f"  {glyph:<6s} {ratio:5.2f} "
    print(line)

print(f"\n{'=' * 65}")
print(f"  TOTAL RUNTIME: {time.time() - T_START:.1f}s")
print(f"{'=' * 65}")
print("""
INTERPRETATION
==============
- Random region with R^2 > 0.3 -> methodology has bias, recheck
- Coding/Repeat with R^2 high + opposition -> detector works on
  sequences, dipole reads DNA structure, proceed to chr22
- All three weak -> wrong observable for sequences; try k-mer
  entropy, codon-position-stratified entropy, or strand asymmetry
""")
