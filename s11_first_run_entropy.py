"""
S11 follow-up (Greg: "analyze the first run with that frame too").
The FIRST real-data run = GW150914 (Session 10, INFO-036). Re-read it
through the INFO-038 entropy-asymmetry frame, and answer the question the
3-event table raised:

  GW150914's event-window entropy asymmetry |H_a-H_b| BARELY moved
  (noise 0.637 -> event 0.620) yet its cos-to-attractor COLLAPSED
  (noise 0.980 -> event 0.231). So what breaks the (-1,-1,+2) attractor
  at the merger -- the marginal entropies (H_a,H_b), or the MI?

Frame-consistent probe: recompute the null direction with the MI column
REMOVED (5-op basis [H_a,H_b,H_a^2,H_b^2,H_a*H_b]) vs the full 6-op basis.
If dropping MI RESTORES the event-block attractor, the merger departure is
MI-driven (shared astrophysical signal), and the per-event |H_a-H_b|
spread (INFO-038) is a SEPARATE, noise-floor mechanism. Two distinct
things, not one.

Speaking posture (Rule C): I THINK the merger departure is MI-driven while
the cross-event asymmetry spread is noise-floor-driven -- two mechanisms --
but we wait on the numbers. Self-contained (kbk pure funcs only).
"""
import json, numpy as np, h5py
from scipy.signal import butter, filtfilt, welch, correlate
from scipy.interpolate import interp1d
from kbk_pipeline import compute_operator_matrix, extract_v1, project_234

FS = 4096; EDGE_S = 2.0; MT = 16.4
ATTR3 = np.array([-1., -1., 2.])/np.sqrt(6.)
H1P = "data/ligo/H-H1_GW150914_32s.hdf5"
L1P = "data/ligo/L-L1_GW150914_32s.hdf5"


def load(p):
    with h5py.File(p, "r") as f: return f["strain/Strain"][:].astype(float)
def bp(x, lo=35., hi=350., o=4):
    b, a = butter(o, [lo/(FS/2), hi/(FS/2)], btype="band"); return filtfilt(b, a, x)
def whiten(x):
    n=len(x); fr, ps = welch(x, fs=FS, nperseg=min(4*FS, n))
    ip=interp1d(fr, ps, bounds_error=False, fill_value=(ps[0], ps[-1]))
    X=np.fft.rfft(x); f=np.fft.rfftfreq(n, 1./FS)
    xw=np.fft.irfft(X/np.sqrt(np.maximum(ip(f), 1e-50)), n=n); return xw/np.std(xw)
def pre(x): return bp(whiten(bp(x)))


def cos_full(M):          # 6-op null projected to (Ha^2,Hb^2,HaHb)
    v, *_ = extract_v1(M); _, c = project_234(v); return abs(float(c))
def cos_noMI(M):          # 5-op null (MI removed), same projection
    M5 = M[:, :5]
    M5c = M5 - M5.mean(0)
    _, _, Vt = np.linalg.svd(M5c, full_matrices=False)
    v5 = Vt[-1]                       # smallest right-singular vector
    sub = v5[2:5]; n = np.linalg.norm(sub)
    return abs(float(sub @ ATTR3 / n)) if n > 0 else 0.0


h = pre(load(H1P)); l = pre(load(L1P))
c = int(MT*FS); half = int(0.25*FS)
xc = correlate(h[c-half:c+half], l[c-half:c+half], mode="full")
lags = np.arange(-half*2+1, half*2); m = np.abs(lags) <= int(0.015*FS)
j = np.argmax(np.abs(xc[m])); la = np.roll(l, lags[m][j])*np.sign(xc[m][j])
crop = int(EDGE_S*FS); a, b = h[crop:-crop], la[crop:-crop]
win, stride = int(0.125*FS), int(0.03125*FS)
M = compute_operator_matrix(a, b, win, stride, mi_bins=16)
centers = (np.arange(M.shape[0])*stride + win/2)/FS + EDGE_S
ev = np.abs(centers - MT) <= 0.5; no = ~ev

print("="*92)
print("FIRST RUN (GW150914) re-read in the INFO-038 entropy frame")
print("="*92)
print(f"{'block':6s} {'H_a':>7s} {'H_b':>7s} {'|Ha-Hb|':>8s} {'MI':>7s} "
      f"{'cos(6op)':>9s} {'cos(noMI)':>10s}")
out = {}
for lbl, mask in [("noise", no), ("event", ev)]:
    Ms = M[mask]
    row = dict(H_a=float(Ms[:,0].mean()), H_b=float(Ms[:,1].mean()),
               asym=float(np.abs(Ms[:,0]-Ms[:,1]).mean()), MI=float(Ms[:,5].mean()),
               cos6=cos_full(Ms), cosNoMI=cos_noMI(Ms))
    out[lbl] = row
    print(f"{lbl:6s} {row['H_a']:7.3f} {row['H_b']:7.3f} {row['asym']:8.4f} "
          f"{row['MI']:7.4f} {row['cos6']:9.3f} {row['cosNoMI']:10.3f}")

print("\nReadout:")
d_asym = out['event']['asym'] - out['noise']['asym']
d_mi   = out['event']['MI']  - out['noise']['MI']
print(f"  entropy asymmetry change at merger: {d_asym:+.4f}  (MI change: {d_mi:+.4f})")
print(f"  cos drop WITH MI    : {out['noise']['cos6']:.3f} -> {out['event']['cos6']:.3f}")
print(f"  cos drop WITHOUT MI : {out['noise']['cosNoMI']:.3f} -> {out['event']['cosNoMI']:.3f}")
restored = out['event']['cosNoMI'] - out['event']['cos6']
print(f"  removing MI changes event cos by {restored:+.3f} "
      f"({'RESTORES attractor -> merger departure is MI-driven' if restored>0.3 else 'little change -> not purely MI'})")

# merger trajectory: per-window asym + MI across +/-0.5s
print("\nper-window trajectory around merger (t, |Ha-Hb|, MI):")
band = np.abs(centers - MT) <= 0.30
for t, ha, hb, mi in zip(centers[band], M[band,0], M[band,1], M[band,5]):
    bar = "#"*int(mi*40)
    print(f"  t={t:6.3f}  |Ha-Hb|={abs(ha-hb):6.3f}  MI={mi:6.3f} {bar}")

json.dump(out, open("s11_first_run_entropy_results.json","w"), indent=2)
print("\nWrote s11_first_run_entropy_results.json")
