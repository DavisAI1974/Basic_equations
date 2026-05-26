"""Ingest all prior-session data into the OD learning store.

Three new problems populated from existing JSONs (no new data runs):

  Problem 1: cardiac_disease_family
    Source: ekg_disease_canary.json + ekg_disease_canary_expand.json
    10 single records, 10 classes (N=1/class) - SEED for future cohort
    Schema: same as cerebro (HR/HRV pair features + temporal)

  Problem 2: simulator_4domain
    Source: per_domain_kbk_results_seed11.json + seed22.json
    4 classes x 2 seeds = 8 samples (Duffing/LV/Brusselator/BK)
    Schema: operator_means + eigenvalues + null direction + AI-Poincare

  Problem 3: four_force_caricature
    Source: four_force_probe_results_seed11.json + seed22.json
    4 classes x 2 seeds = 8 samples (EM/weak/strong/gravity)
    Same schema as simulator_4domain

After ingest, each problem has a trained classifier + history entry.
Speaking posture: cardiac N=1/class is degenerate today (seeds future
cohort); simulator and four-force at N=2/class should show high
accuracy because Sessions 6-7 showed cross-seed cos +0.985+ within
domain.

This is "put all prior datasets in there" per Greg's directive.
NO new data runs - only reading existing JSONs.
"""
import json
import numpy as np
import pandas as pd

from od_learning_store import init_problem, ingest, fit_and_evaluate, show_history
from od_features import temporal_features


# ----------------- Problem 1: cardiac disease family --------------------------

def build_cardiac_payload(rec):
    """Convert ekg_disease_canary record to flat feature dict."""
    f = rec["fits"]
    payload = {
        "subj": f"{rec['db']}_{rec['rec_name']}",
        "label": rec["family"],
        "n_rr": int(rec["n_windows"] * 15),  # rough estimate
        "n_windows": int(rec["n_windows"]),
        "n_sup": int(rec["n_supertimes"]),
        "HR_mean": rec["HR"]["mean"], "HR_std": rec["HR"]["std"],
        "ABP_mean": rec["HRV"]["mean"],  # alias HRV->ABP for cerebro schema parity
        "ABP_std": rec["HRV"]["std"],
        "MI_mean": rec["MI"]["mean"], "MI_std": rec["MI"]["std"],
        "Ha_mean": rec["Ha"]["mean"], "Ha_std": rec["Ha"]["std"],
        "Hb_mean": rec["Hb"]["mean"], "Hb_std": rec["Hb"]["std"],
        "poly1_a": f["poly1"]["coef"]["a"], "poly1_b": f["poly1"]["coef"]["b"],
        "poly1_c": f["poly1"]["coef"]["c"],
        "r2_poly1": f["poly1"]["r2"],
        "poly2_a": f["poly2"]["coef"]["a"], "poly2_b": f["poly2"]["coef"]["b"],
        "poly2_d": f["poly2"]["coef"]["d"], "poly2_e": f["poly2"]["coef"]["e"],
        "poly2_f": f["poly2"]["coef"]["f"], "poly2_c": f["poly2"]["coef"]["c"],
        "r2_poly2": f["poly2"]["r2"],
        "exp_A": f["exp_INFO025"].get("coef", {}).get("A", 0.0),
        "exp_B": f["exp_INFO025"].get("coef", {}).get("B", 0.0),
        "r2_exp": f["exp_INFO025"].get("r2", 0.0),
        "phys_alpha": f["physics_INFO025"]["coef"]["alpha"],
        "phys_c": f["physics_INFO025"]["coef"]["c"],
        "r2_phys": f["physics_INFO025"]["r2"],
    }
    # a_hat unit vector
    v = np.array([payload["poly1_a"], payload["poly1_b"]])
    n = np.linalg.norm(v)
    payload["a_hat"] = float(v[0] / n) if n > 1e-9 else 0.0
    payload["b_hat"] = float(v[1] / n) if n > 1e-9 else 0.0
    # Temporal features from the saved arrays
    payload.update(temporal_features(rec["Ha_arr"], "Ha"))
    payload.update(temporal_features(rec["Hb_arr"], "Hb"))
    payload.update(temporal_features(rec["MI_arr"], "MI"))
    return payload


def ingest_cardiac_disease_family():
    print("\n" + "=" * 78)
    print("PROBLEM 1: cardiac_disease_family")
    print("=" * 78)
    rows = []
    for src in ["ekg_disease_canary.json", "ekg_disease_canary_expand.json"]:
        for rec in json.load(open(src)):
            if rec.get("ok"):
                rows.append(build_cardiac_payload(rec))
    df = pd.DataFrame(rows)
    print(f"Built {len(df)} payloads. Classes:")
    print(df["label"].value_counts().to_string())

    # Feature columns (everything except subj + label)
    feature_cols = [c for c in df.columns if c not in ("subj", "label")]
    init_problem("cardiac_disease_family",
                 feature_cols=feature_cols, label_col="label",
                 subject_id_col="subj", test_frac=0.20,
                 random_state=0, classifier="rf")

    n_added = ingest("cardiac_disease_family", df)
    entry = fit_and_evaluate("cardiac_disease_family",
                             note=f"prior_data_bulk_ingest_{n_added}records")
    print(f"\n  held-out: acc={entry.get('accuracy')}  bal={entry.get('balanced_accuracy')}")
    print(f"  CV-on-train: acc={entry.get('cv_accuracy')}  bal={entry.get('cv_balanced_accuracy')}")
    print(f"  NOTE: N=1/class - classification is degenerate; seeds future "
          f"cohort ingest")
    show_history("cardiac_disease_family")


# ----------------- Problem 2: simulator 4-domain ------------------------------

def build_simulator_payload(r, seed, source):
    """Convert per_domain_kbk result entry to flat feature dict."""
    payload = {
        "subj": f"{source}_{seed}_{r['name']}",
        "label": r["name"],
        "n_points": int(r["n_points"]),
        "kbk_rank_gap": float(r["kbk_rank_gap"]),
        "cos_to_artifact": float(r["cos_to_artifact_+1+1+2"]),
        "ai_two_nn_z": float(r["ai_poincare_two_nn_z"]),
        "ai_levina_bickel": float(r["ai_poincare_levina_bickel"]),
        "ai_local_pca": float(r["ai_poincare_local_pca"]),
    }
    # operator_means (6) and operator_stds (6)
    for i, v in enumerate(r["operator_means"]):
        payload[f"op_mean_{i}"] = float(v)
    for i, v in enumerate(r["operator_stds"]):
        payload[f"op_std_{i}"] = float(v)
    # eigenvalues of operator covariance (top 6)
    eigs = r.get("eigenvalues_op_cov", [])
    for i in range(6):
        payload[f"eig_{i}"] = float(eigs[i]) if i < len(eigs) else 0.0
    # null direction in 6D (V_1)
    vn6 = r.get("extract_v1_v_null_6d", [])
    for i in range(6):
        payload[f"vnull6_{i}"] = float(vn6[i]) if i < len(vn6) else 0.0
    # null direction in [2,3,4] subspace (V_1)
    vn3 = r.get("extract_v1_v_null_234", [])
    for i in range(3):
        payload[f"vnull3_{i}"] = float(vn3[i]) if i < len(vn3) else 0.0
    return payload


def ingest_simulator_4domain():
    print("\n" + "=" * 78)
    print("PROBLEM 2: simulator_4domain")
    print("=" * 78)
    rows = []
    for seed_src in [
        ("per_domain_kbk_results_seed11.json", 11),
        ("per_domain_kbk_results_seed22.json", 22),
    ]:
        src, seed = seed_src
        d = json.load(open(src))
        for r in d["results"]:
            rows.append(build_simulator_payload(r, seed, "kbk"))
    df = pd.DataFrame(rows)
    print(f"Built {len(df)} payloads. Classes:")
    print(df["label"].value_counts().to_string())

    feature_cols = [c for c in df.columns if c not in ("subj", "label")]
    init_problem("simulator_4domain",
                 feature_cols=feature_cols, label_col="label",
                 subject_id_col="subj", test_frac=0.25,
                 random_state=0, classifier="rf")

    n_added = ingest("simulator_4domain", df)
    entry = fit_and_evaluate("simulator_4domain",
                             note=f"prior_data_bulk_ingest_{n_added}records")
    print(f"\n  held-out: acc={entry.get('accuracy')}  bal={entry.get('balanced_accuracy')}")
    print(f"  CV-on-train: acc={entry.get('cv_accuracy')}  bal={entry.get('cv_balanced_accuracy')}")
    show_history("simulator_4domain")


# ----------------- Problem 3: four-force caricature ---------------------------

def ingest_four_force_caricature():
    print("\n" + "=" * 78)
    print("PROBLEM 3: four_force_caricature")
    print("=" * 78)
    rows = []
    for src, seed in [
        ("four_force_probe_results_seed11.json", 11),
        ("four_force_probe_results_seed22.json", 22),
    ]:
        d = json.load(open(src))
        for r in d["results"]:
            rows.append(build_simulator_payload(r, seed, "ffp"))
    df = pd.DataFrame(rows)
    print(f"Built {len(df)} payloads. Classes:")
    print(df["label"].value_counts().to_string())

    feature_cols = [c for c in df.columns if c not in ("subj", "label")]
    init_problem("four_force_caricature",
                 feature_cols=feature_cols, label_col="label",
                 subject_id_col="subj", test_frac=0.25,
                 random_state=0, classifier="rf")

    n_added = ingest("four_force_caricature", df)
    entry = fit_and_evaluate("four_force_caricature",
                             note=f"prior_data_bulk_ingest_{n_added}records")
    print(f"\n  held-out: acc={entry.get('accuracy')}  bal={entry.get('balanced_accuracy')}")
    print(f"  CV-on-train: acc={entry.get('cv_accuracy')}  bal={entry.get('cv_balanced_accuracy')}")
    show_history("four_force_caricature")


if __name__ == "__main__":
    ingest_cardiac_disease_family()
    ingest_simulator_4domain()
    ingest_four_force_caricature()
    print("\n\nAll prior datasets ingested. Store now has 4 problems:")
    print("  cerebro_disease_binary  (44 subjects, real cardiovascular)")
    print("  cardiac_disease_family  (10 subjects, real cardiac, N=1/class)")
    print("  simulator_4domain       (8 samples, simulated, 4 domains x 2 seeds)")
    print("  four_force_caricature   (8 samples, simulated, 4 forces x 2 seeds)")
