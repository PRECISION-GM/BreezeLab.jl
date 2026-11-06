"""Inputs and manifest of the ENA parameterized-forcing ensemble (24 h from 18 July 2017 06 UTC).

One-at-a-time sweeps around the baseline of make_forcing_set.py (--mode parametric) plus combined members,
including droplet-number members (p3_n50, p3_n100) on the baseline forcing for comparing forcing and
microphysics sensitivity. Writes <root>/inputs/<member>/{lsf,sfc,snd,prm,forcing_parameters.json} and
<root>/manifest.tsv (member → parameters → input and output directories).

    python3 make_forcing_ensemble.py <root> --covert <sonde0531 dir> --snd <series snd>
"""
import argparse, json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))

# name: (generator options, microphysics)
MEMBERS = [
    ("base",            {},                                   "p3_n75"),
    ("D0x0.5",          {"D0_scale": 0.5},                    "p3_n75"),
    ("D0x1.5",          {"D0_scale": 1.5},                    "p3_n75"),
    ("aD0.5_P3h",       {"a_D": 0.5, "P_D": 3},               "p3_n75"),
    ("aD0.5_P6h",       {"a_D": 0.5, "P_D": 6},               "p3_n75"),
    ("aD0.5_P12h",      {"a_D": 0.5, "P_D": 12},              "p3_n75"),
    ("aD1_P3h",         {"a_D": 1.0, "P_D": 3},               "p3_n75"),
    ("aD1_P6h",         {"a_D": 1.0, "P_D": 6},               "p3_n75"),
    ("aD1_P12h",        {"a_D": 1.0, "P_D": 12},              "p3_n75"),
    ("tls+1",           {"dT": 1.0},                          "p3_n75"),
    ("tls-1",           {"dT": -1.0},                         "p3_n75"),
    ("qls+1",           {"dq": 1.0},                          "p3_n75"),
    ("qls-1",           {"dq": -1.0},                         "p3_n75"),
    ("flux1.2",         {"flux_scale": 1.2},                  "p3_n75"),
    ("flux0.8",         {"flux_scale": 0.8},                  "p3_n75"),
    ("Nd50",            {},                                   "p3_n50"),
    ("Nd100",           {},                                   "p3_n100"),
    ("warm_dry2",       {"dT": 2.0, "dq": -2.0},              "p3_n75"),
    ("cold_moist2",     {"dT": -2.0, "dq": 2.0},              "p3_n75"),
    ("D0x1.5_Nd50",     {"D0_scale": 1.5},                    "p3_n50"),
    ("D0x0.5_Nd100",    {"D0_scale": 0.5},                    "p3_n100"),
    ("tls_sin2_P12h",   {"A_T": 2.0, "P_T": 12},              "p3_n75"),
    ("qls_sin2_P12h",   {"A_q": 2.0, "P_q": 12},              "p3_n75"),
    ("aD1_P6h_Nd50",    {"a_D": 1.0, "P_D": 6},               "p3_n50"),
]

COLUMNS = ["D0_scale", "a_D", "P_D", "phi_D", "dT", "A_T", "P_T", "dq", "A_q", "P_q", "flux_scale"]
DEFAULTS = {"D0_scale": 1.0, "a_D": 0.0, "P_D": 24.0, "phi_D": 0.0, "dT": 0.0, "A_T": 0.0, "P_T": 24.0,
            "dq": 0.0, "A_q": 0.0, "P_q": 24.0, "flux_scale": 1.0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--covert", required=True)
    ap.add_argument("--snd", required=True)
    ap.add_argument("--end_day", type=float, default=200.5)
    a = ap.parse_args()
    os.makedirs(os.path.join(a.root, "inputs"), exist_ok=True)
    rows = []
    for n, (name, opts, micro) in enumerate(MEMBERS):
        member = f"m{n:02d}_{name}"
        inputs = os.path.join(a.root, "inputs", member)
        cmd = [sys.executable, os.path.join(HERE, "make_forcing_set.py"), inputs, "--mode", "parametric",
               "--covert", a.covert, "--snd", a.snd, "--end_day", str(a.end_day)]
        for k, v in opts.items():
            cmd += [f"--{k}", str(v)]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
        params = json.load(open(os.path.join(inputs, "forcing_parameters.json")))
        values = {**DEFAULTS, **opts}
        rows.append([member, micro, f"{params['D0_effective']:.4e}"] + [str(values[c]) for c in COLUMNS] +
                    [f"{params['H']:.3f}", f"{params['LE']:.3f}", inputs, os.path.join(a.root, member)])
    header = (["member", "microphysics", "D0_per_s"] + [
        "D0_scale", "a_D", "P_D_h", "phi_D_rad", "dT_K_per_day", "A_T_K_per_day", "P_T_h",
        "dq_g_per_kg_per_day", "A_q_g_per_kg_per_day", "P_q_h", "flux_scale"] + ["H_W_m2", "LE_W_m2", "input_dir", "output_dir"])
    with open(os.path.join(a.root, "manifest.tsv"), "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print(f"{len(rows)} members → {os.path.join(a.root, 'manifest.tsv')}")


if __name__ == "__main__":
    main()
