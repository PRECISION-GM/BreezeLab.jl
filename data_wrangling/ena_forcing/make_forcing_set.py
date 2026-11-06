"""Write one ENA forcing set (SAM lsf, sfc, snd, prm) for the Covert-protocol LES, 24–48 h from 18 July 2017 06 UTC.

Two modes:

  --mode covert_mean   the plain sonde-nudging set: the Covert 06–18 UTC *time-mean* lsf (trapezoidal mean of
                       the five 3-hourly records, on the Covert pressure levels) and time-mean prescribed surface
                       fluxes, both held constant (two identical records at day0 and --end_day); SST = --sst.
  --mode parametric    the parameterized forcing model (height-coordinate lsf every --record_minutes):

      subsidence   wls(z, t) = −D(t) z below z_D, blended smoothly (cos² over [z_D, z_D + taper]) into the Covert
                   time-mean wls above;  D(t) = D0 [1 + a_D sin(2π (t − t0)/P_D + φ_D)]
      BL advection tls(z, t) = tls_cov(z) + m(z) [ΔT + A_T sin(2π (t − t0)/P_T + φ_T)]
                   qls(z, t) = qls_cov(z) + m(z) [Δq + A_q sin(2π (t − t0)/P_q + φ_q)]
                   with m = 1 below z_BL and the same cos² taper to 0 above (tls_cov, qls_cov: Covert time means)
      surface      H = flux_scale_H × H̄_cov, LE = flux_scale_LE × LĒ_cov (prescribed), τ and SST as covert_mean
      winds        uls/vls (= the geostrophic wind, the Covert lsf has no ug/vg columns) are the Covert means

  The baseline (D0 = −w̄_cov(z_ref)/z_ref, a_D = 0, ΔT = Δq = A = 0, flux scales 1) differs from covert_mean only in
  the boundary-layer *shape* of wls (linear in z instead of the Covert profile).

Times t are hours since day0 = 199.25 (18 July 06 UTC); t0 = 0 by default; phases in radians.
Units on the command line: D in s⁻¹, ΔT and A_T in K/day, Δq and A_q in g/kg/day, heights in m, periods in h.
The free-troposphere sonde nudging is a model option (--nudging_tq ...), not part of these files; the snd written
here is the multi-record sonde series (--snd), identical for every set.

    python3 make_forcing_set.py <out dir> --mode covert_mean --covert <sonde0531 dir> --snd <series snd>
    python3 make_forcing_set.py <out dir> --mode parametric --covert ... --snd ... --D0_scale 1.5 --a_D 0.5 --P_D 6
"""
import argparse, json, os, shutil
import numpy as np

R, CP, G = 287.0, 1004.0, 9.81
DAY0 = 199.25


def read_records(path):
    lines = [l for l in open(path).read().splitlines() if l.strip()]
    header, i, recs = lines[0], 1, []
    while i < len(lines):
        t = lines[i].split()
        day, nlev, pres0 = float(t[0]), int(t[1]), float(t[2])
        rows = np.array([[float(x.replace("D", "e").replace("d", "e")) for x in l.split()] for l in lines[i + 1:i + 1 + nlev]])
        recs.append((day, pres0, rows))
        i += nlev + 1
    return header, recs


def trapezoid_mean(days, arrays):
    days = np.asarray(days)
    w = np.zeros(len(days))
    dt = np.diff(days)
    w[:-1] += dt / 2
    w[1:] += dt / 2
    return sum(wi * a for wi, a in zip(w, arrays)) / w.sum()


def sam_heights(p, θ, p0):
    T = θ * (p / 1000.0) ** (R / CP)
    z = np.zeros_like(p)
    z[0] = R / G * T[0] * np.log(p0 / p[0])
    for i in range(1, len(p)):
        z[i] = z[i - 1] + R / (2 * G) * (T[i] + T[i - 1]) * np.log(p[i - 1] / p[i])
    return z


def smooth_step(z, z0, depth):
    """0 below z0, 1 above z0 + depth, sin² in between."""
    ξ = np.clip((np.asarray(z, float) - z0) / depth, 0, 1)
    return np.sin(0.5 * np.pi * ξ) ** 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--mode", choices=["covert_mean", "parametric"], required=True)
    ap.add_argument("--covert", required=True, help="directory with the Covert lsf, sfc, prm and the sonde0531 snd")
    ap.add_argument("--snd", required=True, help="the multi-record sonde series snd")
    ap.add_argument("--end_day", type=float, default=201.5)
    ap.add_argument("--record_minutes", type=float, default=15.0)
    ap.add_argument("--sst", type=float, default=294.9)
    # subsidence
    ap.add_argument("--D0", type=float, default=None, help="s⁻¹; default −w̄_cov(z_ref)/z_ref")
    ap.add_argument("--D0_scale", type=float, default=1.0)
    ap.add_argument("--z_ref", type=float, default=1000.0)
    ap.add_argument("--z_D", type=float, default=1500.0)
    ap.add_argument("--a_D", type=float, default=0.0)
    ap.add_argument("--P_D", type=float, default=24.0)
    ap.add_argument("--phi_D", type=float, default=0.0)
    # BL horizontal advection
    ap.add_argument("--z_BL", type=float, default=1500.0)
    ap.add_argument("--taper", type=float, default=500.0)
    ap.add_argument("--dT", type=float, default=0.0)
    ap.add_argument("--A_T", type=float, default=0.0)
    ap.add_argument("--P_T", type=float, default=24.0)
    ap.add_argument("--phi_T", type=float, default=0.0)
    ap.add_argument("--dq", type=float, default=0.0)
    ap.add_argument("--A_q", type=float, default=0.0)
    ap.add_argument("--P_q", type=float, default=24.0)
    ap.add_argument("--phi_q", type=float, default=0.0)
    ap.add_argument("--t0", type=float, default=0.0)
    # surface
    ap.add_argument("--flux_scale", type=float, default=1.0)
    ap.add_argument("--flux_scale_H", type=float, default=None)
    ap.add_argument("--flux_scale_LE", type=float, default=None)
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)
    shutil.copy(os.path.join(a.covert, "prm"), os.path.join(a.out, "prm"))
    shutil.copy(a.snd, os.path.join(a.out, "snd"))

    # Covert time means (06–18 UTC, trapezoidal)
    lsf_header, lsf = read_records(os.path.join(a.covert, "lsf"))
    days = [r[0] for r in lsf]
    p_lsf = lsf[0][2][:, 1]
    assert all(np.allclose(r[2][:, 1], p_lsf) for r in lsf), "Covert lsf records must share pressure levels"
    mean = trapezoid_mean(days, [r[2] for r in lsf])            # columns z p tls qls uls vls wls
    pres0 = lsf[0][1]
    sfc_rows = np.loadtxt(os.path.join(a.covert, "sfc"), skiprows=1)
    sfc_mean = trapezoid_mean(sfc_rows[:, 0], list(sfc_rows[:, 1:]))   # sst H LE tau

    # Heights of the Covert pressure levels (setdata.f90 with the record-1 sounding θ, as BreezeLab's record_heights)
    _, snd0 = read_records(os.path.join(a.covert, "snd"))
    s_day, s_p0, s_rows = snd0[0]
    lp = np.log(s_rows[:, 1])
    θ_lsf = np.interp(np.log(p_lsf)[::-1], lp[::-1], s_rows[::-1, 2])[::-1]
    z_lsf = sam_heights(p_lsf, θ_lsf, pres0)

    def cov(col, z):                                            # Covert time mean, linear in z
        return np.interp(z, z_lsf, mean[:, col])

    w_ref = cov(6, a.z_ref)
    D0_default = -w_ref / a.z_ref
    D0 = (a.D0 if a.D0 is not None else D0_default) * a.D0_scale
    sH = a.flux_scale if a.flux_scale_H is None else a.flux_scale_H
    sLE = a.flux_scale if a.flux_scale_LE is None else a.flux_scale_LE
    H, LE, τ = sfc_mean[1] * sH, sfc_mean[2] * sLE, sfc_mean[3]

    # sfc: constant, two rows
    with open(os.path.join(a.out, "sfc"), "w") as fh:
        fh.write("day sst(K) H(W/m2) LE(W/m2) TAU(m2/s2)\n")
        for d in (DAY0, a.end_day):
            fh.write(f"{d:.5f}      {a.sst:.3f}      {H:.4f}      {LE:.4f}   {τ:.4f}\n")

    fmt = lambda x: f"{x: .5E}"
    with open(os.path.join(a.out, "lsf"), "w") as fh:
        fh.write(lsf_header + "\n")
        if a.mode == "covert_mean":
            for d in (DAY0, a.end_day):
                fh.write(f"{d:.5f}          {len(p_lsf)}      {pres0:.2f}      day,levels,pres0\n")
                for r in mean:
                    fh.write(f"-9999.000 {r[1]:6.1f} {fmt(r[2])} {fmt(r[3])} {r[4]:8.3f} {r[5]:8.3f} {r[6]: .6f}\n")
        else:
            z_top_fine = a.z_D + a.taper + 1000
            z_levels = np.r_[np.arange(0.0, z_top_fine + 1, 25.0), z_lsf[z_lsf > z_top_fine + 50]]
            p_levels = np.exp(np.interp(z_levels, z_lsf, np.log(p_lsf)))
            m_BL = 1 - smooth_step(z_levels, a.z_BL, a.taper)
            blend_D = smooth_step(z_levels, a.z_D, a.taper)
            tls0, qls0 = cov(2, z_levels), cov(3, z_levels)
            uls, vls, wcov = cov(4, z_levels), cov(5, z_levels), cov(6, z_levels)
            nt = int(round((a.end_day - DAY0) * 24 * 60 / a.record_minutes)) + 1
            for n in range(nt):
                t = n * a.record_minutes / 60                     # hours since day0
                D = D0 * (1 + a.a_D * np.sin(2 * np.pi * (t - a.t0) / a.P_D + a.phi_D))
                wls = (1 - blend_D) * (-D * z_levels) + blend_D * wcov
                tls = tls0 + m_BL * (a.dT + a.A_T * np.sin(2 * np.pi * (t - a.t0) / a.P_T + a.phi_T)) / 86400
                qls = qls0 + m_BL * 1e-3 * (a.dq + a.A_q * np.sin(2 * np.pi * (t - a.t0) / a.P_q + a.phi_q)) / 86400
                fh.write(f"{DAY0 + t / 24:.6f}          {len(z_levels)}      {pres0:.2f}      day,levels,pres0\n")
                for k in range(len(z_levels)):
                    fh.write(f"{z_levels[k]:9.2f} {p_levels[k]:8.2f} {fmt(tls[k])} {fmt(qls[k])} "
                             f"{uls[k]:8.3f} {vls[k]:8.3f} {wls[k]: .7f}\n")

    params = dict(vars(a))
    params.update(D0_effective=D0, D0_default=D0_default, w_cov_at_z_ref=w_ref, H=H, LE=LE, tau=τ,
                  covert_mean_H=sfc_mean[1], covert_mean_LE=sfc_mean[2], covert_mean_sst=sfc_mean[0])
    with open(os.path.join(a.out, "forcing_parameters.json"), "w") as fh:
        json.dump(params, fh, indent=1)
    print(json.dumps({k: params[k] for k in ("mode", "D0_effective", "D0_default", "w_cov_at_z_ref", "H", "LE", "tau")}))


if __name__ == "__main__":
    main()
