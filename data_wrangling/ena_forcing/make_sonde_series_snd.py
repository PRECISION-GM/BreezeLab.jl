"""Multi-record SAM `snd` from the ARM ENA radiosondes, 18 July 05:31 → 20 July 05:32 UTC 2017.

Record 1 (day 199.25, 18 July 06 UTC) is copied verbatim from the sonde0531 initial sounding
(`make_sonde0531_snd.py`), so the initial condition of the sonde-start runs does not change.
Every later record is one ARM sonde (ascending part only), binned in *height* above mean sea
level (the model surface): 25 m bins below 3 km, 50 m to 10 km, 100 m above, written as a
height-coordinate SAM record (z given; the p column carries the bin-mean pressure). Above the
sonde top (or in empty bins at the top) the public Covert profile of record 1 is appended, its
pressure levels converted to heights with SAM's `setdata.f90` hydrostatic recipe.

    θ = T (1000/p)^0.2857,  q = 622 e / (p − e) [g/kg] with e from rh over water (Bolton),
    as in make_sonde0531_snd.py.

Record times are the launch times as SAM day of year (18 July 06 UTC = 199.25). A final copy
of the last sonde is appended at `--hold_day` (default 202.0) so the series explicitly holds
the last observed profile past 20 July 05:32 (the model would hold it anyway: profile time
series clamp to the end records).

    python3 make_sonde_series_snd.py <sonde0531 dir> <sonde dir> <out snd> [--first 20170718.1138 --last 20170720.0532]
"""
import argparse, glob, os, re
import numpy as np
import scipy.io

R, CP, G = 287.0, 1004.0, 9.81


def read_snd(path):
    lines = open(path).read().splitlines()
    header, i, blocks = lines[0], 1, []
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        t = lines[i].split()
        nlev = int(t[1])
        blocks.append((lines[i], lines[i + 1:i + 1 + nlev]))
        i += nlev + 1
    return header, blocks


def sam_heights(p, θ, p0):
    """setdata.f90: heights of pressure levels p [hPa] with potential temperature θ, surface pressure p0."""
    T = θ * (p / 1000.0) ** (R / CP)
    z = np.zeros_like(p)
    z[0] = R / G * T[0] * np.log(p0 / p[0])
    for i in range(1, len(p)):
        z[i] = z[i - 1] + R / (2 * G) * (T[i] + T[i - 1]) * np.log(p[i - 1] / p[i])
    return z


def sonde_day(path):
    m = re.search(r"(\d{8})\.(\d{6})\.cdf$", path)
    ymd, hms = m.group(1), m.group(2)
    month, day = int(ymd[4:6]), int(ymd[6:8])
    assert ymd[:4] == "2017" and month == 7
    doy = 181 + day                          # 1 July 2017 = day 182 (SAM: 1.0 = 1 Jan 00 UTC)
    return doy + (int(hms[:2]) + int(hms[2:4]) / 60 + int(hms[4:]) / 3600) / 24


def binned_sonde(path):
    d = scipy.io.netcdf_file(path, "r", mmap=False)
    g = lambda v: np.array(d.variables[v].data, float)
    z, p, T, rh, u, v = g("alt"), g("pres"), g("tdry") + 273.15, g("rh"), g("u_wind"), g("v_wind")
    ok = np.isfinite(z) & np.isfinite(p) & np.isfinite(T) & np.isfinite(rh) & np.isfinite(u) & np.isfinite(v)
    ok &= (p > 0) & (T > 150) & (rh >= 0) & (np.abs(u) < 150) & (np.abs(v) < 150)
    z, p, T, rh, u, v = (a[ok] for a in (z, p, T, rh, u, v))
    up = np.r_[True, p[1:] < np.minimum.accumulate(p)[:-1]]     # ascending part: new minimum pressure only
    z, p, T, rh, u, v = (a[up] for a in (z, p, T, rh, u, v))
    θ = T * (1000 / p) ** 0.2857
    e = rh / 100 * 6.112 * np.exp(17.67 * (T - 273.15) / (T - 29.65))
    q = 622 * e / (p - e)
    edges = np.r_[np.arange(0, 3000, 25.0), np.arange(3000, 10000, 50.0), np.arange(10000, 30001, 100.0)]
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (z >= lo) & (z < hi)
        if m.sum() == 0:
            continue
        rows.append((0.5 * (lo + hi), p[m].mean(), θ[m].mean(), q[m].mean(), u[m].mean(), v[m].mean()))
    return np.array(rows), d.variables["pres"].data[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("init_dir")
    ap.add_argument("sonde_dir")
    ap.add_argument("out")
    ap.add_argument("--first", default="20170718.1138")
    ap.add_argument("--last", default="20170720.0532")
    ap.add_argument("--hold_day", type=float, default=202.0)
    a = ap.parse_args()

    header, blocks = read_snd(os.path.join(a.init_dir, "snd"))
    first_header, first_rows = blocks[0]
    assert abs(float(first_header.split()[0]) - 199.25) < 1e-6, "record 1 must be day 199.25"
    pres0 = float(first_header.split()[2])
    r1 = np.array([list(map(float, l.split())) for l in first_rows])          # z p θ q u v (z = -9999)
    z1 = sam_heights(r1[:, 1], r1[:, 2], pres0)

    files = sorted(glob.glob(os.path.join(a.sonde_dir, "enasondewnpnC1.b1.*.cdf")))
    key = lambda f: re.search(r"(\d{8}\.\d{4})\d\d\.cdf$", f).group(1)
    files = [f for f in files if a.first <= key(f) <= a.last]
    assert files, "no sondes selected"

    records = []
    for f in files:
        rows, psfc = binned_sonde(f)
        top = rows[-1, 0]
        fallback = [(z1[n], *r1[n, 1:]) for n in range(len(z1)) if z1[n] > top + 100]
        if fallback:
            rows = np.vstack([rows, np.array(fallback)])
        records.append((sonde_day(f), float(psfc), rows, os.path.basename(f), top))

    # Compare the free troposphere of the first sonde record with record 1 (diagnostic only)
    print("record 1 vs first sonde record (free troposphere):")
    for zc in (1300, 1600, 2000, 2500, 3000, 4000, 5000, 7000, 9000, 12000, 16000, 19000):
        r = records[0][2]
        print(f"  z={zc:6.0f} m  θ {np.interp(zc, z1, r1[:, 2]):7.2f} → {np.interp(zc, r[:, 0], r[:, 2]):7.2f} K   "
              f"q {np.interp(zc, z1, r1[:, 3]):6.3f} → {np.interp(zc, r[:, 0], r[:, 3]):6.3f} g/kg   "
              f"u {np.interp(zc, z1, r1[:, 4]):6.2f} → {np.interp(zc, r[:, 0], r[:, 4]):6.2f}   "
              f"v {np.interp(zc, z1, r1[:, 5]):6.2f} → {np.interp(zc, r[:, 0], r[:, 5]):6.2f}")

    last = records[-1]
    if a.hold_day > last[0]:
        records.append((a.hold_day, last[1], last[2], last[3] + " (hold)", last[4]))

    with open(a.out, "w") as fh:
        fh.write(header + "\n")
        fh.write(first_header + "\n")
        for l in first_rows:
            fh.write(l + "\n")
        for day, psfc, rows, name, top in records:
            fh.write(f"{day:13.5f} {len(rows):11d} {psfc:12.2f}      day,levels,pres0  ARM {name}\n")
            for zc, pp, th, qq, uu, vv in rows:
                fh.write(f"{zc:9.2f} {pp:12.3f} {th:12.3f} {qq:12.4f} {uu:12.5f} {vv:12.5f}\n")
    for day, psfc, rows, name, top in records:
        print(f"  day {day:9.4f}  {name:45s} levels {len(rows):4d}  sonde top {top:7.0f} m  psfc {psfc:7.1f} hPa")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
