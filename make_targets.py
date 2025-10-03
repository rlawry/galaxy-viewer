#!/usr/bin/env python3
# Build 30 SDSS GALAXY targets across 0<=z<1 with IAU names
from astroquery.sdss import SDSS
from astropy.coordinates import SkyCoord
import astropy.units as u
import pathlib

# 10 bins of width 0.1, 3 per bin = 30
BINS = [(i/10, (i+1)/10, 3) for i in range(10)]

def iau_sdss_name(ra_deg: float, dec_deg: float) -> str:
    c = SkyCoord(ra_deg*u.deg, dec_deg*u.deg, frame="icrs")
    ra_s  = c.ra.to_string(unit=u.hour, sep="", precision=2, pad=True)                       # hhmmss.ss
    dec_s = c.dec.to_string(unit=u.deg,  sep="", precision=1, pad=True, alwayssign=True)     # ±ddmmss.s
    return f"SDSS J{ra_s}{dec_s}"

def query_bin(zmin, zmax, k):
    sql = f"""
        SELECT TOP {k}
            ra, dec, z
        FROM SpecObj
        WHERE class = 'GALAXY' AND z >= {zmin} AND z < {zmax}
        ORDER BY z
    """
    return SDSS.query_sql(sql)

def widen_and_fill(zmin, zmax, k):
    sql = f"""
        SELECT TOP {k}
            ra, dec, z
        FROM SpecObj
        WHERE class = 'GALAXY' AND z >= {max(0.0, zmin-0.02)} AND z < {min(1.0, zmax+0.02)}
        ORDER BY z
    """
    return SDSS.query_sql(sql)

def top_up(need, have_set):
    if need <= 0: return []
    t = SDSS.query_sql(f"""
        SELECT TOP {need*2}
            ra, dec, z
        FROM SpecObj
        WHERE class = 'GALAXY' AND z >= 0.0 AND z < 1.0
        ORDER BY z
    """)
    rows = []
    if t is not None:
        for row in t:
            ra = float(row["ra"]); dec = float(row["dec"])
            key = (round(ra,6), round(dec,6))
            if key in have_set: continue
            rows.append((iau_sdss_name(ra,dec), ra, dec, float(row["z"])))
            if len(rows) >= need: break
    return rows

def main():
    rows = []
    for zmin, zmax, k in BINS:
        t = query_bin(zmin, zmax, k)
        if t is None or len(t) < k:
            t2 = widen_and_fill(zmin, zmax, k)
        else:
            t2 = t
        if t2 is None or len(t2) == 0: continue
        for row in t2[:k]:
            ra = float(row["ra"]); dec = float(row["dec"]); z = float(row["z"])
            rows.append((iau_sdss_name(ra,dec), ra, dec, z))

    # de-dup and top up to 30
    uniq = {}
    for name, ra, dec, z in rows:
        uniq[(round(ra,6), round(dec,6))] = (name, ra, dec, z)
    rows = list(uniq.values())
    rows += top_up(30 - len(rows), set((round(r[1],6), round(r[2],6)) for r in rows))
    rows = rows[:30]

    out = pathlib.Path("targets_autogen.txt")
    with out.open("w", encoding="utf-8") as f:
        f.write("# name, ra_deg, dec_deg, z  (SDSS SpecObj; GALAXY class; 0<=z<1)\n")
        for name, ra, dec, z in rows:
            f.write(f"{name}, {ra:.6f}, {dec:.6f}, {z:.5f}\n")
    print(f"Wrote {out.resolve()} with {len(rows)} targets.")

if __name__ == "__main__":
    main()
