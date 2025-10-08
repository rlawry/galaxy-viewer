#!/usr/bin/env python3
# make_targets.py — 40 SDSS GALAXY targets with 0<=z<0.05 (simple + robust)
from astroquery.sdss import SDSS
from astropy.coordinates import SkyCoord
import astropy.units as u
import pathlib, time, random, re

# 10 bins of width 0.005, 4 per bin = 40 total
BINS = [(i*0.005, (i+1)*0.005, 4) for i in range(10)]  # 0.000–0.050

def iau_sdss_name(ra_deg: float, dec_deg: float) -> str:
    c = SkyCoord(ra_deg*u.deg, dec_deg*u.deg, frame="icrs")
    ra_s  = c.ra.to_string(unit=u.hour, sep="", precision=2, pad=True)
    dec_s = c.dec.to_string(unit=u.deg,  sep="", precision=1, pad=True, alwayssign=True)
    return f"SDSS J{ra_s}{dec_s}"

def _query_sql(sql: str, retries=2):
    m = re.search(r"SELECT\s+TOP\s+(\d+)", sql, re.I)
    tops = [int(m.group(1))] if m else [4]
    tops += [max(1, tops[0]//2)]
    for top in tops:
        sql2 = re.sub(r"(SELECT\s+TOP\s+)\d+", rf"\g<1>{top}", sql, flags=re.I)
        for dr in (17, 16):
            for i in range(retries):
                try:
                    try:
                        return SDSS.query_sql(sql2, data_release=dr)
                    except TypeError:
                        return SDSS.query_sql(sql2)
                except Exception:
                    time.sleep(1.0 + 0.5*i + random.uniform(0,0.3))
    return None

def query_bin(zmin, zmax, k):
    sql = f"""
        SELECT TOP {k}
            ra  AS ra,
            dec AS dec,
            z   AS z
        FROM SpecObj
        WHERE class = 'GALAXY'
          AND z >= {zmin} AND z < {zmax}
        ORDER BY z ASC
    """
    return _query_sql(sql)

def main():
    rows = []
    for zmin, zmax, k in BINS:
        t = query_bin(zmin, zmax, k)
        if not t: continue
        for row in t:
            try:
                ra = float(row["ra"]); dec = float(row["dec"]); z = float(row["z"])
            except Exception:
                continue
            rows.append((iau_sdss_name(ra,dec), ra, dec, z))

    # de-dup and cap at 40
    uniq = {}
    for name, ra, dec, z in rows:
        uniq[(round(ra,6), round(dec,6))] = (name, ra, dec, z)
    rows = list(uniq.values())[:40]

    out = pathlib.Path("targets2_autogen.txt")
    with out.open("w", encoding="utf-8") as f:
        f.write("# label(IAU), ra_deg, dec_deg, z  — SDSS GALAXY 0<=z<0.05\n")
        for name, ra, dec, z in rows:
            f.write(f"{name}, {ra:.6f}, {dec:.6f}, {z:.5f}\n")
    print(f"Wrote {out.resolve()} with {len(rows)} targets.")

if __name__ == "__main__":
    main()
