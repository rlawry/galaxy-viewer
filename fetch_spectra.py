#!/usr/bin/env python3
# fetch_spectra.py — SDSS via query_region, coords or SDSS IAU names
import sys, json, pathlib, re, traceback
from typing import List, Dict, Any, Optional
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.sdss import SDSS

MAX_RADIUS_ARCMIN = 2.9  # stay under 3' cap on some astroquery builds

# --- Input parsing ------------------------------------------------------------

def load_targets(path: str) -> List[Dict[str, Any]]:
    """targets.txt lines: name, ra, dec[, z]
       name may be an SDSS IAU name: 'SDSS Jhhmmss.ss+ddmmss.s'
       ra/dec may be sexagesimal (HH:MM:SS, ±DD:MM:SS) or decimal deg."""
    rows = []
    for ln in pathlib.Path(path).read_text(encoding="utf-8").splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        parts = [p.strip() for p in ln.split(",")]
        name = parts[0]
        ra   = parts[1] if len(parts) > 1 and parts[1] else None
        dec  = parts[2] if len(parts) > 2 and parts[2] else None
        z    = float(parts[3]) if len(parts) > 3 and parts[3] else None
        if (ra is None or dec is None):
            ra_dec = parse_sdss_j_name(name)  # try SDSS IAU name
            if ra_dec:
                ra, dec = ra_dec
        rows.append({"name": name, "ra": ra, "dec": dec, "z": z})
    return rows

def parse_sdss_j_name(name: str) -> Optional[tuple]:
    """
    Parse 'SDSS Jhhmmss.ss+ddmmss.s' or 'Jhhmmss.ss+ddmmss.s' to (ra_deg, dec_deg).
    """
    m = re.search(r'(?:SDSS\s+)?J(?P<rah>\d{2})(?P<ram>\d{2})(?P<ras>\d{2}(?:\.\d+)?)(?P<sgn>[+-])(?P<dd>\d{2})(?P<dm>\d{2})(?P<ds>\d{2}(?:\.\d+)?)', name, re.IGNORECASE)
    if not m:
        return None
    rah, ram, ras = float(m["rah"]), float(m["ram"]), float(m["ras"])
    dd, dm, ds = float(m["dd"]), float(m["dm"]), float(m["ds"])
    sgn = -1.0 if m["sgn"] == "-" else 1.0
    ra_hours = rah + ram/60.0 + ras/3600.0
    ra_deg = 15.0 * ra_hours
    dec_deg = sgn * (dd + dm/60.0 + ds/3600.0)
    return (f"{ra_deg:.8f}", f"{dec_deg:.8f}")  # return as strings to pass into parse_coord

def parse_coord(ra_str: str, dec_str: str) -> SkyCoord:
    sex = (":" in ra_str) or (" " in ra_str) or ("h" in ra_str.lower())
    if sex:
        return SkyCoord(ra_str, dec_str, unit=(u.hourangle, u.deg), frame="icrs")
    return SkyCoord(float(ra_str)*u.deg, float(dec_str)*u.deg, frame="icrs")

# --- SDSS wrappers ------------------------------------------------------------

def _query_region(coord: SkyCoord, radius):
    # Support old/new astroquery signatures
    try:
        return SDSS.query_region(coord, spectro=True, radius=radius)
    except TypeError:
        return SDSS.query_region(coordinates=coord, spectro=True, radius=radius)

def _get_spectra(matches):
    try:
        return SDSS.get_spectra(matches=matches)
    except TypeError:
        return SDSS.get_spectra(matches)

def _row_coords(table):
    # Try lowercase or uppercase RA/DEC
    if "ra" in table.colnames and "dec" in table.colnames:
        return np.array(table["ra"], dtype=float), np.array(table["dec"], dtype=float)
    if "RA" in table.colnames and "DEC" in table.colnames:
        return np.array(table["RA"], dtype=float), np.array(table["DEC"], dtype=float)
    return None, None

def fetch_sdss_spectrum(coord: SkyCoord) -> Optional[Dict[str, Any]]:
    # radii: 1′, 2′, 2.9′
    for rad_arcmin in (1.0, 2.0, MAX_RADIUS_ARCMIN):
        try:
            xid = _query_region(coord, radius=rad_arcmin*u.arcmin)
        except ValueError:
            continue
        if xid is None or len(xid) == 0:
            continue

        # choose nearest row if possible
        idx = 0
        ra_col, dec_col = _row_coords(xid)
        if ra_col is not None:
            cand = SkyCoord(ra_col*u.deg, dec_col*u.deg)
            idx = np.argmin(coord.separation(cand).arcsec)

        one = xid[idx:idx+1]
        try:
            sp_list = _get_spectra(one)
        except Exception:
            sp_list = None
        if not sp_list:
            continue

        hdul = sp_list[0]
        try:
            data = hdul[1].data
            wave = np.power(10.0, data["loglam"])  # Å
            flux = data["flux"]
            z_hdr = None
            for k in ("Z","Z_NOQSO","Z_VI","ZPIPE"):
                if k in hdul[0].header:
                    try:
                        z_hdr = float(hdul[0].header[k]); break
                    except Exception:
                        pass

            meta = {}
            for k in ("PLATE","MJD","FIBERID","SPECOBJID"):
                if k in hdul[0].header:
                    meta[k.lower()] = hdul[0].header[k]

            return {
                "wavelength_A": wave.astype(float).tolist(),
                "flux": flux.astype(float).tolist(),
                "z_hdr": z_hdr,
                "sdss_meta": meta
            }
        finally:
            hdul.close()
    return None

# --- Main ---------------------------------------------------------------------

def main():
    if len(sys.argv) < 2:
        print("Usage: python fetch_spectra.py targets.txt"); sys.exit(1)
    rows = load_targets(sys.argv[1])
    out = []
    for t in rows:
        try:
            if not t["ra"] or not t["dec"]:
                print(f"[skip] no coordinates for: {t['name']}")
                continue
            c = parse_coord(t["ra"], t["dec"])
            spec = fetch_sdss_spectrum(c)
            if spec is None:
                print(f"[skip] no SDSS spectrum near: {t['name']}")
                continue
            out.append({
                "name": t["name"],                     # preserve provided/IAU name
                "ra_deg": c.ra.deg,
                "dec_deg": c.dec.deg,
                "z": t["z"] if t["z"] is not None else spec.get("z_hdr", None),
                "wavelength_A": spec["wavelength_A"],
                "flux": spec["flux"],
                "sdss": spec.get("sdss_meta", {})
            })
            print(f"[ok] {t['name']}  points={len(spec['wavelength_A'])}")
        except Exception as e:
            print(f"[err] {t['name']}: {e}")
            traceback.print_exc(limit=1)
    if not out:
        print("No spectra fetched. Nothing to write."); return
    pathlib.Path("spectra.json").write_text(json.dumps(out), encoding="utf-8")
    print("Wrote spectra.json")

if __name__ == "__main__":
    main()
