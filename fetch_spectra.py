#!/usr/bin/env python3
# fetch_spectra.py — SDSS via query_region, supports: label,ra,dec[,z][,id]
import sys, json, pathlib, re, traceback
from typing import List, Dict, Any, Optional, Tuple
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.sdss import SDSS
import requests
from PIL import Image
from io import BytesIO

MAX_RADIUS_ARCMIN = 2.9  # stay under 3'

# ----------------------- cutout helpers -----------------------

def _safe_slug(s: str) -> str:
    import re
    return re.sub(r"[^A-Za-z0-9._-]+", "_", (s or "target")).strip("_")

CUTOUT_SIZE_PIX = 512
CUTOUT_PXSCALE = 0.262        # "/pix (Legacy Surveys native scale)
LEGACY_LAYER = "ls-dr10-grz"  # color grz composite

def fetch_color_cutout_jpg(coord: SkyCoord, label: str) -> Optional[Tuple[str, str]]:
    """Return (path, source) or None on failure."""
    ra, dec = coord.ra.deg, coord.dec.deg
    url_legacy = (f"https://www.legacysurvey.org/viewer/cutout.jpg?"
                  f"ra={ra:.8f}&dec={dec:.8f}&layer={LEGACY_LAYER}"
                  f"&pixscale={CUTOUT_PXSCALE:.3f}&size={CUTOUT_SIZE_PIX}")
    try:
        r = requests.get(url_legacy, timeout=30); r.raise_for_status()
        img = Image.open(BytesIO(r.content)).convert("RGB")
        outdir = pathlib.Path("cutouts2"); outdir.mkdir(exist_ok=True)
        path = outdir / f"{_safe_slug(label)}_legacy_{CUTOUT_SIZE_PIX}px.jpg"
        img.save(path, quality=90)
        return (str(path), "Legacy DR10 grz")
    except Exception as e:
        print(f"[warn] Legacy cutout failed: {e}; trying SDSS")
        try:
            url_sdss = (f"https://skyserver.sdss.org/dr19/SkyServerWS/ImgCutout/getjpeg?"
                        f"ra={ra:.8f}&dec={dec:.8f}&scale=0.4&width={CUTOUT_SIZE_PIX}&height={CUTOUT_SIZE_PIX}")
            r2 = requests.get(url_sdss, timeout=30); r2.raise_for_status()
            outdir = pathlib.Path("cutouts2"); outdir.mkdir(exist_ok=True)
            path = outdir / f"{_safe_slug(label)}_sdss_{CUTOUT_SIZE_PIX}px.jpg"
            path.write_bytes(r2.content)
            return (str(path), "SDSS DR19 JPEG")
        except Exception as e2:
            print(f"[warn] SDSS cutout failed: {e2}")
            return None

# ----------------------- input parsing -----------------------

SDSS_J_RE = re.compile(
    r'(?:SDSS\s+)?J(?P<rah>\d{2})(?P<ram>\d{2})(?P<ras>\d{2}(?:\.\d+)?)(?P<sgn>[+-])(?P<dd>\d{2})(?P<dm>\d{2})(?P<ds>\d{2}(?:\.\d+)?)',
    re.IGNORECASE
)

def parse_sdss_j_name(s: Optional[str]) -> Optional[tuple[str, str]]:
    if not s: return None
    m = SDSS_J_RE.search(s)
    if not m: return None
    rah, ram, ras = float(m["rah"]), float(m["ram"]), float(m["ras"])
    dd, dm, ds = float(m["dd"]), float(m["dm"]), float(m["ds"])
    sgn = -1.0 if m["sgn"] == "-" else 1.0
    ra_deg = 15.0 * (rah + ram/60.0 + ras/3600.0)
    dec_deg = sgn * (dd + dm/60.0 + ds/3600.0)
    return (f"{ra_deg:.8f}", f"{dec_deg:.8f}")

def load_targets(path: str) -> List[Dict[str, Any]]:
    """
    CSV rows: label, ra, dec[, z][, id]
    - label: UI label (e.g., 'NGC 5548 nucleus')
    - id: precise identifier (e.g., 'SDSS J141759.53+250812.4')
    - ra/dec sexagesimal or degrees; if missing, derived from id/label if SDSS J...
    """
    rows = []
    for ln in pathlib.Path(path).read_text(encoding="utf-8").splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#"): continue
        parts = [p.strip() for p in ln.split(",")]
        label = parts[0] if len(parts) > 0 else None
        ra    = parts[1] if len(parts) > 1 and parts[1] else None
        dec   = parts[2] if len(parts) > 2 and parts[2] else None
        z     = float(parts[3]) if len(parts) > 3 and parts[3] else None
        qid   = parts[4] if len(parts) > 4 and parts[4] else None

        if (ra is None or dec is None):
            ra_dec = parse_sdss_j_name(qid) or parse_sdss_j_name(label)
            if ra_dec:
                ra, dec = ra_dec

        rows.append({"label": label, "ra": ra, "dec": dec, "z": z, "id": qid})
    return rows

def parse_coord(ra_str: str, dec_str: str) -> SkyCoord:
    sex = (":" in ra_str) or (" " in ra_str) or ("h" in ra_str.lower())
    if sex:
        return SkyCoord(ra_str, dec_str, unit=(u.hourangle, u.deg), frame="icrs")
    return SkyCoord(float(ra_str)*u.deg, float(dec_str)*u.deg, frame="icrs")

# ----------------------- SDSS helpers -----------------------

def _query_region(coord: SkyCoord, radius):
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
    if "ra" in table.colnames and "dec" in table.colnames:
        return np.array(table["ra"], dtype=float), np.array(table["dec"], dtype=float)
    if "RA" in table.colnames and "DEC" in table.colnames:
        return np.array(table["RA"], dtype=float), np.array(table["DEC"], dtype=float)
    return None, None

def fetch_sdss_spectrum(coord: SkyCoord) -> Optional[Dict[str, Any]]:
    # try growing radii; for each, test up to N nearest spectra until one has non-zero flux
    for rad_arcmin in (1.0, 2.0, MAX_RADIUS_ARCMIN):
        try:
            xid = _query_region(coord, radius=rad_arcmin*u.arcmin)
        except ValueError:
            continue
        if xid is None or len(xid) == 0:
            continue

        # sort candidates by distance and keep first N
        ra_col, dec_col = _row_coords(xid)
        order = np.arange(len(xid))
        if ra_col is not None:
            cand = SkyCoord(ra_col*u.deg, dec_col*u.deg)
            order = np.argsort(coord.separation(cand).arcsec)
        NTRY = min(8, len(order))
        sel = xid[np.array(order[:NTRY])]

        # pull spectra for the selected candidates
        try:
            sp_list = _get_spectra(sel)
        except Exception:
            sp_list = None
        if not sp_list:
            continue

        for hdul in sp_list:
            try:
                data = hdul[1].data
                wave = np.power(10.0, data["loglam"])  # Å
                flux = np.array(data["flux"], dtype=float)
                if not _valid_flux(flux):
                    continue  # try next candidate

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
                    "sdss_meta": meta,
                }
            finally:
                hdul.close()
    return None


def _valid_flux(a) -> bool:
    a = np.asarray(a, dtype=float)
    if a.size == 0: return False
    if not np.isfinite(a).any(): return False
    if np.nanmax(np.abs(a)) < 1e-12: return False
    if np.all(np.isclose(a, 0.0, atol=1e-6, rtol=0.0)): return False
    if np.nanstd(a) < 1e-12: return False
    return True

# ----------------------- main -----------------------

def main():
    if len(sys.argv) < 2:
        print("Usage: python fetch_spectra.py targets.txt"); sys.exit(1)
    rows = load_targets(sys.argv[1])
    out = []
    for t in rows:
        label = t.get("label") or t.get("id") or "Object"
        try:
            if not t["ra"] or not t["dec"]:
                print(f"[skip] no coordinates for: {label}")
                continue
            c = parse_coord(t["ra"], t["dec"])
            spec = fetch_sdss_spectrum(c)
            if spec is None:
                print(f"[skip] no SDSS spectrum near: {label}")
                continue

            img = None
            try:
                img = fetch_color_cutout_jpg(c, label)  # (path, source) or None
            except Exception as e:
                print(f"[warn] cutout failed for {label}: {e}")

            out.append({
                "name": label,
                "id": t.get("id"),
                "ra_deg": c.ra.deg,
                "dec_deg": c.dec.deg,
                "z": t["z"] if t["z"] is not None else spec.get("z_hdr", None),
                "wavelength_A": spec["wavelength_A"],
                "flux": spec["flux"],
                "sdss": spec.get("sdss_meta", {}),
                "image": ({"path": img[0], "source": img[1], "size_px": CUTOUT_SIZE_PIX} if img else None),
            })
            print(f"[ok] {label}  points={len(spec['wavelength_A'])}")
        except Exception as e:
            print(f"[err] {label}: {e}")
            traceback.print_exc(limit=1)
    if not out:
        print("No spectra fetched. Nothing to write."); return
    pathlib.Path("spectra.json").write_text(json.dumps(out), encoding="utf-8")
    print("Wrote spectra.json")

if __name__ == "__main__":
    main()
