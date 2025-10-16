# pip install -U astroquery astropy numpy requests
from astroquery.mast import Observations
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np, json, os, requests, math

pos = SkyCoord("00h42m44.3 +41d16m09s", unit=(u.hourangle, u.deg), frame="icrs")
out_dir = "m31_stis_optical_json"
os.makedirs(out_dir, exist_ok=True)

# Query HST/STIS near nucleus
obs = Observations.query_criteria(
    coordinates=pos, radius=20*u.arcsec,
    obs_collection="HST", instrument_name="STIS",
    calib_level=[2,3]
)
if len(obs)==0: raise SystemExit("No STIS results.")

# Products → X1D/SX1 only
products = Observations.get_product_list(obs)
def is_x1d(r):
    fn = str(r.get("productFilename") or "").lower()
    psgd = str(r.get("productSubGroupDescription") or "").upper()
    return (fn.endswith(("x1d.fits","sx1.fits")) or psgd in ("X1D","SX1")) \
           and str(r["productType"]).upper()=="SCIENCE"

keep = products[[is_x1d(r) for r in products]]
if len(keep)==0: raise SystemExit("No X1D/SX1.")

# Download all candidate spectra
manifest = Observations.download_products(keep, download_dir=out_dir, mrp_only=False)

# Collect optical CCD long-slit segments
segments = []
for p in manifest:
    path = p["Local Path"]
    try:
        with fits.open(path) as hdul:
            h0 = hdul[0].header
            if h0.get("DETECTOR","").upper()!="CCD":  # skip MAMA
                continue
            gr = (h0.get("OPT_ELEM") or "").upper()
            if gr not in {"G430L","G750L","G430M","G750M"}:
                continue
            sci = hdul["SCI"].data
            w = np.asarray(sci["WAVELENGTH"][0], float)  # Å
            f = np.asarray(sci["FLUX"][0], float)
            # Keep optical window 350–950 nm
            nm = w/10.0
            m = (nm>=350.0) & (nm<=950.0) & np.isfinite(f) & np.isfinite(nm)
            if not np.any(m): 
                continue
            seg = {
                "grating": gr,
                "cenwave_A": h0.get("CENWAVE"),
                "obsid": str(h0.get("OBSID") or h0.get("ROOTNAME")),
                "file": os.path.basename(path),
                "wavelength_nm": nm[m].tolist(),
                "flux": f[m].tolist(),
                "flux_unit": "erg s^-1 cm^-2 Å^-1"
            }
            segments.append(seg)
            ra = float(h0.get("RA_TARG") or h0.get("RA_APER") or pos.ra.deg)
            dec = float(h0.get("DEC_TARG") or h0.get("DEC_APER") or pos.dec.deg)
    except Exception:
        continue

if not segments:
    raise SystemExit("Downloaded STIS files, but none in optical CCD gratings.")

# Optional stitch to one array (simple concat+sort)
STITCH = True
if STITCH:
    all_nm = np.concatenate([np.array(s["wavelength_nm"]) for s in segments])
    all_fx = np.concatenate([np.array(s["flux"]) for s in segments])
    o = np.argsort(all_nm)
    stitched = {"wavelength_nm": all_nm[o].tolist(),
                "flux": all_fx[o].tolist(),
                "flux_unit": "erg s^-1 cm^-2 Å^-1"}
else:
    stitched = None

# SDSS cutout
sdss_url = ("https://skyserver.sdss.org/dr19/SkyServerWS/ImgCutout/getjpeg"
            f"?ra={ra:.7f}&dec={dec:.7f}&scale=0.396&width=512&height=512")
img_path = os.path.join(out_dir, "m31_sdss_dr19.jpg")
try:
    r = requests.get(sdss_url, timeout=60); 
    open(img_path,"wb").write(r.content) if r.ok else None
except Exception:
    img_path = None

# JSON payload: keep segments + stitched
record = {
    "name": "M31 nucleus",
    "ra_deg": ra, "dec_deg": dec,
    "instrument": "STIS/CCD",
    "segments": segments,           # per-grating chunks
    "stitched": stitched,           # wide band if STITCH=True
    "sdss_image": os.path.basename(img_path) if img_path else None,
    "credits": {"spectrum":"HST/STIS via MAST","image":"SDSS DR19"}
}

with open(os.path.join(out_dir, "m31_stis.json"), "w", encoding="utf-8") as f:
    json.dump(record, f, separators=(",",":"), ensure_ascii=False)

print("segments:", len(segments), "stitched:", bool(stitched))
