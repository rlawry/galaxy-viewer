# pip install -U astroquery astropy numpy requests
from astroquery.mast import Observations
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np, json, os, requests

# -------- settings --------
POS = SkyCoord("00h42m44.3 +41d16m09s", unit=(u.hourangle, u.deg), frame="icrs")
OUT_DIR = "m31_stis_optical_json"
JSON_PATH = os.path.join(OUT_DIR, "m31_stis.json")
IMG_NAME = "m31_sdss_dr19.jpg"
IMG_WEB_PATH = f"{OUT_DIR}/{IMG_NAME}"  # path to use inside the JSON (relative to your HTML)
RADIUS = 30*u.arcsec                    # wider to catch G750L
REB_A = 3.0                             # median-bin width in Angstrom
GR_OK = {"G430L","G750L","G430M","G750M"}  # CCD optical gratings
VIS_NM_MIN, VIS_NM_MAX = 350.0, 1030.0

os.makedirs(OUT_DIR, exist_ok=True)

def is_x1d(row):
    fn = str(row.get("productFilename") or "").lower()
    psgd = str(row.get("productSubGroupDescription") or "").upper()
    ptype = str(row.get("productType") or "").upper()
    return ptype=="SCIENCE" and (fn.endswith(("x1d.fits","sx1.fits")) or psgd in ("X1D","SX1"))

def rebin_median(xA, y, stepA):
    """Median rebin onto ~constant step (Å). xA must be sorted asc."""
    xA = np.asarray(xA, float); y = np.asarray(y, float)
    if xA.size==0: return xA, y
    edges = np.arange(xA.min(), xA.max()+stepA, stepA)
    idx = np.digitize(xA, edges)
    xb, yb = [], []
    for k in range(1, len(edges)):
        m = (idx==k)
        if np.any(m):
            xb.append(np.median(xA[m]))
            yb.append(np.median(y[m]))
    return np.array(xb), np.array(yb)

def process_file(path):
    """Return list of segment dicts from one X1D/SX1 file; already filtered and rebinned."""
    out = []
    with fits.open(path) as hdul:
        h0 = hdul[0].header
        det = (h0.get("DETECTOR") or "").upper()
        gr  = (h0.get("OPT_ELEM")  or "").upper()
        if det != "CCD" or gr not in GR_OK:
            return out

        sci = hdul["SCI"]
        tab = sci.data
        rows = range(len(tab))
        wA = np.concatenate([np.asarray(tab["WAVELENGTH"][i], float) for i in rows])
        f  = np.concatenate([np.asarray(tab["FLUX"][i],       float) for i in rows])
        dq = np.concatenate([np.asarray(tab["DQ"][i],           int) for i in rows]) if "DQ" in tab.names else None

        # mask: finite, nonzero, within band
        nm = wA/10.0
        m = np.isfinite(wA) & np.isfinite(f) & (f!=0) & (nm>=VIS_NM_MIN) & (nm<=VIS_NM_MAX)

        # apply serious DQ mask if available
        sdq = int(sci.header.get("SDQFLAGS", 0))
        if dq is not None and sdq:
            m &= ((dq & sdq) == 0)

        if not np.any(m):
            return out

        wA, f = wA[m], f[m]
        o = np.argsort(wA); wA, f = wA[o], f[o]

        # median rebin for visual smoothness
        wAb, fb = rebin_median(wA, f, REB_A)

        seg = {
            "id": (h0.get("ROOTNAME") or h0.get("OBSID") or os.path.basename(path)).strip(),
            "name": f"STIS {gr} segment",
            "instrument": f"STIS/{det} {gr}",
            "grating": gr,
            "cenwave_A": h0.get("CENWAVE"),
            "file": os.path.basename(path),
            "wavelength_A": wAb.tolist(),
            "flux": fb.tolist(),
            "flux_sdss": (fb*1e17).tolist(),
            "flux_unit": "erg s^-1 cm^-2 Å^-1",
            "flux_sdss_unit": "10^-17 erg s^-1 cm^-2 Å^-1",
            "ra_deg": float(h0.get("RA_TARG")  or h0.get("RA_APER")  or POS.ra.deg),
            "dec_deg": float(h0.get("DEC_TARG") or h0.get("DEC_APER") or POS.dec.deg),
        }
        out.append(seg)
    return out

# ---- query + download ----
obs = Observations.query_criteria(
    coordinates=POS, radius=RADIUS,
    obs_collection="HST", instrument_name="STIS",
    calib_level=[2,3]
)
if len(obs) == 0:
    raise SystemExit("No STIS results.")

products = Observations.get_product_list(obs)
keep = products[[is_x1d(r) for r in products]]
if len(keep) == 0:
    raise SystemExit("No X1D/SX1 products.")

manifest = Observations.download_products(keep, download_dir=OUT_DIR, mrp_only=False)

# ---- extract segments ----
segments = []
for row in manifest:
    path = row.get("Local Path")
    try:
        segments.extend(process_file(path))
    except Exception:
        continue

if not segments:
    raise SystemExit("Downloaded STIS files, but none in optical CCD gratings after filtering.")

# ---- stitched spectrum across all segments ----
all_wA = np.concatenate([np.asarray(s["wavelength_A"], float) for s in segments])
all_f  = np.concatenate([np.asarray(s["flux"],         float) for s in segments])
o = np.argsort(all_wA); all_wA, all_f = all_wA[o], all_f[o]
wA_st, f_st = rebin_median(all_wA, all_f, REB_A)

# position for image
ra = float(segments[0]["ra_deg"])
dec = float(segments[0]["dec_deg"])

stitched = {
    "id": "M31_STIS_STITCHED",
    "name": "M31 nucleus — HST/STIS stitched",
    "instrument": "HST/STIS CCD (G430L+G750L/M)",
    "wavelength_A": wA_st.tolist(),
    "flux": f_st.tolist(),
    "flux_sdss": (f_st*1e17).tolist(),
    "flux_unit": "erg s^-1 cm^-2 Å^-1",
    "flux_sdss_unit": "10^-17 erg s^-1 cm^-2 Å^-1",
    "ra_deg": ra, "dec_deg": dec,
    "image": {"path": IMG_WEB_PATH},
    "credits": {"spectrum":"HST/STIS via MAST","image":"SDSS DR19"}
}

# ---- fetch SDSS cutout (optional) ----
sdss_url = ("https://skyserver.sdss.org/dr19/SkyServerWS/ImgCutout/getjpeg"
            f"?ra={ra:.7f}&dec={dec:.7f}&scale=0.396&width=512&height=512")
img_path = os.path.join(OUT_DIR, IMG_NAME)
try:
    r = requests.get(sdss_url, timeout=60)
    if r.ok:
        with open(img_path, "wb") as fh:
            fh.write(r.content)
except Exception:
    pass  # image is optional

# attach image path to each segment too
for s in segments:
    s["image"] = {"path": IMG_WEB_PATH}

# ---- write ARRAY JSON compatible with the viewer ----
payload = [stitched] + segments
with open(JSON_PATH, "w", encoding="utf-8") as f:
    json.dump(payload, f, separators=(",", ":"), ensure_ascii=False)

print(f"wrote {JSON_PATH} with {len(payload)} objects "
      f"(stitched + {len(segments)} segments).")
