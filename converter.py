#!/usr/bin/env python3
# Reformat STIS/GMOS JSON -> Galaxy Spectra Viewer format
# Usage: python reformat_for_viewer.py input.json spectra.json

import json, sys, math, os

def as_list(x):
    return [] if x is None else list(x)

def finite_pair(nm, fx):
    return (isinstance(nm, (int, float)) and isinstance(fx, (int, float))
            and math.isfinite(nm) and math.isfinite(fx))

def pick_wl_flux(rec):
    """Return wavelength in Å and flux arrays."""
    # 1) stitched block preferred
    st = rec.get("stitched") or {}
    nm = st.get("wavelength_nm"); fx = st.get("flux")
    if isinstance(nm, list) and isinstance(fx, list) and len(nm) == len(fx):
        wlA = [v * 10.0 for v in nm]
        return wlA, fx

    # 2) root wavelength_nm/flux
    nm = rec.get("wavelength_nm"); fx = rec.get("flux")
    if isinstance(nm, list) and isinstance(fx, list) and len(nm) == len(fx):
        wlA = [v * 10.0 for v in nm]
        return wlA, fx

    # 3) root wavelength_A/flux already in Å
    wa = rec.get("wavelength_A"); fx = rec.get("flux")
    if isinstance(wa, list) and isinstance(fx, list) and len(wa) == len(fx):
        return wa, fx

    # 4) merge all segments (if present)
    segs = rec.get("segments") or []
    nm_all, fx_all = [], []
    for s in segs:
        nm = s.get("wavelength_nm")
        fx = s.get("flux")
        if isinstance(nm, list) and isinstance(fx, list) and len(nm) == len(fx):
            nm_all.extend(nm); fx_all.extend(fx)
    if nm_all and fx_all and len(nm_all) == len(fx_all):
        wlA = [v * 10.0 for v in nm_all]
        # sort by wavelength
        idx = sorted(range(len(wlA)), key=lambda i: wlA[i])
        wlA = [wlA[i] for i in idx]; fx_all = [fx_all[i] for i in idx]
        return wlA, fx_all

    # nothing usable
    return [], []

def to_viewer_obj(rec, idx=0):
    wlA, fx = pick_wl_flux(rec)

    # keep only finite pairs
    pairs = [(w, f) for w, f in zip(wlA, fx) if finite_pair(w, f)]
    wlA = [w for w, _ in pairs]
    fx  = [f for _, f in pairs]

    name = (rec.get("name") or rec.get("targname") or rec.get("id") or "Object").strip() if isinstance(rec.get("name") or rec.get("targname") or rec.get("id"), str) else (rec.get("id") or f"Object {idx+1}")
    img = None
    if isinstance(rec.get("image"), dict) and isinstance(rec["image"].get("path"), str):
        img = {"path": rec["image"]["path"]}
    elif isinstance(rec.get("sdss_image"), str):
        img = {"path": rec["sdss_image"]}

    out = {
        "id": rec.get("id") or name,
        "name": name,
        "wavelength_A": wlA,
        "flux": fx,
        "image": img,
        "z": float(rec.get("z") or 0.0),
        "ra_deg": rec.get("ra_deg"),
        "dec_deg": rec.get("dec_deg"),
    }
    return out

def main():
    inp = sys.argv[1] if len(sys.argv) > 1 else "m31_stis_optical_json/m31_stis.json"
    outp = sys.argv[2] if len(sys.argv) > 2 else "spectra2.json"

    with open(inp, "r", encoding="utf-8") as f:
        src = json.load(f)

    # Accept single object or array of objects
    records = src if isinstance(src, list) else [src]
    out = [to_viewer_obj(r, i) for i, r in enumerate(records)]

    # optional: drop empties with no spectrum
    out = [o for o in out if o.get("wavelength_A") and o.get("flux")]

    with open(outp, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, separators=(",", ":"))

    print(f"Wrote {outp} with {len(out)} object(s).")

if __name__ == "__main__":
    main()
