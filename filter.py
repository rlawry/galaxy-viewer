#!/usr/bin/env python3
# Usage: python filter_visible_band.py spectra2.json spectra.json
import json, sys, math

MIN_NM, MAX_NM = 380.0, 740.0   # keep the original visible range

def filter_obj(o):
    wlA = o.get("wavelength_A") or []
    fx  = o.get("flux") or []
    if not (isinstance(wlA, list) and isinstance(fx, list) and len(wlA) == len(fx)):
      return None
    wlA_out, fx_out = [], []
    for w, f in zip(wlA, fx):
      if isinstance(w, (int, float)) and isinstance(f, (int, float)) \
         and math.isfinite(w) and math.isfinite(f):
        nm = 0.1 * w   # Å → nm
        if MIN_NM <= nm <= MAX_NM:
          wlA_out.append(w); fx_out.append(f)
    if not wlA_out:
      return None
    o2 = dict(o)
    o2["wavelength_A"] = wlA_out
    o2["flux"] = fx_out
    return o2

def main():
    inp = sys.argv[1] if len(sys.argv) > 1 else "spectra2.json"
    out = sys.argv[2] if len(sys.argv) > 2 else "spectra2_1.json"
    with open(inp, "r", encoding="utf-8") as f:
        src = json.load(f)
    objs = src if isinstance(src, list) else [src]
    out_objs = [fo for o in objs if (fo := filter_obj(o))]
    with open(out, "w", encoding="utf-8") as f:
        json.dump(out_objs, f, ensure_ascii=False, separators=(",", ":"))
    print(f"Wrote {out} with {len(out_objs)} object(s).")

if __name__ == "__main__":
    main()
