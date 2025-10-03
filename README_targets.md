# Build a 30-target list (z=0.5..5.0)

Run locally:

```bash
pip install --upgrade astroquery astropy
python make_targets.py
python fetch_spectra.py targets_autogen.txt
```

Notes:
- Uses SDSS SpecObj QSOs for high-z coverage.
- Falls back by widening bins if any are empty.
