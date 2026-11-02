"""Download Ear-SAAD EEG files (no stimuli) from the Zenodo BIDS zip via range requests.

Resumable: skips files already present at the right size. Stimuli (.wav, ~1.4 GB)
are only needed for the auditory-attention experiment (E3) and are skipped here.
"""
import sys
from pathlib import Path

from remotezip import RemoteZip

from earsim.head import DATA_DIR

URL = "https://zenodo.org/records/16536441/files/bids_dataset.zip"
OUT = DATA_DIR / "earsaad"
with_stimuli = "--stimuli" in sys.argv

with RemoteZip(URL) as z:
    infos = [i for i in z.infolist() if not i.is_dir() and (with_stimuli or "stimuli" not in i.filename)]
    for n, i in enumerate(infos, 1):
        dest = OUT / i.filename
        if dest.exists() and dest.stat().st_size == i.file_size:
            continue
        print(f"[{n}/{len(infos)}] {i.filename} ({i.file_size / 1e6:.0f} MB)", flush=True)
        z.extract(i, OUT)
print("done", flush=True)
