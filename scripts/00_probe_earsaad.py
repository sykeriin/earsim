"""Probe the Ear-SAAD BIDS zip (Zenodo 16536441) without downloading all 12 GB.

Lists the archive via HTTP range requests, pulls the small metadata files and one
subject's EEG, and prints channel names/types, sampling rate and whether
electrode positions exist. Output goes to data/earsaad/.
"""
from __future__ import annotations

import sys
from pathlib import Path

from remotezip import RemoteZip

from earsim.head import DATA_DIR

URL = "https://zenodo.org/records/16536441/files/bids_dataset.zip"
OUT = DATA_DIR / "earsaad"
META = ("participants", "dataset_description", "README", "_electrodes.tsv",
        "_coordsystem.json", "_channels.tsv", "_eeg.json")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with RemoteZip(URL) as z:
        infos = z.infolist()
        print(f"{len(infos)} entries, {sum(i.file_size for i in infos) / 1e9:.1f} GB uncompressed")
        for i in infos[:25]:
            print(f"  {i.filename}  {i.file_size / 1e6:.1f} MB")
        small = [i for i in infos if any(m in i.filename for m in META) and i.file_size < 5e6]
        # keep metadata for the first subject only plus top-level files
        subs = sorted({p for i in small for p in i.filename.split("/") if p.startswith("sub-")})
        keep = [i for i in small if not any(p.startswith("sub-") and p != (subs[0] if subs else "")
                                            for p in i.filename.split("/"))]
        for i in keep:
            z.extract(i, OUT)
        print(f"extracted {len(keep)} metadata files to {OUT}")
        eeg = [i for i in infos if i.filename.endswith(("_eeg.edf", "_eeg.set", "_eeg.bdf", "_eeg.vhdr"))
               and (subs[0] if subs else "sub-01") in i.filename]
        print("subject-1 EEG files:", [(i.filename, round(i.file_size / 1e6)) for i in eeg])
        if "--eeg" in sys.argv and eeg:
            z.extract(sorted(eeg, key=lambda i: i.file_size)[0], OUT)
    for f in sorted(OUT.rglob("*_electrodes.tsv"))[:1] + sorted(OUT.rglob("*_channels.tsv"))[:1]:
        print(f"\n== {f.relative_to(OUT)}")
        print("\n".join(f.read_text().splitlines()[:70]))


if __name__ == "__main__":
    main()
