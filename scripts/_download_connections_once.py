"""One-shot resumable download of proofread_connections_783.feather (no curl)."""

from __future__ import annotations

import sys
from pathlib import Path

import requests
from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "data" / "raw" / "proofread_connections_783.feather"
URL = (
    "https://zenodo.org/records/10676866/files/"
    "proofread_connections_783.feather?download=1"
)
EXPECTED = 852_022_274
STATUS = ROOT / "data" / "raw" / "_dl_status.txt"


def log(msg: str) -> None:
    STATUS.write_text(msg + "\n", encoding="utf-8")
    print(msg, flush=True)


def main() -> int:
    DEST.parent.mkdir(parents=True, exist_ok=True)
    attempt = 0
    while True:
        attempt += 1
        have = DEST.stat().st_size if DEST.exists() else 0
        if have >= EXPECTED:
            log(f"DONE size={have}")
            return 0
        if attempt > 40:
            log(f"FAIL after {attempt} attempts size={have}")
            return 1

        headers = {"Range": f"bytes={have}-"} if have else {}
        log(f"attempt={attempt} resume_from={have}")
        try:
            with requests.get(
                URL, stream=True, timeout=180, allow_redirects=True, headers=headers
            ) as r:
                if have and r.status_code == 200:
                    # Server ignored Range; restart
                    have = 0
                    mode = "wb"
                else:
                    mode = "ab" if have else "wb"
                r.raise_for_status()
                total = int(r.headers.get("content-length") or 0)
                bar_total = (have + total) if (have and r.status_code == 206) else (
                    total or EXPECTED
                )
                with DEST.open(mode) as f, tqdm(
                    total=bar_total,
                    initial=have if mode == "ab" else 0,
                    unit="B",
                    unit_scale=True,
                    desc="connections",
                ) as bar:
                    for chunk in r.iter_content(chunk_size=1 << 20):
                        if not chunk:
                            continue
                        f.write(chunk)
                        bar.update(len(chunk))
        except Exception as exc:  # noqa: BLE001
            log(f"error attempt={attempt}: {exc!r}")
            continue


if __name__ == "__main__":
    raise SystemExit(main())
