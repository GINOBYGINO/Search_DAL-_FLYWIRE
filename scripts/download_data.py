"""Download FlyWire v783 raw tables into data/raw/ and refresh data/README.md.

Codex 程式下載需 api_token；本專案改以公開來源取得同等欄位：
- Zenodo：proofread_connections
- flywire_annotations v2.1.0：含 top_nt、classification、cell_type 的神經元註解表
下載後再從註解表拆出 neurons／classification／cell_types 三份便於使用的表到 data/processed/。
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import subprocess
import sys
import zipfile
from datetime import date
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
DATA_README = ROOT / "data" / "README.md"

MATERIALIZATION = "783"
ANNOTATION_TAG = "v2.1.0"

RAW_FILES = [
    {
        "filename": "proofread_connections_783.feather",
        "url": (
            "https://zenodo.org/records/10676866/files/"
            "proofread_connections_783.feather?download=1"
        ),
        "source": "Zenodo record 10676866 (FlyWire Whole-brain Connectome Connectivity Data)",
        "version": "783.0",
        "purpose": "校對後神經元對神經元連線（neuropil、突觸數、傳導物質平均機率）",
        "expected_bytes": 852_022_274,
    },
    {
        "filename": "Supplemental_file1_neuron_annotations.tsv",
        "url": (
            "https://raw.githubusercontent.com/flyconnectome/flywire_annotations/"
            f"{ANNOTATION_TAG}/supplemental_files/Supplemental_file1_neuron_annotations.tsv"
        ),
        "source": f"flyconnectome/flywire_annotations {ANNOTATION_TAG}",
        "version": ANNOTATION_TAG,
        "purpose": (
            "神經元註解母表（含 top_nt、flow／super_class／cell_class、cell_type 等；"
            "等同 Codex neurons + classification + cell_types 的公開對應）"
        ),
        "expected_bytes": 27_015_208,
    },
]


def md5_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        while True:
            block = f.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def download_with_curl(url: str, dest: Path) -> None:
    """Prefer curl for large files: native resume, fewer Windows file locks."""
    cmd = [
        "curl.exe",
        "-L",
        "--retry",
        "20",
        "--retry-delay",
        "2",
        "--retry-all-errors",
        "-C",
        "-",
        "-o",
        str(dest),
        url,
    ]
    print(" ", " ".join(cmd))
    subprocess.run(cmd, check=False)


def download(
    url: str,
    dest: Path,
    force: bool = False,
    expected_bytes: int | None = None,
) -> bool:
    """Download url to dest. Returns True if a new download happened.

    Large / resumable downloads use curl.exe when available. Re-downloads if
    existing file is smaller than expected_bytes.
    """
    if (
        dest.exists()
        and not force
        and (expected_bytes is None or dest.stat().st_size >= expected_bytes)
    ):
        print(f"skip (exists): {dest.name}")
        return False

    if dest.exists() and expected_bytes and dest.stat().st_size < expected_bytes:
        print(
            f"incomplete {dest.name}: {dest.stat().st_size:,} < {expected_bytes:,}; resume"
        )

    dest.parent.mkdir(parents=True, exist_ok=True)
    if force and dest.exists():
        dest.unlink()

    # Large feather: curl resume is more reliable on Windows than requests+rename
    use_curl = expected_bytes is not None and expected_bytes >= 50_000_000
    if use_curl and shutil.which("curl.exe"):
        for attempt in range(1, 16):
            print(f"curl attempt {attempt} -> {dest.name}")
            download_with_curl(url, dest)
            size = dest.stat().st_size if dest.exists() else 0
            if expected_bytes is None or size >= expected_bytes:
                print(f"saved: {dest.name} ({size:,} bytes)")
                return True
            print(f"still incomplete ({size:,}); retrying…")
        raise RuntimeError(
            f"download incomplete for {dest.name}: "
            f"{dest.stat().st_size if dest.exists() else 0:,} < {expected_bytes:,} bytes"
        )

    tmp = dest.with_suffix(dest.suffix + ".part")
    resume_from = tmp.stat().st_size if tmp.exists() and not force else 0
    headers = {"Range": f"bytes={resume_from}-"} if resume_from else {}
    mode = "ab" if resume_from else "wb"

    with requests.get(
        url, stream=True, timeout=120, allow_redirects=True, headers=headers
    ) as r:
        if resume_from and r.status_code == 200:
            resume_from = 0
            mode = "wb"
        r.raise_for_status()
        total = int(r.headers.get("content-length") or 0)
        if resume_from and r.status_code == 206 and total:
            total = resume_from + total
        with tmp.open(mode) as f, tqdm(
            total=total or expected_bytes,
            initial=resume_from if mode == "ab" else 0,
            unit="B",
            unit_scale=True,
            desc=dest.name,
        ) as bar:
            for chunk in r.iter_content(chunk_size=1 << 20):
                if chunk:
                    f.write(chunk)
                    bar.update(len(chunk))

    with tmp.open("rb") as f:
        magic = f.read(4)
    if magic[:2] == b"PK":
        extract_dir = dest.parent / f"_extract_{dest.stem}"
        if extract_dir.exists():
            shutil.rmtree(extract_dir)
        extract_dir.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(tmp, "r") as zf:
            zf.extractall(extract_dir)
        candidates = list(extract_dir.rglob("*.csv")) + list(extract_dir.rglob("*.tsv"))
        if not candidates:
            raise RuntimeError(f"zip from {url} contained no csv/tsv")
        chosen = next(
            (p for p in candidates if "neuron" in p.name.lower()),
            candidates[0],
        )
        if dest.exists():
            dest.unlink()
        chosen.replace(dest)
        shutil.rmtree(extract_dir)
        tmp.unlink(missing_ok=True)
        print(f"extracted -> {dest.name}")
        return True

    if expected_bytes is not None and tmp.stat().st_size < expected_bytes:
        raise RuntimeError(
            f"download incomplete for {dest.name}: "
            f"{tmp.stat().st_size:,} < {expected_bytes:,} bytes"
        )

    if dest.exists():
        dest.unlink()
    tmp.replace(dest)
    print(f"saved: {dest.name} ({dest.stat().st_size:,} bytes)")
    return True


def export_tables_from_annotations(annot_path: Path) -> dict[str, Path]:
    """Split annotation dump into neurons / classification / cell_types tables."""
    PROCESSED.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(annot_path, sep="\t", low_memory=False)

    neuron_cols = [
        c
        for c in [
            "root_id",
            "supervoxel_id",
            "pos_x",
            "pos_y",
            "pos_z",
            "soma_x",
            "soma_y",
            "soma_z",
            "nucleus_id",
            "top_nt",
            "top_nt_conf",
            "known_nt",
            "known_nt_source",
            "side",
            "status",
        ]
        if c in df.columns
    ]
    class_cols = [
        c
        for c in [
            "root_id",
            "flow",
            "super_class",
            "cell_class",
            "cell_sub_class",
            "ito_lee_hemilineage",
            "hartenstein_hemilineage",
            "morphology_group",
            "side",
            "nerve",
            "status",
        ]
        if c in df.columns
    ]
    type_cols = [
        c
        for c in [
            "root_id",
            "cell_type",
            "hemibrain_type",
            "side",
            "vfb_id",
            "fbbt_id",
        ]
        if c in df.columns
    ]

    out = {
        "neurons_783.csv": PROCESSED / "neurons_783.csv",
        "classification_783.csv": PROCESSED / "classification_783.csv",
        "cell_types_783.csv": PROCESSED / "cell_types_783.csv",
    }
    df[neuron_cols].to_csv(out["neurons_783.csv"], index=False)
    df[class_cols].to_csv(out["classification_783.csv"], index=False)
    df[type_cols].to_csv(out["cell_types_783.csv"], index=False)
    for name, path in out.items():
        print(f"wrote processed/{name} ({path.stat().st_size:,} bytes)")
    return out


def file_row(path: Path, version: str, download_date: str, source: str, url: str, purpose: str) -> str:
    exists = path.exists()
    size = f"{path.stat().st_size:,} bytes" if exists else "尚未下載"
    checksum = md5_file(path) if exists else "—"
    when = download_date if exists else "—"
    return (
        f"| `{path.relative_to(ROOT / 'data').as_posix()}` ({size}; md5 `{checksum}`) | "
        f"{version} | {when} | {source}<br>`{url}` | {purpose} |"
    )


def write_data_readme(download_date: str) -> None:
    lines = [
        "# 資料出處紀錄",
        "",
        f"專案 materialization：**v{MATERIALIZATION}**",
        f"註解鎖定：**flywire_annotations {ANNOTATION_TAG}**（對齊 Nature 2024 公開版）",
        f"本 README 最近更新下載日：**{download_date}**",
        "",
        "原始檔放在 `raw/`，**不要手改**。由註解母表拆出的表在 `processed/`。報告引用請依下表。",
        "",
        "## 為什麼不用 Codex API 直接抓 neurons？",
        "",
        "Codex `download_resource` 現在需要帳號 `api_token`。"
        f"公開的 `flywire_annotations {ANNOTATION_TAG}` 註解表已含 `top_nt`／`top_nt_conf`"
        "以及 classification、cell_type 欄位，與論文公開版一致，故以此為準。",
        "",
        "| 本機路徑 | 版本 | 下載日期 | 來源／網址 | 用途 |",
        "|---|---|---|---|---|",
    ]

    for meta in RAW_FILES:
        path = RAW / meta["filename"]
        lines.append(
            file_row(
                path,
                meta["version"],
                download_date,
                meta["source"],
                meta["url"],
                meta["purpose"],
            )
        )

    annot_url = RAW_FILES[1]["url"]
    annot_source = RAW_FILES[1]["source"]
    derived = [
        (
            PROCESSED / "neurons_783.csv",
            "由註解母表匯出（含 top_nt）",
        ),
        (
            PROCESSED / "classification_783.csv",
            "由註解母表匯出（flow／super_class／cell_class 等）",
        ),
        (
            PROCESSED / "cell_types_783.csv",
            "由註解母表匯出（cell_type／hemibrain_type 等）",
        ),
    ]
    for path, purpose in derived:
        lines.append(
            file_row(
                path,
                ANNOTATION_TAG,
                download_date,
                annot_source + "（本機拆表）",
                annot_url,
                purpose,
            )
        )

    lines.extend(
        [
            "",
            "## 重跑下載",
            "",
            "```bash",
            "conda activate flywire-dal",
            "python scripts/download_data.py",
            "# 強制重抓：",
            "python scripts/download_data.py --force",
            "```",
            "",
            "## 未下載",
            "",
            "- 完整突觸點表 `flywire_synapses_783.feather`（Zenodo，約 9.5 GB）：第 3 步若需要再補。",
            "",
        ]
    )
    DATA_README.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {DATA_README}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--force",
        action="store_true",
        help="重新下載即使檔案已存在",
    )
    args = parser.parse_args()

    RAW.mkdir(parents=True, exist_ok=True)
    today = date.today().isoformat()

    for meta in RAW_FILES:
        dest = RAW / meta["filename"]
        try:
            download(
                meta["url"],
                dest,
                force=args.force,
                expected_bytes=meta.get("expected_bytes"),
            )
        except Exception as exc:  # noqa: BLE001
            print(f"ERROR downloading {meta['filename']}: {exc}", file=sys.stderr)
            return 1

    annot = RAW / "Supplemental_file1_neuron_annotations.tsv"
    try:
        export_tables_from_annotations(annot)
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR exporting processed tables: {exc}", file=sys.stderr)
        return 1

    write_data_readme(today)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
