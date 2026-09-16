
from pathlib import Path
import shutil, os

home = Path.home()
downloads = home / "Downloads"
here = Path.cwd()

def newest_file(name):
    hits = [p for p in downloads.rglob(name) if p.is_file()]
    if not hits:
        return None
    return max(hits, key=lambda p: p.stat().st_mtime)

def newest_dir(name):
    hits = [p for p in downloads.rglob(name) if p.is_dir()]
    if not hits:
        return None
    return max(hits, key=lambda p: p.stat().st_mtime)

for name in ["vina.exe", "LZX-2-73.pdbqt"]:
    src = newest_file(name)
    if src is None:
        raise FileNotFoundError(f"Could not find {name} anywhere under Downloads")
    dst = here / name
    shutil.copy2(src, dst)
    print("Copied", src, "->", dst)

for dname in ["slahs_pdbqt", "lvtkl_pdbqt"]:
    src = newest_dir(dname)
    if src is None:
        raise FileNotFoundError(f"Could not find folder {dname} anywhere under Downloads")
    dst = here / dname
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    print("Copied folder", src, "->", dst)

print("\nInput collection complete.")
