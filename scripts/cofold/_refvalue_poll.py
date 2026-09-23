import json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "cofold"))
import importlib.util
spec = importlib.util.spec_from_file_location(
    "opc", Path(__file__).resolve().parent / "openprotein_cofold.py")
opc = importlib.util.module_from_spec(spec); spec.loader.exec_module(opc)

deadline = time.time() + float(sys.argv[1] if len(sys.argv) > 1 else 1800)
while time.time() < deadline:
    tot = {"collected": 0, "pending": 0, "failed": 0}
    for eng in ("protenix_v2", "protenix"):
        r = opc.collect(eng, "refvalue")
        for k in tot:
            tot[k] += r.get(k, 0)
        print(f"  {eng}: {r['collected']} ok / {r['pending']} pending / {r['failed']} failed",
              flush=True)
    print(f"[{time.strftime('%H:%M:%S')}] TOTAL {tot}", flush=True)
    if tot["pending"] == 0:
        print("ALL TERMINAL"); break
    time.sleep(90)
