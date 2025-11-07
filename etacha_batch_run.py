# This script automates batch processing of ETACHA calculations
# using a template ETACHA input file and a BETACHA parameter file.
import re, csv, time, subprocess
from pathlib import Path

# Executable path for ETACHA4
ETACHA_EXE   = r"C:\Program Files\LISEcute\ETACHA4.exe"

# Directory where ETACHA outputs results. Don't forget to clean up later.
RESULTS_DIR  = Path(r"My Documents\LISEcute\results")

OUTPUT_CSV   = "test1_results.csv"
BETACHA      = "test1.betacha"

QUIET_SEC    = 2.0
MAX_WAIT_SEC = 300

OUT_COLUMNS = [
    "Ab","Zb","Qb","Eb",
    "At","Zt","thick","density",
    "QM","QF","sig(QF)","iQmax",
    "EnergyResidual",
    "EqThick-CS","EqThick-Slope",
    "finalFile",
]

PLACEHOLDERS = ["Ap","Zp","q","E","At","Zt","Thick","d"]

def read_betacha(path: Path):
    rows = []
    for ln in path.read_text().splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("!"): continue
        parts = ln.split()
        rows.append(dict(zip(PLACEHOLDERS, parts)))
    return rows

def create_case_files(betacha_rows, betacha_path: Path):
    base = betacha_path.stem
    case_files = []
    template_text = Path("template.etacha").read_text()

    outdir = Path("etacha_cases")
    outdir.mkdir(exist_ok=True)

    for i, params in enumerate(betacha_rows, 1):
        prefix = f"{base}_case{i:03d}"
        case_path = outdir / (prefix + ".etacha")
        txt = template_text
        for k in PLACEHOLDERS:
            txt = txt.replace("{{"+k+"}}", str(params[k]))
        case_path.write_text(txt)
        case_files.append((prefix, case_path, params))
    return case_files

def kill_etacha():
    subprocess.run(["taskkill", "/IM", "ETACHA4.exe", "/F"], capture_output=True)

def snap_results():
    snap = {}
    for p in RESULTS_DIR.glob("*.txt"):
        try:
            st = p.stat()
            snap[p] = (st.st_size, st.st_mtime)
        except: pass
    return snap

def changed_since(prev):
    changed = set()
    for p in RESULTS_DIR.glob("*.txt"):
        try:
            st = p.stat()
            sig = (st.st_size, st.st_mtime)
        except: continue
        if prev.get(p) != sig:
            changed.add(p)
    return changed

def choose_final(prefix):
    f1 = list(RESULTS_DIR.glob(f"{prefix}_forExcel.txt"))
    if f1: return f1[0]
    etas = []
    for p in RESULTS_DIR.glob(f"{prefix}_Eta*.txt"):
        m = re.search(r"_Eta(\d+)\.txt$", p.name)
        if m: etas.append((int(m.group(1)), p))
    if etas:
        etas.sort()
        return etas[-1][1]
    return None

def parse_forExcel(path):
    txt = path.read_text(errors="ignore")
    lines = [ln for ln in txt.splitlines() if ln.strip()]
    header = re.split(r"\s+", lines[0])
    last   = re.split(r"\s+", lines[-1])
    row = dict(zip(header, last))
    energy = float(row["Energy"]) if "Energy" in row else None
    qmap = {}
    for h,val in row.items():
        if h.startswith("q") and re.fullmatch(r"q\d+", h):
            try: qmap[int(h[1:])] = float(val)
            except: pass
    s = sum(qmap.values())
    if s > 1.5: qmap = {q:v/100.0 for q,v in qmap.items()}
    return energy, qmap

def parse_eta(path):
    qmap = {}
    for ln in path.read_text(errors="ignore").splitlines():
        parts = ln.split()
        if len(parts)>=2 and parts[0].isdigit():
            try: qmap[int(parts[0])] = float(parts[1])
            except: pass
    s = sum(qmap.values())
    if s > 1.5: qmap = {q:v/100.0 for q,v in qmap.items()}
    return None, qmap

def compute_moments(qmap):
    if not qmap: return None,None,None,None
    s = sum(qmap.values())
    if s == 0: return None,None,None,None
    qmap = {q:v/s for q,v in qmap.items()}
    qm = sum(q*p for q,p in qmap.items())
    var= sum(p*(q-qm)**2 for q,p in qmap.items())
    sig= var**0.5
    qf = max(qmap, key=qmap.get)
    return qm, qf, sig, qf

def run_case(prefix, case_path):
    print(f"Running {prefix}")
    before = snap_results()
    proc = subprocess.Popen([ETACHA_EXE, str(case_path), "-r"])
    start = last = time.time()
    saw = False

    while True:
        time.sleep(0.5)
        now = time.time()
        changed = [p for p in changed_since(before) if p.name.startswith(prefix+"_")]
        if changed:
            saw = True
            last = now
            for p in changed:
                try: st=p.stat(); before[p]=(st.st_size, st.st_mtime)
                except: pass
        if saw and now-last>=QUIET_SEC: break
        if proc.poll() is not None: break
        if now-start>MAX_WAIT_SEC:
            print("timeout")
            break

    kill_etacha()
    return choose_final(prefix)

def main():
    if not Path(ETACHA_EXE).exists():
        print("ETACHA not found"); return
    if not Path("template.etacha").exists():
        print("Missing template.etacha"); return

    betacha_path = Path(BETACHA)
    rows_in = read_betacha(betacha_path)
    case_files = create_case_files(rows_in, betacha_path)

    results = []
    for (prefix, case_path, params) in case_files:
        final_path = run_case(prefix, case_path)
        energy_res = None; qmap = {}

        if final_path:
            if final_path.name.endswith("_forExcel.txt"):
                energy_res, qmap = parse_forExcel(final_path)
            else:
                energy_res, qmap = parse_eta(final_path)

        QM,QF,SIG,IQMAX = compute_moments(qmap)

        results.append({
            "Ab":params["Ap"],"Zb":params["Zp"],"Qb":params["q"],"Eb":params["E"],
            "At":params["At"],"Zt":params["Zt"],
            "thick":params["Thick"],"density":params["d"],
            "QM":QM,"QF":QF,"sig(QF)":SIG,"iQmax":IQMAX,
            "EnergyResidual":energy_res,
            "EqThick-CS":"","EqThick-Slope":"",
            "finalFile":str(final_path) if final_path else ""
        })

    with open(OUTPUT_CSV,"w",newline="",encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=OUT_COLUMNS)
        w.writeheader()
        for r in results:
            w.writerow(r)

    print(f"\nDONE → {OUTPUT_CSV}")

if __name__ == "__main__":
    main()
