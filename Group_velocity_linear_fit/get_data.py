import yaml
import numpy as np

INFILE  = "band.yaml"
OUTFILE = "qdistance_freq_p-43m_GX_20.txt"
N_POINTS = 20

with open(INFILE, "r") as f:
    data = yaml.safe_load(f)

rows = []
for q in data["phonon"][:N_POINTS]:
    dist = q["distance"]
    freqs = [b["frequency"] for b in q["band"][:3]]   # first 3 bands
    rows.append([dist] + freqs)

rows = np.array(rows)

np.savetxt(
    OUTFILE,
    rows,
    fmt="%.8f",
    header=f"{'distance':>15} {'band1':>15} {'band2':>15} {'band3':>15}",
)

print(f"Wrote {len(rows)} rows to {OUTFILE}")
