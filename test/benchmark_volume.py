import time, sys

sys.path.insert(0,'/Users/cwagen/code/cctk')
import cctk

#### 71.3 seconds (6/11/2020)

files = [
    "test/static/eliminationTS.out",
    "test/static/LSD_custom.out",
    "test/static/diiron_complex.out",
    "test/static/ibuprofen_solvated.out",
    "test/static/ibuprofen_solvated2.out",
    "test/static/gaussian_file.out",
    "test/static/acetone-couplings1.out",
    "test/static/acetone-couplings2.out",
    "test/static/acetone-couplings4.out",
    "test/static/acetone-couplings5.out",
    "test/static/acetone-couplings6.out",
    "test/static/pentane_conformation_1.out",
    "test/static/pentane_conformation_2.out",
    "test/static/pentane_conformation_3.out",
    "test/static/pentane_conformation_4.out",
    "test/static/glycosylation_TS.out",
]

mols = list()

for f in files:
    file = cctk.GaussianFile.read_file(f)
    if isinstance(file, list):
        file = file[-1]
    mols.append(file.get_molecule())

print("start")

w_start = time.time()
p_start = time.process_time()

for m in mols:
    m.volume()

w_end = time.time()
p_end = time.process_time()
print(f"Elapsed time {w_end-w_start:.2f} s (CPU: {p_end-p_start:.2f} s)")
