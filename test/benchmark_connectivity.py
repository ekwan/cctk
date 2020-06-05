import time, sys

sys.path.insert(0,'/Users/cwagen/code/cctk')
import cctk

#### 24.8 seconds before refactoring (6/5/2020)
#### 0.63 seconds after refactoring (6/5/2020)

file = "test/static/glycosylation_TS.out"
mol = cctk.GaussianFile.read_file(file).get_molecule()

w_start = time.time()
p_start = time.process_time()

for n in range(100):
    mol.assign_connectivity()

w_end = time.time()
p_end = time.process_time()

print(f"Elapsed time {w_end-w_start:.2f} s (CPU: {p_end-p_start:.2f} s)")
