import time, sys
import numpy as np

sys.path.insert(0,'/Users/cwagen/code/cctk')
import cctk

#### 24.8 seconds before refactoring (6/5/2020)
#### 0.63 seconds after refactoring (6/5/2020)

file = "test/static/glycosylation_TS.out"
mol = cctk.GaussianFile.read_file(file).get_molecule()

w_start = time.time()
p_start = time.process_time()

for n in range(100):
    pass
    mol.assign_connectivity()

w_end = time.time()
p_end = time.process_time()

print(f"Elapsed time {w_end-w_start:.2f} s (CPU: {p_end-p_start:.2f} s)")

#### 29.0 seconds before refactoring (6/5/2020)
#### 0.89 seconds after refactoring (6/6/2020)

w_start = time.time()
p_start = time.process_time()

for n in range(100):
    mol.assign_connectivity(periodic_boundary_conditions=np.array([100, 100, 100]))

w_end = time.time()
p_end = time.process_time()

print("Periodic boundary conditions:")
print(f"Elapsed time {w_end-w_start:.2f} s (CPU: {p_end-p_start:.2f} s)")
