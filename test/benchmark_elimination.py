import time, sys, glob
import numpy as np

sys.path.insert(0,'/Users/cwagen/code/cctk')
import cctk

#### 80.6 seconds before refactoring (6/6/2020) (I think)
#### 26.7 seconds after refactoring (6/6/2020)

files = glob.glob("/Users/cwagen/code/Martin/ts/*.out")

w_start = time.time()
p_start = time.process_time()

e = cctk.ConformationalEnsemble()
for path in files:
    f = cctk.GaussianFile.read_file(path)
    m = f.get_molecule().assign_connectivity()
    if len(e):
        m = m.renumber_to_match(e.molecules[0])
    e.add_molecule(m)

w_end = time.time()
p_end = time.process_time()
print(f"Elapsed time {w_end-w_start:.2f} s (CPU: {p_end-p_start:.2f} s)")

#### 18.6 seconds before refactoring (6/6/2020)
#### 3.88 seconds after refactoring (6/6/2020)

w_start = time.time()
p_start = time.process_time()

print("eliminating... (50x)")
for n in range(50):
    e.eliminate_redundant()

w_end = time.time()
p_end = time.process_time()
print(f"Elapsed time {w_end-w_start:.2f} s (CPU: {p_end-p_start:.2f} s)")

