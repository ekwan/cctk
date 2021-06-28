import cctk, sys, tqdm

# usage: python generate_si.py config.txt output.txt
# config should be tab-separated list of titles and filenames

ensemble = cctk.Ensemble()
titles = list()

with open(sys.argv[1], "r+") as f:
    for line in tqdm.tqdm(f.readlines()):
        try:
            title, path = line.rsplit(' ', 1)
            gf = cctk.GaussianFile.read_file(path.strip())
            m, p = gf.get_molecule(properties=True)
            p["route_card"] = gf.route_card
            p["imaginaries"] = gf.imaginaries()
            ensemble.add_molecule(m, p)
            titles.append(title.strip())
        except Exception as e:
            print(f"skipping {line} - error {e}")

si_file = cctk.SIFile(ensemble, titles)
si_file.write_file(sys.argv[2])
