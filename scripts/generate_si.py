import cctk

ensemble = cctk.GaussianFile.read_file("test/static/gaussian_file.out").ensemble
si_file = cctk.SIFile(ensemble, ["title"] * len(ensemble))

si_file.write_file("si.txt")
