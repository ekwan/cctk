import cctk

file = cctk.XYZFile.read_file("water.xyz")

print(f"{file.title} is title")
print(f"{file.molecule.atomic_numbers} is atomic numbers")

angle = file.molecule.get_angle(2,1,3)
print(f"The angle between atoms 2, 1, and 3 is {angle} degrees")

if angle < 105:
    print("This is a strained water molecule!")
