#### Tutorial 4: Analyzing Charges and NICS Along the Reaction Coordinate

Analysis of higher-order properties along a potential energy surface is frequently complicated by the multitude of jobs required. 
For instance, analyzing the change in atomic charges or aromaticity (as measured by nucleus-independent chemical shift, or NICS)
requires a separate job for each desired point, making manual analysis challenging. 

##### Step 1: Finding Points Along the Intrinsic Reaction Coordinate

The transition state for the reaction was found using conventional techniques (scanning the C1–C5 bond distance) 
and confirmed with a frequency calculation (*v*<sub>i</sub> = -1336 cm<sup>-1</sup>: 

<img src='TS.png' width='350px'>

The intrinsic reaction coordinate was followed backwards and forwards for 50 steps using the following input line:

```
#p irc=(calcfc, forward, maxpoints=50, stepsize=2) m062x/6-31g(d)
```

This resulted in the generation of two IRC `.out` files, each containing 51 distinct structures. 

##### Step 2: Generating Jobs to Calculate NICS(0)/Hirschfeld Charges (`generate_nics.py`)

##### Step 3: Analysis (`analyze_nics.py`)

