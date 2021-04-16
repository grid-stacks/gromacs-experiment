# Gromacs Experiment



## 1) Cleaning protein (removing crystal water)

```
grep -v HOH 1aki.pdb > 1aki_cleaned.pdb
```



## 2) Generating topology (pdb2gmx)

```
gmx pdb2gmx -f 1aki_cleaned.pdb -o 1aki_processed.gro -water spce
```

After running this command force field list will be promted (15 OPLS-AA/L selected).

The purpose of pdb2gmx is to generate three files:

- The topology for the molecule.
- A position restraint file.
- A post-processed structure file.



## 3) Solvation: defining box 

```
gmx editconf -f 1aki_processed.gro -o 1aki_newbox.gro -c -d 1.0 -bt cubic
```

The above command centers the protein in the box (-c), and places it at least 1.0 nm from the box edge (-d 1.0). The box type is defined as a cube (-bt cubic). 



