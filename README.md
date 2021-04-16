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



## 4) Solvation: solvating with water

```
gmx solvate -cp 1aki_newbox.gro -cs spc216.gro -o 1aki_solv.gro -p topol.top
```

The configuration of the protein (-cp) is contained in the output of the previous editconf step, and the configuration of the solvent (-cs) is part of the standard GROMACS installation. 



## 5) Adding ion: assembling tpr file

```
gmx grompp -f ions.mdp -c 1aki_solv.gro -p topol.top -o ions.tpr -maxwarn 5
```

It will generate an atomic-level description of the system in the binary file ions.tpr.



## 6) Adding ion: adding ions to neutralize

In CLI:

```
gmx genion -s ions.tpr -o 1aki_solv_ions.gro -p topol.top -pname NA -nname CL -neutral
```

After prompting 13 Sol was selected.

In Jupyter Notebook:

```
printf "SOL" | gmx genion -s ions.tpr -o 1aki_solv_ions.gro -p topol.top -pname NA -nname CL -neutral
```

In the genion command, we provide the structure/state file (-s) as input, generate a .gro file as output (-o), process the topology (-p) to reflect the removal of water molecules and addition of ions, define positive and negative ion names (-pname and -nname, respectively), and tell genion to add only the ions necessary to neutralize the net charge on the protein by adding the correct number of negative ions (-neutral, which in this case will add 8 Cl- ions to offset the +8 charge on the protein).





