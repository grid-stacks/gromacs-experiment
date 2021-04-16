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



## 7) Energy minimization: assembling the binary input

```
gmx grompp -f minim.mdp -c 1aki_solv_ions.gro -p topol.top -o em.tpr
```



## 8) Energy minimization: carry out the EM

```
gmx mdrun -ntmpi 2 -ntomp 1 -deffnm em -v -pin on [-nb gpu]
```

-nb gpu is for running in GPU.

We will get 4 files:

- em.log: ASCII-text log file of the EM process
- em.edr: Binary energy file
- em.trr: Binary full-precision trajectory
- em.gro: Energy-minimized structure



## 9) Analyzing the energy terms

In CLI:

```
gmx energy -f em.edr -o potential.xvg
```

After prompting "10 0" was selected.

In Jupyter Notebook:

```
printf "10 0" | gmx energy -f em.edr -o potential.xvg
```



## 10) Plotting potential file:

```
%matplotlib inline
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
```

```
potential = np.genfromtxt([i for i in open('potential.xvg').read().splitlines() if not i.startswith(('#','@'))])
```

```
plt.plot(*potential.T)
plt.xlabel('step')
plt.ylabel('potential')
```

![img](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAaAAAAEGCAYAAAAjc0GqAAAAOXRFWHRTb2Z0d2FyZQBNYXRwbG90bGliIHZlcnNpb24zLjMuMiwgaHR0cHM6Ly9tYXRwbG90bGliLm9yZy8vihELAAAACXBIWXMAAAsTAAALEwEAmpwYAAAodklEQVR4nO3deZhcVZ3/8fe3qnrf0uk02TqrhCA7JIaIOrIJGTfQAY2OQ0Z5figuM+osjqKDA+jILGTAGVFGUURHUHDBhYGwjgsGOgQIhCUhhNBZu9Od9L5/f3/cU6G60510SKpuL5/X89ynb517zulTl5BvznLPNXdHREQk1xJxN0BERCYmBSAREYmFApCIiMRCAUhERGKhACQiIrFIxd2AsWTKlCk+d+7cuJshIjKmrFmzpsHdqwenKwAdgrlz51JbWxt3M0RExhQze3modA3BiYhILBSAREQkFgpAIiISCwUgERGJhQKQiIjEQgFIRERioQAkIiKxUADKgZ+treMHfxxyGbyIyISlAJQDv3xyO7c9tiXuZoiIjCoKQDmQn0zQ3dsfdzNEREYVBaAcyE8l6FIAEhEZQAEoBwpS6gGJiAymAJQD+QpAIiL7UQDKAQUgEZH9KQDlQH4qQVefApCISCYFoBwoCKvg3D3upoiIjBoKQDmQn4puc0+fApCISJoCUA6kA1BXb1/MLRERGT0UgHKgIJUE0EIEEZEMCkA5kO4BdWshgojIPrEEIDO72syeMrMnzOxeM5sR0ueaWUdIf8LMvplRZpGZrTOzjWZ2g5lZSC8ws9tD+mozm5tRZoWZbQjHioz0eSHvhlA2P5vfNz8ZApB6QCIi+8TVA/pXdz/J3U8BfgX8Y8a1F939lHB8LCP9RuAyYEE4loX0S4Emdz8aWAlcC2Bmk4ErgdOBJcCVZlYZylwLrHT3BUBTqCNrUkkDtAhBRCRTLAHI3ZszPpYAB/yb2cymA+Xu/ohHa5m/D1wYLl8A3BLO7wDOCb2j84FV7t7o7k3AKmBZuHZ2yEsom64rK/JCD6ivXwFIRCQttjkgM/uKmb0C/DkDe0DzzGytmT1sZm8JaTOBuow8dSEtfe0VAHfvBfYCVZnpg8pUAXtC3sF1DdXOy8ys1sxq6+vrX8M3hWQi3QPSEJyISFrWApCZ3WdmTw9xXADg7le4+yzgh8AnQ7HtwGx3PxX4LPA/ZlYO2BC/It2dGO7aoaYPyd1vcvfF7r64urp6uGwHlBeG4NQDEhF5VSpbFbv7uSPM+j/Ar4Er3b0L6Arl15jZi8AxRL2UmowyNcC2cF4HzALqzCwFVACNIf3MQWUeAhqASWaWCr2gzLqyIpmI4nxvv3pAIiJpca2CW5Dx8d3AcyG92syS4Xw+0WKDTe6+HWgxs6VhDucS4Beh/F1AeoXbRcADYZ7oHuA8M6sMiw/OA+4J1x4MeQll03VlRV5CixBERAbLWg/oIL5mZguBfuBlIL3a7U+Aq8ysF+gDPubujeHa5cD3gCLg7nAAfAe41cw2EvV8lgO4e6OZXQ08FvJdlVHX54DbzOwaYG2oI2vSc0AaghMReVUsAcjd/2yY9DuBO4e5VgucMER6J3DxMGVuBm4eIn0T0dLsnEgl03vBaQhORCRNOyHkgBYhiIjsTwEoB5KaAxIR2Y8CUA7oQVQRkf0pAOVAugekZdgiIq9SAMqBvIReSCciMpgCUA4k9y1CUA9IRCRNASgH9CCqiMj+FIByQA+iiojsTwEoB/QgqojI/hSAckAPooqI7E8BKAdeXYatACQikqYAlAPpZdi9WoQgIrKPAlAOJBKGmR5EFRHJpACUI3mJhJZhi4hkUADKkWTC9CCqiEgGBaAcSSVNPSARkQwKQDmSl0xoGbaISAYFoBxJJkyLEEREMigA5UhewrQMW0QkgwJQjiSTpgdRRUQyKADlSLQMW0NwIiJpCkA5Ei3DVg9IRCRNAShHUkk9iCoikkkBKEdSehBVRGQABaAcSWkRgojIAApAOZKXSGgZtohIBgWgHNGDqCIiAykA5YiG4EREBlIAypGUdkIQERkg1gBkZn9rZm5mUzLSPm9mG83seTM7PyN9kZmtC9duMDML6QVmdntIX21mczPKrDCzDeFYkZE+L+TdEMrmZ/u7RsuwNQQnIpIWWwAys1nA24AtGWnHAcuB44FlwDfMLBku3whcBiwIx7KQfinQ5O5HAyuBa0Ndk4ErgdOBJcCVZlYZylwLrHT3BUBTqCOrUgkNwYmIZIqzB7QS+Hsg82/lC4Db3L3L3V8CNgJLzGw6UO7uj7i7A98HLswoc0s4vwM4J/SOzgdWuXujuzcBq4Bl4drZIS+hbLqurEklE/QrAImI7BNLADKzdwNb3f3JQZdmAq9kfK4LaTPD+eD0AWXcvRfYC1QdoK4qYE/IO7iuodp6mZnVmlltfX39iL/jYKmE0aNVcCIi+6SyVbGZ3QdMG+LSFcAXgPOGKjZEmh8g/bWUOVBd+19wvwm4CWDx4sWvuQuTShh9WoQgIrJP1gKQu587VLqZnQjMA54M6whqgMfNbAlRb2RWRvYaYFtIrxkinYwydWaWAiqAxpB+5qAyDwENwCQzS4VeUGZdWaNl2CIiA+V8CM7d17n7Ue4+193nEgWK09x9B3AXsDysbJtHtNjgUXffDrSY2dIwh3MJ8ItQ5V1AeoXbRcADYZ7oHuA8M6sMiw/OA+4J1x4MeQll03VlTVKLEEREBshaD+i1cPdnzOzHwHqgF/iEu/eFy5cD3wOKgLvDAfAd4FYz20jU81ke6mo0s6uBx0K+q9y9MZx/DrjNzK4B1oY6siqVSNCrZdgiIvvEHoBCLyjz81eArwyRrxY4YYj0TuDiYeq+Gbh5iPRNREuzcyal9wGJiAygnRByRK/kFhEZSAEoR/QgqojIQApAOZJKJOjrd6I1ECIiogCUI6lE9PiR5oFERCIKQDmSTEYBSMNwIiIRBaAcyUtEt1oBSEQkogCUI8n0EJy24xERARSAciYVhuC0IamISEQBKEeSWoQgIjKAAlCOaA5IRGQgBaAc0RyQiMhACkA5ojkgEZGBFIByJBWG4DQHJCISUQDKkfQQXK+G4EREAAWgnElvxdOrITgREUABKGe0FY+IyEAKQDmSpzkgEZEBFIByJD0H1KPXcouIAApAOZNehq0ekIhIRAEoR15dhKAAJCICCkA5s+85IC3DFhEBIHWgi2b22QNdd/frjmxzxq+klmGLiAxwwAAElOWkFRNAnpZhi4gMcMAA5O7/lKuGjHd6HYOIyEAH6wEBYGaFwKXA8UBhOt3dP5Kldo076TmgHs0BiYgAI1+EcCswDTgfeBioAVqy1ajxKLlvGbbmgEREYOQB6Gh3/xLQ5u63AO8ATsxes8afPC3DFhEZYKQBqCf83GNmJwAVwNystGic0m7YIiIDjWgOCLjJzCqBLwF3AaXAP2atVeNQSq/kFhEZYEQByN2/HU4fBuZnrznjV0pzQCIiAxxwCM7MPhR+fnao43B/uZn9rZm5mU0Jn+eaWYeZPRGOb2bkXWRm68xso5ndYGYW0gvM7PaQvtrM5maUWWFmG8KxIiN9Xsi7IZTNP9zvcjBJzQGJiAxwsDmgkvCzbIij9HB+sZnNAt4GbBl06UV3PyUcH8tIvxG4DFgQjmUh/VKgyd2PBlYC14b6JwNXAqcDS4ArwzAiIc9Kd18ANIU6siqlOSARkQEO9iDqt8Lpfe7++8xrZvamw/zdK4G/B35xsIxmNh0od/dHwufvAxcCdwMXAF8OWe8A/jP0js4HVrl7YyizClhmZrcBZwMfDGVuCeVvPMzvc0DqAYmIDDTSVXBfH2HaiJjZu4Gt7v7kEJfnmdlaM3vYzN4S0mYCdRl56kJa+torAO7eC+wFqjLTB5WpAvaEvIPrGqqtl5lZrZnV1tfXH8rXHFwPqYRpDkhEJDjYZqRvBM4AqgfN+ZQDyYOUvY/o4dXBrgC+AJw3xLXtwGx3321mi4Cfm9nxgA2RN92VGO7aoaYPyd1vAm4CWLx48WF1X/JTCTp7FIBERODgq+DyieZ6UgzcmLQZuOhABd393KHSzexEYB7wZFhHUAM8bmZL3H0H0BXKrzGzF4FjiHopNRnV1ADbwnkdMAuoM7MU0TNKjSH9zEFlHgIagElmlgq9oMy6sqqyOJ+m9u5c/CoRkVHvYHNADwMPm9n33P3lI/EL3X0dcFT6s5ltBha7e4OZVQON7t5nZvOJFhtscvdGM2sxs6XAauASXh0CvAtYATxCFBQfcHc3s3uAr2YsPDgP+Hy49mDIe1soe9B5qCNhUnEee9p7Dp5RRGQCGOmDqAVmdhPR7gf7yrj72Ue4PX8CXGVmvUAf8LH0IgLgcuB7QBHR4oO7Q/p3gFvNbCNRz2d5aFujmV0NPBbyXZVR1+eA28zsGmBtqCPr1AMSEXnVSAPQT4BvAt8mCgxHjLvPzTi/E7hzmHy1wAlDpHcCFw9T5mbg5iHSNxEtzc6pypJ86prac/1rRURGpZEGoF53z+oy5YmgsjiPJg3BiYgAI1+G/Usz+7iZTTezyekjqy0bhyYV59Pc2UNvn1bCiYiMtAeU3sbm7zLSHO0Ld0gmF+fhDns7eqgqLYi7OSIisRrpZqTzst2QiaCyJNpyrqldAUhEZERDcGZWbGZfDCvhMLMFZvbO7DZt/JlWHr3N/KHnd8XcEhGR+I10Dui7QDfRrggQPeR5TVZaNI69Ye5kppTm88y25ribIiISu5EGoNe5+78Q3ozq7h0MvaWNHEAiYdRUFtPQ2hV3U0REYjfSANRtZkWEPdPM7HWELXPk0EwpzWd3qx5GFREZaQD6MvC/wCwz+yFwP9FuAnKIJpfkqwckIsLIV8Hda2ZrgKVEQ29/7e4NWW3ZODV7cjG7Wrpo6eyhrDAv7uaIiMRmpKvg7nf33e7+a3f/Vdg49P5sN248Om5GOQDP72iJuSUiIvE6YAAys8Kw48EUM6vM2AVhLjAjJy0cZ14/PQpA67drJZyITGwHG4L7KPBpomDzeEZ6M/BfWWrTuDatvJDK4jyeqtsbd1NERGJ1sPcBXQ9cb2afcvfX/ApueZWZ8ZYF1Tzw3C7cnfBSPhGRCWekq+C+ZWZ/ZWZ3hOOTZqYZ9NfoxJkVNLZ109zZG3dTRERiM9IA9A1gUfiZPtfrGV6jaRXRljw79nbG3BIRkfiMdDfsN7j7yRmfHzCzJ7PRoIlgeghA2/Z0sHBaWcytERGJx0h7QH1h9wMAzGw+R/jNqBPJwmll5KcSPPxCfdxNERGJzUgD0N8BD5rZQ2b2EPAA8DdZa9U4V1aYx1kLq/n1uu309XvczRERicVIA9DvgW8B/eH4FvBItho1Ebz75JnUt3Sx+qXdcTdFRCQWIw1A3wfmAVeHYx5wa7YaNRG8dWE1ZvDYS01xN0VEJBYjXYSwcNAihAe1COHwlBakmFdVwvX3v8BH3zqfwrxk3E0SEcmpkfaA1prZ0vQHMzudaFhODsM7TppOv8Nz2hdORCagkQag04E/mNlmM9tMNP/zVjNbZ2ZPZa1149z73zALgCe2aBhORCaekQ7BLctqKyaomspi5k0p4cHn6/nLN82LuzkiIjk10vcBvZzthkxU5xx7FN/+3Us8vqWJ02ZXxt0cEZGcGekQnGTJh98c9XxWrd8Zc0tERHJLAShmMycVMXtyMXVNHXE3RUQkpxSARoGZk4p4YUcL/doVQUQmkFgCkJl92cy2mtkT4Xh7xrXPm9lGM3vezM7PSF8UVt1tNLMbLLxIx8wKzOz2kL46vK01XWaFmW0Ix4qM9Hkh74ZQNj9HX31IF546g+d3tnDz71+KsxkiIjkVZw9opbufEo7fAJjZccBy4HiilXffMLP0E5o3ApcBC8KRXpl3KdDk7kcDK4FrQ12TgSuJlpAvAa40s/Qs/7Xh9y8AmkIdsXnf4lmctbCa6+/fQHdvf5xNERHJmdE2BHcBcJu7d7n7S8BGYImZTQfK3f0Rd3eirYEuzChzSzi/Azgn9I7OB1a5e6O7NwGrgGXh2tkhL6Fsuq5YmBmXnDGXls5e7n56e5xNERHJmTgD0CfN7CkzuzmjZzITeCUjT11ImxnOB6cPKOPuvcBeoOoAdVUBe0LewXXtx8wuM7NaM6utr8/e6xPeuqCa+dUlfOvhTUQxVkRkfMtaADKz+8zs6SGOC4iG014HnAJsB/49XWyIqvwA6a+lzIHq2v+C+03uvtjdF1dXVw+X7bAlEsZH/2Q+67c38/iWPVn7PSIio8VId0I4ZO5+7kjymdl/A78KH+uAWRmXa4BtIb1miPTMMnVmlgIqgMaQfuagMg8BDcAkM0uFXlBmXbE6c+FRAKyr28OiOXooVUTGt7hWwU3P+Pge4OlwfhewPKxsm0e02OBRd98OtJjZ0jCHcwnwi4wy6RVuFwEPhHmie4DzzKwyDPGdB9wTrj0Y8hLKpuuK1VFlBUwtL+COx+vo7NELZ0VkfItrDuhfMjYyPQv4DIC7PwP8GFgP/C/wCXdP/018OfBtooUJLwJ3h/TvAFVmthH4LPAPoa5GoncXPRaOq0IawOeAz4YyVaGO2JkZV19wAk9vbeaaX6+PuzkiIlllmvAeucWLF3ttbW3Wf88VP1vHD1dv4frlp3DBKcOujxARGRPMbI27Lx6cPtqWYQvwpXcex3HTy7lu1QvaHUFExi0FoFGoMC/Jx858HS/vbueXT42K9REiIkecAtAo9acnTGN+dQk3/07b84jI+KQANErlJRO89ZhqXtjZSm+ftucRkfFHAWgUWzq/io6ePv7f92tp6+o9eAERkTFEAWgUO//4aVxz4Qk89EI9Z/3bQ+zY2xl3k0REjhgFoFHuQ0vn8J0Vi9nV0sWPHt0Sd3NERI4YBaAx4Oxjp/Lmo6fwXw9uZOsevTlVRMYHBaAx4uNnvo7efqd2c+PBM4uIjAEKQGPEaXMqSRj84I8v06oFCSIyDigAjRGFeUlWnDGX2pebeOM/389TdXvibpKIyGFRABpDrnzX8fzko28E4Gt3P6cX14nImKYANMYsnjuZz5x7DH94cTf3Pbsr7uaIiLxmCkBj0IeWzqEgleAnta/Qp81KRWSMUgAag/JTCc5aeBT3rt/Ju77+O9Zva467SSIih0wBaIy64QOnct37Tmbb3g4+/sM1eoOqiIw5CkBjVH4qwXtPq+E/P3Aam3e385M1dXE3SUTkkCgAjXFvXjCF+VNKuPHBjbxY3xp3c0RERkwBaBz48Jvnsauli3Ove5jr7n1eCxNEZExIxd0AOXx/sXQOZy2s5l/veZ4bHthIfWs3X33PCZhZ3E0TERmWAtA4UVNZzPXLT6W8MI9b//gyj760m4+8eR4fXDJbgUhERiUNwY0zX3rncfzbxSeTl0xwxc+e5pM/Wkt3r96oKiKjj3pA40x+KsFFi2r4s9Nm8s93P8dN/7eJzu4+Ljh1JucdN5XCvGTcTRQRARSAxi0z4/N/eiwAd66p4/7ndnFSTQU3/cViplUUxtw6ERENwY1rZsYX3v56Hr3iXK5ffgpP1e3lo7fW8tjmRnr6NCwnIvFSAJoAkgnjglNm8pX3nMCTdXu5+JuP8K6v/44dezvjbpqITGAKQBPIn58+h8e/9Daue9/J1DV18OHvPUZDa1fczRKRCUoBaIKZXJLPe0+r4e/OX8iz25tZfM19XHLzozy+pSnuponIBKNFCBPUijPmcvKsSdz7zA5uf+wVLrrxD7zjpBmc8boq3nXyDEoL9EdDRLLL9FbNkVu8eLHX1tbG3YwjrqG1i6/++ll+t7GBXS1dVJXk89X3nsjZxx5FXlKdZBE5PGa2xt0XD06P5W8XM/uymW01syfC8faQPtfMOjLSv5lRZpGZrTOzjWZ2g4XH+82swMxuD+mrzWxuRpkVZrYhHCsy0ueFvBtC2fwcfv1RZ0ppAde9/xRWf+EcfvrxM8hPJfjorWs4/h/v4cPffZRXGtvjbqKIjENx/vN2pbufEo7fZKS/mJH+sYz0G4HLgAXhWBbSLwWa3P1oYCVwLYCZTQauBE4HlgBXmlllKHNt+P0LgKZQx4RnZpw2u5JffurNrHz/yaw4Yw6PbW7i7Tf8lv9ZvYXndujFdyJy5IyJ8RUzmw6Uu/sjHo0Zfh+4MFy+ALglnN8BnBN6R+cDq9y90d2bgFXAsnDt7JCXUDZdlxD1iN5zag1XvOM4bvnIEhJmfOFn61j2H7/lz7/9R17Y2aIdt0XksMU50/xJM7sEqAX+JgQJgHlmthZoBr7o7r8FZgKZb1yrC2mEn68AuHuvme0FqjLTB5WpAva4e+8Qde3HzC4j6nkxe/bs1/hVx65Fcyp5/Etvo66pnbue2MZ1973AeSv/j7KCFBctruHT5xxDRXFe3M0UkTEoawHIzO4Dpg1x6Qqi4bSrAQ8//x34CLAdmO3uu81sEfBzMzseGGo75/Q/wYe7dqjpQ3L3m4CbIFqEMFy+8SyZMOZUlfCpcxaw7IRpPPHKHn67oYHv/n4z3/39ZuZNKeFDS+ew4o1zSGnRgoiMUNYCkLufO5J8ZvbfwK9CmS6gK5yvMbMXgWOIeik1GcVqgG3hvA6YBdSZWQqoABpD+pmDyjwENACTzCwVekGZdclBLJhaxoKpZVy8eBYfPH02j77UyO82NnD1r9bzjQc3clJNBWcfexQXLZpFUb42PhWR4cWyDNvMprv79nD+GeB0d19uZtVAo7v3mdl84LfAie7eaGaPAZ8CVgO/Ab7u7r8xs0+EPB8zs+XAe939fWERwhrgtPBrHwcWhbp+Atzp7reFlXZPufs3Dtbu8boM+0i468ltPPDsTp7aupdN9W2kEsbsycXMry5l+Rtmcc7rj9J7iUQmqOGWYccVgG4FTiEa+toMfNTdt5vZnwFXAb1AH3Clu/8ylFkMfA8oAu4GPuXubmaFwK3AqUQ9n+XuvimU+QjwhfBrv+Lu3w3p84HbgMnAWuBDofd1QApAI/OHFxt4+Pl6Xmlqp3ZzE7taujh+RjlvWVDN+cdP5bgZ5RSk1DsSmShGVQAaqxSADl1Xbx/X3fsCj2zazfptzfT2O2YwrbyQs489ik+efTTTK4ribqaIZJEC0BGgAHR4mtq6efD5XWxpbGfDzlZWrd9Jd18/0ysKOX5GNHe0dP5k5laVkEhouE5kvBguAGnDL8mZyrARatrzO1r45ZPb2Lang9UvNXLfszsBKC1I8bbjpvLG+VXMqy7hpJoKDdmJjEMKQBKbhdPKWDhtIQDuzov1rTz+8h4e3dzIz9du5WdrtwJQXpjihJkVzJhUxJK5k1l24jTKC/XskchYpyG4Q6AhuNzp6u1jV3MX67c3c/e67WxpbGdLYzsNrd0AlBWmOGZqGRcvquH0+VVMryikME+9JJHRSHNAR4ACULzcnUdfauSxzY3Ut3Txx02NPL+zBYgelj1xZgWnzJrElNJ85k0p5YSZ5cyYVKQdvUVipjkgGfPMjNPnV3H6/CoA+vudxzY3UtfUwUsNbTzw3C7uXFNHS1fvvjLJhFFTWcSCo0qZMamI2ZOLmVtVwtwpJcyeXEx+SsFJJC7qAR0C9YDGhs6ePp7d3syGXa1s2d3OpoZWNtW3sXVPBy2drwanhEFNZTEnzqxg6euqmFpWwJyqKDBpFweRI0c9IJkwCvOSnDq7klNnV+53ramtm5d2t7G5ITpebGjjwed28et12wfkm1Scx7TyQk6cWcHU8kKqSvM5qWYSx88o11yTyBGiACQTSmVJPpUl+ZyWEZx6+vppbOtmx95OXm5sZ8vuNnY0d/Ly7nbuXb+Tls4e0m+fSBhMLS9k9uRiTpk1iUnF+cyYVMisycXMqixmSmm+thwSGSEFIJnw8pIJppYXMrW8kJNnTdrvurtT39LF2lf28MzWvWzd08m6rXv47u83093XPyBvcX6SOVUlzK8uYV5VCQunlTFvSglVpflMKS3QggiRDJoDOgSaA5JM7k5bdx/b9nRQ19TOK40dbA7De5sa2qhr6hjw4r70FkTTK6JgN786mm+qLivgqLJCjp1WptdZyLikOSCRI8zMKC2Inkc6ZmrZftd7+vpZv62Z7Xs7oyG+5k62NnWwo7mD53e2cO/6nQMCVEEqwZTSAqpK86kuLWBKaQFTyvJDWgFTQnpVaQGTivK0XZGMeQpAIlmSl0xw8qxJnDxr6Ovdvf3saumkobWbTfWtPLOtmab2bna3drNtbyfrtu5ld1v3kK8/TyWMySX5AwLTlLICqkryQ+CKzqvLCphckq+hPxmVFIBEYpKfSlBTWUxNZbSg4b2n7Z+nv9/Z29FDQ2sX9a1dNLR2s7u1i4bWLhpautnd1kV9azeb6ttoaO2iq7d//0qAOVXFzJxUREVRHjMmFVFTGZ1PKS1gekUhFUV5lBflaYWf5JQCkMgolkjYvpV7C4YY5suUnpNqaOmKAlNLdxS4Wrr29a427Grl/md37bd4Iq2iKFp+PrWikGnl0dxU1NPKj36WRMOC1aUFWu0nh00BSGScSM9JlRakmDulZNh8vX39tHT2srejh+17O6lv7WJvRw9727vZ2dzFjuZOdjZ38tz2ZhpauxhiBJD8VILywjzKC1OUFeVRUZTHpPAz8yjP/Fwc/SzJTyp4CaAAJDLhpJKJfb2qAwUqeHUIcHdbNPSXXkyxo7mT5o5eWjp79gWvLbvbovOOniGD1r7fn7AhA9TkknxmTS6mtCBJWWEe0ysKKSvMo6wwxZTSApJadDHuKACJyLAyhwCPPqp0RGXcndauqIe1p72H5hCUhjv2tHfz8u426lu6aOvuG7JOMygrSFFelEdVaQHVpQWUFab2Ba7i/CQlBSkqi/M5qryAsoIUk4rz9WDwKKcAJCJHlJmFnkseNfvvhjSs/n6npauX1q5e9rR3s7O5k5bOXpo7e9mVPu/oYWdLJ1v3dNDa1UNTWw+tGZvPDpYwKMpLUlyQYlLocZUVpihLDx8W5jG5JI/JJQX7hi9LCpLRebiuIcPsUQASkVEhkTE0N3NSEcfPqBhRuZ6+ftq7+2jv7mVncxdNbd20dvXS0Bqdt3X30dbVy572Hlq6etjd2s3mhrYQ3Hro6Tvww/gJi97Smx4OTAewssLUgPTywiholeRH6UWhV1acn6QkP0VJQUq7rw+iACQiY1peMkFFUYKKojymVxQdUtn0cGFjWzctnb20dfXS1t1La1cfrZ3RHFdLZ9Qraw7nLZ097Grp5MX63n2fDxbEXm2rUVIQBan0sGFJwasBqjg/6n1NLsmnoiiP4oIUxXlJikOvrKQgRVnBq3nHes9MAUhEJqzM4cLD0dnTty9QtYVhxI7uPtq6Q1Drinpo6d5YW1ffvmDX3t3H7tb26Lyrj9au3mGf5xrY9tAzK0hRmJ+kOD9JcV7U8yrKiz4XhfSivCRFIegNTEtSnE7Pe/VaYSqZk502FIBERA5TYV6Swrwk1WUFR6S+5s4eWjt79w0tZgasdE+ttevV8/aePjq6o2NPRw/b93bQ3t1HZ08f7d19dPT0cajbfhbmJSjOT+0LZv99yeKDrpo8VApAIiKjTPSM1eH1yjK5O509/bR399IRglV7ONJBavC1jp6Q1t1PR08vxVl4SaMCkIjIOGdm0dDcKHvTr5ZkiIhILBSAREQkFgpAIiISCwUgERGJRWwByMw+ZWbPm9kzZvYvGemfN7ON4dr5GemLzGxduHaDhSewzKzAzG4P6avNbG5GmRVmtiEcKzLS54W8G0LZ/Bx9bRERCWIJQGZ2FnABcJK7Hw/8W0g/DlgOHA8sA75hZullGzcClwELwrEspF8KNLn70cBK4NpQ12TgSuB0YAlwpZmld6a6Fljp7guAplCHiIjkUFw9oMuBr7l7F4C77wrpFwC3uXuXu78EbASWmNl0oNzdH3F3B74PXJhR5pZwfgdwTugdnQ+scvdGd28CVgHLwrWzQ15C2XRdIiKSI3EFoGOAt4RhsIfN7A0hfSbwSka+upA2M5wPTh9Qxt17gb1A1QHqqgL2hLyD69qPmV1mZrVmVltfX3/IX1RERIaWtQdRzew+YNoQl64Iv7cSWAq8Afixmc0Hhtp8yA+Qzmsoc6C69r/gfhNwE4CZ1ZvZy8PlPYgpQMNrLDue6b7sT/dkaLovQxsL92XOUIlZC0Dufu5w18zscuCnYTjtUTPrJ7qJdcCsjKw1wLaQXjNEOhll6swsBVQAjSH9zEFlHiL6DzXJzFKhF5RZ18G+U/VI8g3FzGrdffFrLT9e6b7sT/dkaLovQxvL9yWuIbifE83DYGbHAPlEgeEuYHlY2TaPaLHBo+6+HWgxs6VhDucS4BehrruA9Aq3i4AHQmC7BzjPzCrD4oPzgHvCtQdDXkLZdF0iIpIjce0FdzNws5k9DXQDK0JgeMbMfgysB3qBT7h7+h29lwPfA4qAu8MB8B3gVjPbSNTzWQ7g7o1mdjXwWMh3lbs3hvPPAbeZ2TXA2lCHiIjkkPmh7tEtr4mZXRbmkySD7sv+dE+GpvsytLF8XxSAREQkFtqKR0REYqEAJCIisVAAyjIzWxb2tdtoZv8Qd3tyycxmmdmDZvZs2PPvr0P6ZDNbFfbiW5WxRdKwewGON2aWNLO1Zvar8Fn3xGySmd1hZs+FPzNv1H0BM/tM+P/naTP7kZkVjpf7ogCURWEfu/8C/hQ4DvhA2O9uougF/sbdX0/00PEnwvf/B+D+sBff/eHzwfYCHG/+Gng247PuCVwP/K+7HwucTHR/JvR9MbOZwF8Bi939BCBJ9L3HxX1RAMquJcBGd9/k7t3AbUR7100I7r7d3R8P5y1Ef6HMZOD+fZl78Q25F2BOG50DZlYDvAP4dkbyRL8n5cCfEB6JcPdud9/DBL8vQQooCg/aFxM9OD8u7osCUHYNtx/dhGPRazJOBVYDU8PDxYSfR4VsE+V+/Qfw90B/RtpEvyfzgXrgu2Fo8ttmVsIEvy/uvpXobQFbgO3AXne/l3FyXxSAsuuQ9p0br8ysFLgT+LS7Nx8o6xBp4+p+mdk7gV3uvmakRYZIG1f3JEgBpwE3uvupQBthWGkYE+K+hLmdC4B5wAygxMw+dKAiQ6SN2vuiAJRdw+1tN2GYWR5R8Pmhu/80JO8Mr9gg/Ey/jmMi3K83Ae82s81EQ7Jnm9kPmNj3BKLvWefuq8PnO4gC0kS/L+cCL7l7vbv3AD8FzmCc3BcFoOx6DFhg0RtY84kmB++KuU05E/bt+w7wrLtfl3Epc/++zL34htwLMFftzQV3/7y717j7XKI/Dw+4+4eYwPcEwN13AK+Y2cKQdA7RllwT+r4QDb0tNbPi8P/TOURzqePivsS1F9yE4O69ZvZJoo1Rk8DN7v5MzM3KpTcBfwGsM7MnQtoXgK8RvYLjUqL/wS4GcPcD7QU43umewKeAH4Z/rG0CPkz0j+QJe1/cfbWZ3QE8TvQ91xK9HqaUcXBftBWPiIjEQkNwIiISCwUgERGJhQKQiIjEQgFIRERioQAkIiKxUAASGWPM7NNmVhx3O0QOl5Zhi4wxYReFxe7eEHdbRA6HHkQVGcXChpw/JtpSJQn8hGhPsAfNrMHdzzKz84B/AgqAF4EPu3trCFS3A2eF6j7o7htz/R1EhqMhOJHRbRmwzd1PDu+D+Q+ivb3OCsFnCvBF4Fx3Pw2oBT6bUb7Z3ZcA/xnKiowaCkAio9s64Fwzu9bM3uLuewddX0r0ssPfh+2OVgBzMq7/KOPnG7PdWJFDoSE4kVHM3V8ws0XA24F/NrN7B2UxYJW7f2C4KoY5F4mdekAio5iZzQDa3f0HRC8mOw1oAcpClj8CbzKzo0P+YjM7JqOK92f8fCQ3rRYZGfWAREa3E4F/NbN+oAe4nGgo7W4z2x7mgf4S+JGZFYQyXwReCOcFZraa6B+bw/WSRGKhZdgi45SWa8topyE4ERGJhXpAIiISC/WAREQkFgpAIiISCwUgERGJhQKQiIjEQgFIRERi8f8BUTQRy0A2wRMAAAAASUVORK5CYII=)



## 11) Equilibration: assembling the binary input

```
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
```







