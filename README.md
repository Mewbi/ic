# 🧪 Scientific Research

This repository will store some stuff about my scientific research to [UFABC](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwiiu-eFr4CFAxXrHLkGHSqNCxcQFnoECAgQAQ&url=https%3A%2F%2Fwww.ufabc.edu.br%2F&usg=AOvVaw1AKiSt62gZiblUScBnWY7c&opi=89978449).

Some codes and notes will not have a previous context.

## ⚗️ Thematic

This research is related to applied math about Optimization of Molecular Geometry. 

The main thematic is Geometric Optimization.

Some previous research and general notes could be found in my [Notion](https://energetic-blinker-147.notion.site/IC-5f13b1deac5f4073808ff39dc48a2302).

## 🎯 Objetive

The main objetive of this research is found a way to receive a generic function
with `n` parameters, an initial point and return a point where the function
converge, this point could be a local maximum, local minimum or a saddle.

To define a convergence, we check if `norm` of actual point (by default is Euclidian Norm)
is equal os smaller than tolerance value (by default is 10e-5).

This method was called CBPD (Convergence Based in Partial Derivatives) and is
bases in Newton's Method and Secant Method.

## 📝 [Notes](notes/)

Notes are formal research documents that explain the project step by step.

They are written using LaTeX and delve into the project's context, explain the algorithm logic, and present a performance comparison between CPBD and Newton's method to converge the SEP function of the F + H2O reaction.

The PDF file could be read [here](notes/00_main.pdf)

## 🔬 [CNMAC](cnmac/)

[CNMAC](https://www.cnmac.org.br/novo/index.php/CNMAC) is a congress about applied math and computation.

The content in this directory is an abstract to submit to 2024 congress.

The abstract PDF could be read [here](cnmac/resumo.pdf)

## 📊 Optimization

This module contains the code related to optmization:

- `base.py`: base code common to every optimization method
- `cbpd.py`: contains the logic of CBPD method
- `scipy.py`: contains some optimization methods from SciPy and my implementation of Newton's method
- `result.py`: contains result data of convergence, which is used by all methods

### Example

The following code is an example of a convergence process of SEP used in this project

```python
from fh2o_module import li_dawes_guo as ldg
from optimization import cbpd

ldg.init() # Required to SEP function works

point = [0.9901, 1.8261, 1.0003, 112.1962, 170.6405, -8.7494] # Initial Point
func = cbpd.FunctionCBPD(ldg.pes, point) # Create an object of optimization

result = func.converge_numerically(tolerance=0.00001,
                                max_iterations=100,
                                norm="euclidian")

print(result)
```

Output

```bash
Converge: True
Iterations: 24
Init Point: [0.9901, 1.8261, 1.0003, 112.1962, 170.6405, -8.7494]
Final Point: [0.9728234257540236, 1.7700969608929964, 0.9354300701098185, 108.64886591218642, 173.61687633352184, 0.10001533816143876]
Init Value: -19.63909180874393
Final Value: -22.337468737466477
```

