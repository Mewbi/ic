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

### Example

The following code is and example of a function `f(x,y) = -sin(x)*sin(y)`
converging using numerical .

```python
```

Output

```bash
[ Numerically] - Converge point: [1.5707915451117773, 1.5707913195622745] 
[ Analytically ] - Converge point: [1.5707965450768195, 1.5707963195694736]
```

