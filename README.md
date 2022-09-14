# Scientific Research

This repository will store some stuff about my scientific research.

Some codes and notes will not have a previous context.

## Thematic

The main thematic is Geometric Optimization.

My notes could be found [here](https://energetic-blinker-147.notion.site/IC-5f13b1deac5f4073808ff39dc48a2302)

## Objetive

The main objetive of this research is found a way to receive a generic function
with `n` parameters, an initial point and a variation value and minimize this
function returning the energy of molecule, the position and the gradient value.

Something like

```python
def fn(x1, x2, ..., xn):
...

p = (x1, x2, ..., xn)

h = 0.1

energy, position, grad = minimize(fn, p, h)
```

Some cases (probably most), we cannot found the exactly point when `grad == 0`,
so this function will have a tolerance value.

Return `gradient` value is important to check how much distance we are from the
real minimum.

## Steps

Probably we will use the newthon method in this project.

But to direct our studies, the following steps will be followed:

### Root of Function

To understand better the newthon method, the first step is implement a way to find
roots of a generic function.

> Why start searching root of a function?

Simple, it's the simpler way to understand this method. We can avoid necessity to
calculate the Hessian and others annoying steps.

### Minimize Function

Understading the previous step, now we should found the minimum points of a function.

It's an "improve" of the previous method. But in this case we probably gonna use
the first and second derivate of the function.

### Applying in Molecule

With the previous algorithm in hands, apply in some real cases.

Use real molecule energy functions to find the stable points.

### Apply to a non Analytic Function

This step is a little bit abstract.

But in many cases, we will not have an analytic function of molecule energy.

The algorithm should work in these case, receiving the energy value in some
point, and calculate the minimum without the function.

## Benchmark

Just a simples bech with newtons_method 2 write in `python` and `golang`

```bash
# Python
time python3 ./newtons-method-2.py 
Converge in (0.6234994049308765, 0.02803775852868567)

real 0m0,024s
user 0m0,024s
sys  0m0,000s
```

```bash
# Golang
time ./newtons-method-2 
Converge in [0.6234994049308766 0.028037758528685664] 

real 0m0,003s
user 0m0,000s
sys  0m0,003s
```
