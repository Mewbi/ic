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

Will be used the secant method to converge the function.

But to direct our studies, the following steps will be followed:

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

