# Scientific Research

This repository will store some stuff about my scientific research to [UFABC].

Some codes and notes will not have a previous context.

## Thematic

The main thematic is Geometric Optimization.

My notes could be found [here](https://energetic-blinker-147.notion.site/IC-5f13b1deac5f4073808ff39dc48a2302)

## Objetive

The main objetive of this research is found a way to receive a generic function
with `n` parameters, an initial point and return a point where the function
converge, this point could be a local maximum, local minimum or a saddle.

Some cases (probably most), we cannot found the exactly point where all partial
derivatives is equal to zero. So this method have a tolerance value close to zero.

## Method

This method is based to secant method to converge the function.

## Example

The following code is and example of a function `f(x,y) = -sin(x)*sin(y)`
converging using numerical and analytical method.

The difference between both method is, in analytical we must specify the
function, an initial point and all partial derivaties of function. In numerical
method, we just need specify the function and initial point.

```python
import math
from optimization import optimization as op

# The function
def fn(p):
    x = p[0]
    y = p[1]
    return -math.sin(x)*math.sin(y)

# Partial derivative to X
def dx_fn(p):
    x = p[0]
    y = p[1]
    return -math.cos(x)*math.sin(y)

# Partial derivative to Y
def dy_fn(p):
    x = p[0]
    y = p[1]
    return -math.cos(y)*math.sin(x)


initial_point = [2.5, 0.8]
func = op.Function(fn, initial_point)

# Converge numerically
try:
    point = func.converge_numerically()
    print("[ Numerically] - Converge point: {} ".format(point))
except Exception as err:
    print(err)


# Converge analytically
func.derivatives = [dx_fn, dy_fn]
try:
    point = func.converge_analytically()
    print("[ Analytically ] - Converge point: {} ".format(point))
except Exception as err:
    print(err)
```

Output

```bash
[ Numerically] - Converge point: [1.5707915451117773, 1.5707913195622745] 
[ Analytically ] - Converge point: [1.5707965450768195, 1.5707963195694736]
```

## Applying in Geometric Optimization
