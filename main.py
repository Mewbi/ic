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

