#!/usr/bin/env python3

from fh2o_module import li_dawes_guo as ldg
from optimization import cbpd
from optimization import scipy as scp

ldg.init()

def f(x):
    return x[0] ** 5 + 2*x[0]**4 - 48*x[0]**3 + 126*x[0]**2 - 81 * x[0]

#point = [0.9609, 0.9609, 10, 98.940315, 0.0001, 0.0001]
# point = [0.9728, 1.77, 0.9354, 86.92032, 173.6165, -0.0737]
point = [0.9901, 1.8261, 1.0003, 112.1962, 170.6405, -8.7494]
# point = [2.0]
# optimize_vars = [True, True, True, True, True, True]
# func = scp.FunctionScipy(f, point, [True])
func = cbpd.FunctionCBPD(ldg.pes, point)

result = func.converge_numerically(
                                tolerance=10e-5,
                                max_iterations=50,
                                norm="euclidian"
                            )
# result = func.converge_newton(tolerance = 0.000001,
#                               max_iterations= 100,
#                               details=True)

print(result)
# print(result.convergence_steps)
# result.plot()
# result.save_steps()
