#!/usr/bin/env python3

from fh2o_module import li_dawes_guo as ldg
from optimization import cbpd
from optimization import scipy as scp

ldg.init()

#point = [0.9609, 0.9609, 10, 98.940315, 0.0001, 0.0001]
point = [0.726975, 1.025, 1.3536, 103.2989, 118.2898, 68.3691]
optimize_vars = [True, False, True, False, False, False]
func = scp.FunctionScipy(ldg.pes, point, optimize_vars)
#func = cbpd.FunctionCBPD(ldg.pes, point)

result = func.converge_newtown(tolerance=0.00001,
                               max_iterations=100,
                               details=True)

print(result.final_point)
print(result.convergence_steps)
result.plot()
result.save_steps()
