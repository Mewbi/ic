#!/usr/bin/env python3

from fh2o_module import li_dawes_guo as ldg
from optimization import cbpd
from optimization import scipy as scp

ldg.init()

#point = [0.9609, 0.9609, 10, 98.940315, 0.0001, 0.0001]
point = [0.9728, 10, 0.9223, 300, 300, 300]
optimize_vars = [True, False, True, False, False, False]
func = scp.FunctionScipy(ldg.pes, point, optimize_vars)
#func = cbpd.FunctionCBPD(ldg.pes, point)

result = func.converge_newtown(tolerance=0.00001,
                               max_iterations=100)

print(result)
