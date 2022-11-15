import functions as fn

import optimization.optimization as op

func = op.Optimize(fn.fn3, [2.5, 0.8])
func.derivates = [fn.dx_fn3, fn.dy_fn3]

try:
    point = func.converge_analytical()
    print(point)
except:
    print("Deu ruim")

