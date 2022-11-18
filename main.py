import functions as fn

import optimization.optimization as op

func = op.Optimize(fn.fn3, [2.5, 0.8])
func.derivates = [fn.dx_fn3, fn.dy_fn3]

try:
    point = func.converge_analytical()
    print("Ponto Calculado Analiticamente: {}".format(point))
    point = func.converge_numerical()
    print("Ponto Calculado Numericamente: {} ".format(point))
except Exception as err:
    print(err)
