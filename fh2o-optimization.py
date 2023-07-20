#!/usr/bin/env python3

import numpy as np
from fh2o_module import li_dawes_guo as ldg
from optimization import base
from optimization import cbpd
from optimization import scipy as scp
from geometries import *

def f(x):
    return x[0]**3 - (3/2)*x[0]*x[1] + x[1]**2

def try_converge(point, relevant_vars, geometry, var_name, variation):

    print("Initial Point: {}".format(point))
    #func = cbpd.FunctionCBPD(ldg.pes, point)
    func = scp.FunctionScipy(ldg.pes, point, relevant_vars)
    #func = scp.FunctionScipy(f, point)

    #result = func.converge(method='Newton-CG',
    #                              tolerance=0.00001,
    #                              max_iterations=1000)
    result = func.converge_newtown(tolerance=0.00001,
                                  max_iterations=50)
    #p, iterations, init_val, final_val = func.converge_numerically(tolerance=0.00001, max_iterations=10000)

    if result.converge:
        print("Converge: Success - {}".format(result.final_point))
    else:
        print("Converge: Fail")

    result.specie = geometry["specie"]
    result.variable = var_name
    result.variation = variation
    return result

def gradual_converge(geometry, variation, limit):

    results = base.Results()

    point = geometry["stationary"]
    relevant_vars = geometry["relevant_vars"]

    for var_idx, is_relevant in enumerate(relevant_vars):
        if not is_relevant:
            continue

        var_name = vars_name[var_idx]

        # -25 to -5 percent
        for l in range(limit):
            p = point[:]
            p[var_idx] -= point[var_idx] * (variation * (limit - l))
            v = int( (variation * (limit - l)) * 100)
            r = try_converge(p, relevant_vars, geometry, var_name, v)
            results.add_single_result(r)

        # Try converge in 0%
        r = try_converge(point, relevant_vars, geometry, var_name, 0)
        results.add_single_result(r)

        # 5 to 25 percent
        for l in range(limit):
            p = point[:]
            p[var_idx] += point[var_idx] * (variation * (l + 1))
            v = int( (variation * (l + 1)) * 100)
            r = try_converge(p, relevant_vars, geometry, var_name, v)
            results.add_single_result(r)

    return results

vars_name = [
    "R HO",
    "R OH'",
    "R H'F",
    "θ HOH'",
    "θ OH'F",
    "φ HOH'F"
]

ldg.init()


results = base.Results()

for geometry in geometries:
    r = gradual_converge(geometry, 0.05, 5)
    r.normalize_final_points(geometry["energy"])
    results.add_multiple_results(r.results)

results.csv('results/result-optimization-1-var.csv')
