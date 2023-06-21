#!/usr/bin/env python3

import csv
import json
import numpy as np
from fh2o_module import li_dawes_guo as ldg
from optimization import cbpd
from optimization import scipy as scp
from geometries import *

def f(x):
    return x[0]**3 - (3/2)*x[0]*x[1] + x[1]**2

def try_converge(point):

    print("Initial Point: {}".format(point))
    #func = cbpd.FunctionCBPD(ldg.pes, point)
    func = scp.FunctionScipy(ldg.pes, point)
    #func = scp.FunctionScipy(f, point)

    try:
        p, iterations, init_val, final_val = func.converge(method='Newton-CG',
                                      tolerance=0.00001,
                                      max_iterations=1000)
        #p, iterations, init_val, final_val = func.converge_numerically(tolerance=0.00001, max_iterations=10000)
        print("Converge: Success")
        print("Converge in: {}".format(p))
        return True, p, iterations, init_val, final_val
    except Exception as err:
        print("Converge: Fail: {}".format(err))
        return False, [], 0, 0, 0

def gradual_converge(geometry, variation, limit):

    point = geometry["stationary"]
    relevant_vars = geometry["relevant_vars"]

    data = {
            "specie": geometry["specie"]
            }
    for var_idx, is_relevant in enumerate(relevant_vars):
        if not is_relevant:
            continue

        data_var = []
        var_name = vars_name[var_idx]

        # -25 to -5 percent
        for l in range(limit):
            p = point[:]
            p[var_idx] -= point[var_idx] * (variation * (limit - l))
            success, converge_point, iterations, init_val, final_val = try_converge(p)
            data_var.append({
                "initial point": "{}".format(p),
                "converge": success,
                "converge point": "{}".format(converge_point),
                "iterations": iterations,
                "init value": init_val,
                "final value": final_val
            })

        # Try converge in 0%
        success, converge_point, iterations, init_val, final_val = try_converge(point)
        data_var.append({
            "initial point": "{}".format(point),
            "converge": success,
            "converge point": "{}".format(converge_point),
            "iterations": iterations,
            "init value": init_val,
            "final value": final_val
        })

        # 5 to 25 percent
        for l in range(limit):
            p = point[:]
            p[var_idx] += point[var_idx] * (variation * (l + 1))
            success, converge_point, iterations, init_val, final_val = try_converge(p)
            data_var.append({
                "initial point": "{}".format(p),
                "converge": success,
                "converge point": "{}".format(converge_point),
                "iterations": iterations,
                "init value": init_val,
                "final value": final_val
            })

        data[var_name] = data_var

    return data

vars_name = [
    "R HO",
    "R OH'",
    "R H'F",
    "θ HOH'",
    "θ OH'F",
    "φ HOH'F"
]

ldg.init()

# x0 = geometries[1]['stationary']
# x0[0] *= 1.0
# success, converge_point, iterations, init_val, final_val = try_converge(x0)
# print('\n\n\n\n')

# x0 = converge_point
# x0[0] *= 1.01

# success, converge_point, iterations, init_val, final_val = try_converge([1.1, 2.2])
# print(success, converge_point, iterations, init_val, final_val)

# exit()

result = []

for geometry in geometries:
    print("\nSpecie: {}".format(geometry['specie']))
    print("Stationary Geometry: {}".format(geometry['stationary']))
    data = gradual_converge(geometry, 0.05, 5)
    result.append(data)

#with open('result-optimization-1-var.json', 'w') as f:
#    json.dump(result, f)


# Save Result in CSV
header = ["specie", 'variable', "variation", "converge", "iterations", "init_value", "final_value", "init_point", "final_point"]
variations = [-25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25]
with open('result-optimization-1-var.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(header)

    for d in result:
        specie = d["specie"]

        for var in vars_name:

            if not var in d:
                continue
        
            for i, v in enumerate(d[var]):
                variation = variations[i]
                converge = v["converge"]
                iterations = v["iterations"]
                init_val = v["init value"]
                final_val = v["final value"]
                init = v["initial point"]
                final = v["converge point"]

                row = [specie, var, variation, converge, iterations, init_val, final_val, init, final]
                writer.writerow(row)
                

