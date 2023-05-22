import json
import numpy as np
from fh2o_module import li_dawes_guo as ldg
from optimization import optimization as op


def try_converge(point):

    print("Initial Point: {}".format(point))
    func = op.Function(ldg.pes, point)

    try:
        p, iterations = func.converge_numerically(tolerance=0.00001, max_iterations=10000)
        print("Converge: Success")
        print("Converge in: {}".format(p))
        return True, p, iterations
    except Exception as err:
        print("Converge: Fail")
        return False, [], 0

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
            success, converge_point, iterations = try_converge(p)
            data_var.append({
                "initial point": "{}".format(p),
                "converge": success,
                "converge point": "{}".format(converge_point),
                "iterations": iterations
            })

        # Try converge in 0%
        success, converge_point, iterations = try_converge(point)
        data_var.append({
            "initial point": "{}".format(point),
            "converge": success,
            "converge point": "{}".format(converge_point),
            "iterations": iterations
        })

        # 5 to 25 percent
        for l in range(limit):
            p = point[:]
            p[var_idx] += point[var_idx] * (variation * (l + 1))
            success, converge_point, iterations = try_converge(p)
            data_var.append({
                "initial point": "{}".format(p),
                "converge": success,
                "converge point": "{}".format(converge_point),
                "iterations": iterations
            })

        data[var_name] = data_var

    return data

geometries = [
    {
        "specie": "F + H2O",
        "stationary": [0.9613, 0.9613, 10, 104.20, 0, 0],
        "relevant_vars": [True, True, False, True, False, False] 
    },
    {
        "specie": "R-vdW (F--H2O)",
        "stationary": [0.9668, 0.9668, 2.2992, 105.07, 70.54, -83.05],
        "relevant_vars": [True, True, True, True, True, True] 
    },
    {
        "specie": "TS",
        "stationary": [0.9705, 1.0389, 1.3128, 103.00, 121.99, 68.72],
        "relevant_vars": [True, True, True, True, True, True] 
    },
    {
        "specie": "P-vdW (HO--HF)",
        "stationary": [0.9738, 1.8037, 0.9330, 111.48, 179.16, 0],
        "relevant_vars": [True, True, True, True, True, True] 
    },
    {
        "specie": "HO+HF",
        "stationary": [0.9730, 10, 0.9202, 0, 0, 0],
        "relevant_vars": [True, False, True, False, False, False] 
    },
]

vars_name = [
    "R HO",
    "R OH'",
    "R H'F",
    "θ HOH'",
    "θ OH'F",
    "φ HOH'F"
]

ldg.init()

result = []

for geometry in geometries:
    print("\nSpecie: {}".format(geometry['specie']))
    print("Stationary Geometry: {}".format(geometry['stationary']))
    data = gradual_converge(geometry, 0.05, 5)
    result.append(data)

with open('result-optimization-1-var.json', 'w') as f:
    json.dump(result, f)

print(result)
