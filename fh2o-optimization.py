import csv
import json
import numpy as np
from fh2o_module import li_dawes_guo as ldg
from optimization import cbpd
from optimization import scipy as scp

def f(x, y):
    return x**2 - (3/2)*x*y + y**2

def try_converge(point):

    print("Initial Point: {}".format(point))
    #func = cbpd.FunctionCBPD(ldg.pes, point)
    #func = scp.FunctionScipy(ldg.pes, point)
    func = scp.FunctionScipy(f, [1.1, 2.2])

    try:
        p, iterations, init_val, final_val = func.converge(method='BFGS',
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

# x0 = geometries[1]['stationary']
# x0[0] *= 1.0
# success, converge_point, iterations, init_val, final_val = try_converge(x0)
# print('\n\n\n\n')

# x0 = converge_point
# x0[0] *= 1.01
# success, converge_point, iterations, init_val, final_val = try_converge(x0)
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
                

