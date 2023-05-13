import numpy as np
from fh2o_module import li_dawes_guo as ldg
from optimization import optimization as op

# Defined points
# points = [
#    [0.900, 0.900, 9, 100, 2, 1],
#    [0.9668, 0.9668, 2.2992, 105.07, 70.54, -83.05],
#]


def try_converge(point):

    print("Initial Point: {}".format(point))
    func = op.Function(ldg.pes, point)

    try:
        p = func.converge_numerically(tolerance=0.00001, max_iterations=10000)
        print("Converge: Success")
        print("Converge in: {}".format(p))
    except Exception as err:
        print("Converge: Fail")

def gradual_converge(point, relevant_vars, variation, limit):

    for var_idx, is_relevant in enumerate(relevant_vars):
        if not is_relevant:
            continue

        for l in range(limit):
            p = point[:]
            p[var_idx] -= point[var_idx] * (variation * l)
            try_converge(p)

        for l in range(limit):
            p = point[:]
            p[var_idx] += point[var_idx] * (variation * l)
            try_converge(p)

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



ldg.init()

for geometry in geometries:
    print("\nSpecie: {}".format(geometry['specie']))
    print("Stationary Geometry: {}".format(geometry['stationary']))
    gradual_converge(geometry['stationary'], geometry['relevant_vars'], 0.05, 5)
