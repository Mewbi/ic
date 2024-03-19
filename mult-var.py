import json
import random as rd
from fh2o_module import li_dawes_guo as ldg
from optimization import cbpd
from optimization import result
from optimization import scipy as scp
from geometries import *

ANG_RANGE = 10
BOND_RANGE = 0.3
MAX_ITERATIONS = 50
TOLERANCE = 1e-5
DIVERGE_VAL = 115.302739153
CASES_PER_GEOMETRY = 100

ldg.init()

def gen_random_conv_case(geometry):
    point = geometry["stationary"][:] # Copy array
    relevant_vars = geometry["relevant_vars"]
    specie = geometry["specie"]
    
    for var_idx, is_relevant in enumerate(relevant_vars):
        if not is_relevant:
            continue

        var_range = BOND_RANGE
        if var_idx > 2: # Means its an angle variation
            var_range = ANG_RANGE

        p = rd.uniform(-var_range, var_range)
        point[var_idx] += p

    return {
        "point": point,
        "variation": -1,
        "variable": "random",
        "specie": specie,
        "relevant_vars": relevant_vars,
    }

def try_converge_cbpd(case):
    p = case["point"]
    v = case["variation"]
    var = case["variable"]
    specie = case["specie"]

    print("[CBPD] - Trying converge init point: {} - ".format(p), end='')
    f = cbpd.FunctionCBPD(ldg.pes, p)
    result = f.converge_numerically(
            tolerance = TOLERANCE,
            max_iterations = MAX_ITERATIONS,
        )
   
    if result.final_value == DIVERGE_VAL:
        result.converge = False

    if result.converge:
        print("Success")
    else:
        print("Fail")

    result.specie = specie
    result.variable = var
    result.variation = v
    return result

def try_converge_newton(case):
    p = case["point"]
    v = case["variation"]
    var = case["variable"]
    specie = case["specie"]
    relevant_vars = case["relevant_vars"]

    print("[Newton] - Trying converge init point: {} - ".format(p), end='')
    f = scp.FunctionScipy(ldg.pes, p, relevant_vars)
    result = f.converge_newton(
            tolerance = TOLERANCE,
            max_iterations = MAX_ITERATIONS,
        )
    
    if result.final_value == DIVERGE_VAL:
        result.converge = False
    
    if result.converge:
        print("Success")
    else:
        print("Fail")

    result.specie = specie
    result.variable = var
    result.variation = v
    return result

# Gen all cases
all_cases = []
for geo in geometries:
    specie = geo['specie']
    energy = geo['energy']
    cases = []
    for i in range(CASES_PER_GEOMETRY):
        c = gen_random_conv_case(geo)
        cases.append(c)

    all_cases.append({
            "cases": cases,
            "energy": energy,
            "specie": specie
        })


# CBPD
results_cbpd = result.Results()
for data in all_cases:
    partial_results = result.Results()
    print("\n--------------\n")
    print("\tCBPD \t" + data["specie"])
    print("\n--------------\n")
    for c in data["cases"]:
        r = try_converge_cbpd(c)
        partial_results.add_single_result(r)

    partial_results.normalize_final_points(data["energy"])
    results_cbpd.add_multiple_results(partial_results.results)

results_cbpd.csv('results/mult_var_cbpd.csv')

# Newton
results_newton = result.Results()
for data in all_cases:
    partial_results = result.Results()
    print("\n--------------\n")
    print("\tNewton \t" + data["specie"])
    print("\n--------------\n")
    for c in data["cases"]:
        r = try_converge_newton(c)
        partial_results.add_single_result(r)

    partial_results.normalize_final_points(data["energy"])
    results_newton.add_multiple_results(partial_results.results)

results_newton.csv('results/mult_var_newton.csv')

# Parse results
parsed_cbpd = results_cbpd.get_results_metrics()
parsed_newton = results_newton.get_results_metrics()

print("\n--------------\n")
print("\tCBPD")
print("\n--------------\n")
print(json.dumps(parsed_cbpd, indent = 2))

print("\n--------------\n")
print("\tNewton")
print("\n--------------\n")
print(json.dumps(parsed_newton, indent = 2))
