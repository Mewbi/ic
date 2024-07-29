import json
from fh2o_module import li_dawes_guo as ldg
from optimization import cbpd
from optimization import result
from optimization import scipy as scp
from geometries import *

vars_name = [
    "R HO",
    "R OH'",
    "R H'F",
    "θ HOH'",
    "θ OH'F",
    "φ HOH'F"
]

VARIATION_PERCENT = 0.05
VARIATION_LIMIT = 5
MAX_ITERATIONS = 50
TOLERANCE = 1e-5
DIVERGE_VAL = 115.302739153
ldg.init()

def gen_conv_cases(geometry, variation, limit):
    cases = []

    point = geometry["stationary"]
    relevant_vars = geometry["relevant_vars"]
    specie = geometry["specie"]

    for var_idx, is_relevant in enumerate(relevant_vars):
        if not is_relevant:
            continue

        var_name = vars_name[var_idx]

        # -25 to -5 percent
        for l in range(limit):
            p = point[:]
            p[var_idx] -= point[var_idx] * (variation * (limit - l))
            v = int( (variation * (limit - l)) * 100)
            case = {
                "point": p,
                "variation": v,
                "variable": var_name,
                "specie": specie,
                "relevant_vars": relevant_vars,
            }
            cases.append(case)

        # Try converge in 0%
        case = {
            "point": point,
            "variation": 0,
            "variable": var_name,
            "specie": specie,
            "relevant_vars": relevant_vars,
        }
        cases.append(case)

        # 5 to 25 percent
        for l in range(limit):
            p = point[:]
            p[var_idx] += point[var_idx] * (variation * (l + 1))
            v = int( (variation * (l + 1)) * 100)
            case = {
                "point": p,
                "variation": v,
                "variable": var_name,
                "specie": specie,
                "relevant_vars": relevant_vars,
            }
            cases.append(case)

    return cases[:]


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

results_cbpd = result.Results()

for geo in geometries:
    partial_results = result.Results()
    print("\n--------------\n")
    print("\tCBPD")
    print("\n--------------\n")
    cases = gen_conv_cases(geo, VARIATION_PERCENT, VARIATION_LIMIT)
    for i, case in enumerate(cases):
        r = try_converge_cbpd(case)
        partial_results.add_single_result(r)

    for p in partial_results.results:
        print(p.specie, p.variable, p.variation, p.init_value)

    # partial_results.normalize_final_points(geo["energy"])
    results_cbpd.add_multiple_results(partial_results.results)
results_cbpd.csv('results/one_var_cbpd.csv')


results_newton = result.Results()

for geo in geometries:
    partial_results = result.Results()
    print("\n--------------\n")
    print("\tNewton")
    print("\n--------------\n")
    cases = gen_conv_cases(geo, VARIATION_PERCENT, VARIATION_LIMIT)
    for i, case in enumerate(cases):
        r = try_converge_newton(case)
        partial_results.add_single_result(r)

    # partial_results.normalize_final_points(geo["energy"])
    results_newton.add_multiple_results(partial_results.results)
results_newton.csv('results/one_var_newton.csv')

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

print("\n\n--------------\n")
print("\tLatex Data")
print("\n--------------\n")

order_conv = ["HO+HF", "R-vdW (F--H2O)", "TS", "P-vdW (HO--HF)", "F + H2O"]
print("\n--------------\n")
print("\tCBPD - Convergence Tax")
print("\n--------------\n")
for i, conv in enumerate(order_conv):
    data = parsed_cbpd["convergence"][conv]
    print("({}, {:.2f})".format(i+1, data["conv_percent"]))

print("\n--------------\n")
print("\tNewton - Convergence Tax")
print("\n--------------\n")
for i, conv in enumerate(order_conv):
    data = parsed_newton["convergence"][conv]
    print("({}, {:.2f})".format(i+1, data["conv_percent"]))

print("\n--------------\n")
print("\tCBPD - Iterations")
print("\n--------------\n")
for it, data in parsed_cbpd["iterations"].items():
    print("({}, {:.2f})".format(it, data["percent"]))

print("\n--------------\n")
print("\tNewton - Iterations")
print("\n--------------\n")
for it, data in parsed_newton["iterations"].items():
    print("({}, {:.2f})".format(it, data["percent"]))
