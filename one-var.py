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
            #print(case)
            cases.append(case)

        # Try converge in 0%
        case = {
            "point": point,
            "variation": 0,
            "variable": var_name,
            "specie": specie,
            "relevant_vars": relevant_vars,
        }
        #print(case)
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
            #print(case)
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
    result = f.converge_newtown(
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

def parse_results_metrics(r: result.Results):
    parsed = {}
    for result in r.results:
        var = result.variable
        if var not in parsed:
            parsed[var] = {
                    "conv_success": 0,
                    "conv_fail": 0,
                    "total": 0,
                    "conv_percent": 0,
                }

        parsed["var"]["total"] += 1
        if result.converge:
            parsed["var"]["conv_success"] += 1
        else:
            parsed["var"]["conv_fail"] += 1

    for var, values in parsed:
        values["conv_percent"] = values["conv_success"] / values["total"]

results_cbpd = result.Results()

for geo in geometries:
    print("\n--------------\n")
    print("\tCBPD")
    print("\n--------------\n")
    cases = gen_conv_cases(geo, VARIATION_PERCENT, VARIATION_LIMIT)
    for i, case in enumerate(cases):
        r = try_converge_cbpd(case)
        results_cbpd.add_single_result(r)

    results_cbpd.normalize_final_points(geo["energy"])
results_cbpd.csv('results/one_var_cbpd.csv')


results_newton = result.Results()

for geo in geometries:
    print("\n--------------\n")
    print("\tNewton")
    print("\n--------------\n")
    cases = gen_conv_cases(geo, VARIATION_PERCENT, VARIATION_LIMIT)
    for i, case in enumerate(cases):
        r = try_converge_newton(case)
        results_newton.add_single_result(r)

    results_newton.normalize_final_points(geo["energy"])
results_newton.csv('results/one_var_newton.csv')


# Parse results

