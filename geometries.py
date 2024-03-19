'''
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
'''

geometries = [
    {
        "specie": "F + H2O",
        "stationary": [0.9609, 0.9609, 20, 104.1477, 300, 300],
        "relevant_vars": [True, True, False, True, False, False],
        "energy": -0.00207
    },
    {
        "specie": "R-vdW (F--H2O)",
        "stationary": [0.9641, 0.9641, 2.2981, 104.3416, 67.1171, -87.5402],
        "relevant_vars": [True, True, True, True, True, True],
        "energy": -3.37121
    },
    {
        "specie": "TS",
        "stationary": [0.9693, 1.0250, 1.3536, 103.2989, 118.2898, 68.3691],
        "relevant_vars": [True, True, True, True, True, True],
        "energy": 1.49907
    },
    {
        "specie": "P-vdW (HO--HF)",
        "stationary": [0.9728, 1.7700, 0.9354, 108.6504, 173.6165, -0.0737],
        "relevant_vars": [True, True, True, True, True, True],
        "energy": -22.28683
    },
    {
        "specie": "HO+HF",
        "stationary": [0.9728, 20, 0.9223, 300, 300, 300],
        "relevant_vars": [True, False, True, False, False, False],
        "energy": -16.15134
    },
    # Ignoring this to avoid incorrct metrics
    # {
    #     "specie": "P'-vdW (HO--FH)",
    #     "stationary": [0.9735, 3.9974, 0.9220, 17.51, 46.76, 0.0001],
    #     "relevant_vars": [True, True, True, True, True, True],
    #     "energy": -20.29986
    # },
]
