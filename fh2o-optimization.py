import numpy as np
import li_dawes_guo as ldg
from optimization import optimization as op
import random

# Defined points
points = [
    #[0, 0, 0, 0, 0, 0],
    [0.900, 0.900, 9, 100, 2, 1],
    [0.800, 0.800, 11, 90, 2, 1],
    [0.700, 0.700, 11, 90, 2, 1],
    #[0.9668, 0.9668, 2.2992, 105.07, 70.54, -83.05],
    #[0.234, 0.234, 0.234, 0.234, 0.234, 0.234],
    #[0.152, 0.464, 0.687, 0.664, 0.753, 0.678],
    #[1.124, 2.234, 0.234, 0.654, 0.616, 1.523],
    #[0.734, 0.842, 0.876, 0.123, 0.535, 0.345],
]

ldg.init()

print("---------------DEFINED POINTS---------------")
for point in points:
    print("\n\nTrying to converge from: {} - Initial energy of f_h2o: {}".format(point, ldg.pes(point)))
    func = op.Function(ldg.pes, point)

    try:
        p = func.converge_numerically(tolerance=0.00001)
        print("Converge in point: {}\n Energy of f_h2o: {}".format(p, ldg.pes(p)))
    except Exception as err:
        print(err)

#print("\n\n---------------RANDOM POINTS---------------")
#for i in range(5):
#    point = [random.uniform(0.1, 5.0) for _ in range(6)]
#    print("\n\nTrying to converge from: {} - Initial energy of f_h2o: {}".format(point, ldg.pes(point)))
#    func = op.Function(ldg.pes, point)

#    try:
#        p = func.converge_numerically(gap = 10.0)
#        print("Converge in point: {}\n Energy of f_h2o: {}".format(p, ldg.pes(p)))
#        with open("converge.txt", "a") as file1:
            # Writing data to a file
#            file1.write("Initial Point: {}\nConverge in: {}\nEnergy F2HO: {}\n\n".format(point, p, ldg.pes(p)))
#    except Exception as err:
#        print(err)
