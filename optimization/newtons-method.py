# Simple implemation of newtons method to find root of a function
# Uses only the concept function and your derivate

import functions as fn
import graph


def secant_step(f, old_p, new_p, old_pts, new_pts):
    # Secant
    return old_p - f(old_pts) * ((old_p-new_p) / (f(old_pts) - f(new_pts) ))

def check_converge(functions, pts, tolerance):
    for f in functions:
        if abs(f(pts)) > tolerance:
            return False
    return True

def newtons_method(function, derivates, pts, gap, n = 10, tolerance=1e-5):
    x = [pts[0]]
    y = [pts[1]]
    z = [function(pts)]

    # Define a next point
    old = pts
    new = []
    for p in pts:
        new.append(p+gap)

    for _ in range(n):
        new_tmp = []
        d_values = []
        for p, d in enumerate(derivates):
            d_values.append(d(old))
            v = secant_step(d, old[p], new[p], old, new)
            new_tmp.append(v)

        old = new[:]
        new = new_tmp[:]

        print("Function = {:.4f} - Points {} - Derivates {}".format( function(new) , new, d_values))
        y.append(new[1])
        x.append(new[0])
        z.append(function(new))


        if check_converge(derivates, new, tolerance):
            return new, [x, y, z]


    raise Exception("Newton's Method didn't converge.")


try:
    f = fn.functions[0] # Main function f(x,y)
    d = fn.derivates[0] # Partial derivates df/dx, df/dy
    p = fn.points[0] # Initial point

    print("Starting !!!!!")

    c, pts = newtons_method(f, d, p(), 0.1, 1000000, 1e-3)
        
    print("Converge in {}  - Function Value = {:.4f}".format(c, f(c)))
    graph.plot_graph(f, pts[0], pts[1], pts[2])
except Exception as e:
    print(e)
