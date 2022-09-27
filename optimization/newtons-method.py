# Simple implemation of newtons method to find root of a function
# Uses only the concept function and your derivate

import functions as fn
import graph

def secant_step(f, pts, gap):
    pts_tmp = pts[:]
    pts_new = []
    for i, p in enumerate(pts):
        pts_tmp[i] = p+gap
        print(pts, pts_tmp, f(pts), f(pts_tmp))

        v = p - f(pts)*((sum(pts)-sum(pts_tmp))/(f(pts)-f(pts_tmp)))
        pts_new.append(v)
        pts_tmp = pts[:]

    return pts_new
    # Secant 
#    return x1 - f(x1)*((x1-x2)/(f(x1)-f(x2)))

def secant_step2(f, pts, gap):
    pts_tmp = pts[:]
    pts_tmp = [x+gap for x in pts_tmp]
    v = f(pts)*((sum(pts)-sum(pts_tmp))/(f(pts)-f(pts_tmp)))
    pts_new = []
    for _, p in enumerate(pts):        
        pts_new.append(p-v)

    return pts_new
    # Secant 
#    return x1 - f(x1)*((x1-x2)/(f(x1)-f(x2)))

def newtons_method(f, point, gap, n = 10, tolerance=1e-10):
    x = [point[0]]
    y = [point[1]]
    z = [f(point)]
    for _ in range(n):
        new = secant_step2(f, point, gap)
        print(new, f(new))
        x.append(new[0])
        y.append(new[1])
        z.append(f(new))
        if abs(f(new)) < tolerance:
            return new, [x, y, z]
        point = new[:]
    raise Exception("Newton's Method didn't converge.")


try:
    c, pts = newtons_method(fn.fn2, fn.pt2(), 0.1, 1000000)
    graph.plot_graph(fn.fn2, pts[0], pts[1], pts[2])
    print("Converge in {}".format(c))
except Exception as e:
    print(e)
