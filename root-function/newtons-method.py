# Simple implemation of newtons method to find root of a function
# Uses only the concept function and your derivate

import functions as fn


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
    for _ in range(n):
        new = secant_step2(f, point, gap)
        print(new, f(new))
        if abs(f(new)) < tolerance:
            return new
        point = new[:]
    raise Exception("Newton's Method didn't converge.")



try:
    c = newtons_method(fn.fn2, fn.pt2(), 0.1, 1000000)
    print("Converge in {}".format(c))
except Exception as e:
    print(e)
