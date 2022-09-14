# Simple implemation of newtons method to find root of a function
# Uses only the concept function and your derivate

import functions as fn

def f(x): # Example function
    return x**2 - x - 1

def fp(x): # f'(x)
    return 2*x - 1

def newthon_step(f, fp, x):
    return x - f(x)/fp(x)

def newthon_step2(f, x1, x2):
    # Secant method
    return x1 - f(x1)*((x1-x2)/(f(x1)-f(x2)))

def newtons_method(f, fp, start, n = 10, tolerance=1e-10):
    prev = start
    for _ in range(n):
        #new = newthon_step(f, fp, prev)
        new = newthon_step2(f, prev, prev+1)
        if abs(f(new)) < tolerance:
            return new
        prev = new
    raise Exception("Newton's Method didn't converge.")

try:
    c = newtons_method(fn.fn1, fp, fn.pt1(), 1000)
    print("Converge in {}".format(c))
except Exception as e:
    print(e)
