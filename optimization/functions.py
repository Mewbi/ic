import math

# fn
# will define a function

# pt
# an suggestion of an initial point

def fn1(p):
    return p[0]**2 - p[0]- 1

def pt1():
    return [1]

def fn2(p):
    return 9 - p[0]**2 - p[1]**2

def pt2():
    #return [2, 2]
    return [0.5, 2]
    #return [0.5, 2]

def fn3(p):
    return p[0]**2 + p[1]**2

def pt3():
    return [0.5, 2]

def fn4(p):
    return 2*math.cos(p[0]) - p[0] + (p[0]**2)/10 

def pt4():
    return [-5]


functions = [fn1, fn2]
points = [pt1, pt2]
