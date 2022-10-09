import math

# fn
# will define a function

# pt
# an suggestion of an initial point

def fn1(p):
    x = p[0]
    y = p[1]
    return y - (1/12)*(y**3) - (1/4)*(x**2)

def dx_fn1(p):
    x = p[0]
    y = p[1]
    return -(1/2)*x

def dy_fn1(p):
    x = p[0]
    y = p[1]
    return 1 - (1/4)*(y**2)

def pt_fn1():
    x = 1124
    y = -1345
    return [x, y]




def fn2(p):
    x = p[0]
    y = p[1]
    return y**2 - x**2

def dx_fn2(p):
    x = p[0]
    y = p[1]
    return -2*x

def dy_fn2(p):
    x = p[0]
    y = p[1]
    return 2*y

def pt_fn2():
    x = 23
    y = -31
    return [x, y]




def fn3(p):
    x = p[0]
    y = p[1]
    return -math.sin(x)*math.sin(y)

def dx_fn3(p):
    x = p[0]
    y = p[1]
    return -math.cos(x)*math.sin(y)

def dy_fn3(p):
    x = p[0]
    y = p[1]
    return -math.cos(y)*math.sin(x)

def pt_fn3():
    x = -21
    y = 51
    return [x, y]




def fn4(p):
    x = p[0]
    y = p[1]
    return (10*math.cos(x*y))/(1+2*(y**2))

def dx_fn4(p):
    x = p[0]
    y = p[1]
    return (-10*y*math.sin(x*y))/(2*(y**2)+1)

def dy_fn4(p):
    x = p[0]
    y = p[1]
    return (-10*((2*x*(y**2) + x) * math.sin(x*y) + 4*y*math.cos(x*y))) / ((2*(y**2) + 1)**2)

def pt_fn4():
    x = -10
    y = -10
    return [x, y]




def fn5(p):
    x = p[0]
    y = p[1]
    return (3 * math.cos(x+y)) / (1 + (x - y)**2 )

def dx_fn5(p):
    x = p[0]
    y = p[1]
    aux_1 = (-3 * math.sin(x+y)) / ( (x-y)**2 + 1 )
    aux_2 = ( 6*(x-y) * math.cos(x+y) ) / ( ( (x-y)**2 + 1 )**2 )
    return aux_1 - aux_2

def dy_fn5(p):
    x = p[0]
    y = p[1]
    aux_1 = (6*(x-y) * math.cos(y+x)) / (((x-y)**2 + 1)**2)
    aux_2 = (3 * math.sin(y+x)) / ( (x-y)**2 + 1)
    return aux_1 - aux_2

def pt_fn5():
    x = 5
    y = -5
    return [x, y]



functions = [fn1, fn2, fn3, fn4, fn5]
derivates = [[dx_fn1, dy_fn1], [dx_fn2, dy_fn2], [dx_fn3, dy_fn3], [dx_fn4, dy_fn4], [dx_fn5, dy_fn5]]
points = [pt_fn1, pt_fn2, pt_fn3, pt_fn4, pt_fn5]
