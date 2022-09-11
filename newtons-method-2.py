from math import exp

A = ( -200 , -100 , -170 , 15)
a = ( -1 , -1 , -6.5 , 0.7)
b = (0 , 0 , 11 , 0.6)
c = ( -10 , -10 , -6.5 , 0.7)
x_0 = (1 , 0 , -0.5 , -1)
y_0 = (0 , 0.5 , 1.5 , 1)

def inv_matrix(A):
    det = A[0]*A[3] - A[1]*A[2]
    return (A[3]/det, -A[1]/det, -A[2]/det, A[0]/det)

def mult_matrix(A, v):
    return (A[0]*v[0] + A[2]*v[1], A[1]*v[0] + A[3]*v[1])

def f(i, x, y):
    aux_x = x - x_0[i]
    aux_y = y - y_0[i]
    return a[i]*aux_x*aux_x + b[i]*aux_x*aux_y + c[i]*aux_y*aux_y

def partial_f_x(i, x, y):
    aux_x = x - x_0[i]
    aux_y = y - y_0[i]
    return 2*a[i]*aux_x + b[i]*aux_y

def partial_f_y(i, x, y):
    aux_x = x - x_0[i]
    aux_y = y - y_0[i]
    return b[i]*aux_x + 2*c[i]*aux_y

def grad_V(x, y):
    D_x = D_y = 0
    for i in range(4):
        D_x += A[i]*exp(f(i, x, y))*partial_f_x(i, x, y)
        D_y += A[i]*exp(f(i, x, y))*partial_f_y(i, x, y)
    return ( D_x , D_y )

def hess_V(x, y):
    D_xx = D_xy = D_yy = 0
    for i in range(4):
        aux_x = partial_f_x(i, x, y)
        aux_y = partial_f_y(i, x, y)
        D_xx += A[i]*exp(f(i, x, y))*(aux_x*aux_x + 2*a[i])
        D_yy += A[i]*exp(f(i, x, y))*(aux_y*aux_y + 2*c[i])
        D_xy += A[i]*exp(f(i, x, y))*(aux_x*aux_y + b[i])
    return ( D_xx , D_xy , D_xy , D_yy )

def newton_step(f, fp, x, y):
    (step_x, step_y) = mult_matrix(inv_matrix(fp(x ,y)), f(x, y))
    return (x - step_x, y - step_y)

def newtons_method(f, fp, start, n_iter=10, tolerance=1e-10):
    ( prev_x , prev_y ) = start
    for _ in range(n_iter):
        (new_x, new_y) = newton_step(f, fp, prev_x, prev_y)
        ( grad_x , grad_y ) = f(new_x ,new_y)
        if abs(grad_x) + abs(grad_y) < tolerance:
            return (new_x, new_y)
        (prev_x, prev_y) = (new_x, new_y)
    raise Exception ("Newton's Method didn't converge.")


try:
    c = newtons_method(grad_V, hess_V, (0, 0), 10000000)
    print("Converge in {}".format(c))
except Exception as e:
    print(e)
