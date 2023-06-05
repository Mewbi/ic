import numpy as np
from scipy import optimize

from optimization import base

class FunctionScipy(base.Function):
    '''
    Find criticial points of a function using Newtons Method

    ...

    Attributes
    ----------
    function: function
        Function that will be optimized

    point: float list
        Initial point of otimization


    Methods
    ----------
    converge():
        Try to convert a function using Newtons Method
    '''

    def __init__(self, function, point):
        super().__init__(function, point)

        self.hessian_numerical = []
        for i, _ in enumerate(self.point):
            for j, _ in enumerate(self.point):
                self.hessian_numerical.append(self.__second_derivate(i, j))

    def __second_derivate(self, i, j):
        f = self.derivatives_numerical[i]
        return lambda x: self.__derivative_numerical(f, x, j)
        
    
    def __derivative_numerical(self, f, point, index):
        '''
        Calculate the partial derivative numerically of function
        '''
        old = np.copy(point)
        new = np.copy(point)
        new[index] += base.Function.DEFAULT_GAP_NUMERICAL
        d = ( f(new) - f(old) ) / base.Function.DEFAULT_GAP_NUMERICAL
        return d


    def converge(self, method='CG',
                 max_iterations = base.Function.DEFAULT_MAX_ITERATIONS,
                 tolerance = base.Function.DEFAULT_TOLERANCE):
        def grad(x):
            '''
            x -> np.array

            return: np.array (same size of x0)
            '''
            derivatives = self.derivatives_numerical
            return np.array([d(x) for d in derivatives])

        def hess(x):
            '''
            x -> np.array

            return: np.array
            '''

            hessian = self.hessian_numerical
            h = np.array()
            print(hessian)
            for line in hessian:
                l = np.array()
                for f in line:
                    np.append(l, f(x))
                np.append(l, h)
            return h
                

        print('inside converge: before minimize')
        init_val = self.function(self.point)
        print('hess: {}'.format(hess(point)))

        res = optimize.minimize(self.function, self.point,
                                      method=method, jac=grad, tol=tolerance,
                                      callback=None, options={'maxiter':max_iterations})
        #print('inside converge: after')

        converge = res.success
        point = res.x
        iterations = res.nit
        final_val = self.function(point)

        print(f'Results from scipy: {res}')
        print(f'Message from scipy: {res.message}')

        if not converge:
            raise ValueError("Method didn't converge.")

        return point, iterations, init_val, final_val
