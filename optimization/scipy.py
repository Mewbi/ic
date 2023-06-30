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

        n = self.dimension
        self.hessian_numerical = []

        for i in range(n):
            self.hessian_numerical.append([])
            for j in range(n):
                self.hessian_numerical[i].append(self.__second_derivate(i, j))

    def __second_derivate(self, i, j):
        f = self.derivatives_numerical[i]
        return lambda x: self.__derivative_numerical(f, x, j)


    def __derivative_numerical(self, f, point, index):
        '''
        Calculate the partial derivative numerically of function
        '''
        old = np.copy(point)
        new = np.copy(point)
        old[index] -= base.Function.DEFAULT_GAP_NUMERICAL
        new[index] += base.Function.DEFAULT_GAP_NUMERICAL
        d = ( f(new) - f(old) ) / (2 * base.Function.DEFAULT_GAP_NUMERICAL)
        return d
    
    def __grad(self, x):
        '''
        x -> np.array

        return: np.array (same size of x0)
        '''
        derivatives = self.derivatives_numerical
        return np.array([d(x) for d in derivatives])


    def __hess(self, x):
        '''
        x -> np.array

        return: np.array
        '''

        hessian = self.hessian_numerical
        n = self.dimension
        h = np.empty((n,n))
        
        for i, line in enumerate(hessian):
            for j, f in enumerate(line):
                h[i,j] = f(x)
        return h

    def converge(self, method='CG',
                 max_iterations = base.Function.DEFAULT_MAX_ITERATIONS,
                 tolerance = base.Function.DEFAULT_TOLERANCE):
        '''
        Finds the closest convergence point based on the initial point using scipy methods 

            Parameters
            ----------
                max_iterations: int, optional
                    Max iterations to try converge the function
                tolerance float, optional
                    How close gradient should be to converge

            Returns
            ----------
                A dict mapping the result of convergence processs

                {
                    'converge': bool, # True if function has converged and False if not

                    'point': float list, # Point where function has converged or stopped

                    'iterations': int, # Number of iterations

                    'init_value': float, # Init value of function

                    'final_value': float # Final value of function
                }
        '''

        init_value = self.function(self.point)

        res = optimize.minimize(self.function, self.point,
                                      method=method, jac=self.__grad, tol=tolerance, hess=self.__hess,
                                      callback=None, options={'maxiter':max_iterations})

        #print('inside converge: after')

        converge = res.success
        point = res.x.tolist()
        iterations = res.nit
        final_value = self.function(point)

        result = {
            'converge': converge,
            'point': point,
            'iterations': iterations,
            'init_value': init_value,
            'final_value': final_value
        }

        return result
