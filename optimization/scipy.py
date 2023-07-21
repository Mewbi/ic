import numpy as np
from numpy import linalg
from scipy import optimize

from optimization import base
from optimization import result

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

    def __init__(self, function, point, optimize_vars = []):
        super().__init__(function, point, optimize_vars)

        n = self.dimension
        self.hessian_numerical = []

        optimize = self.optimize_vars
        pos_x = 0
        for i in range(n):
            if optimize[i] is not True:
                continue

            self.hessian_numerical.append([])
            
            for j in range(n):
                if optimize[j] is not True:
                    continue
                self.hessian_numerical[pos_x].append(self.__second_derivate(i, j))
            pos_x += 1
    


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
        optimize = self.optimize_vars
        grad = []
        for i, d in enumerate(derivatives):
            if optimize[i] is not True:
                continue
            grad.append(d(x))
        return np.array(grad)


    def __hess(self, x):
        '''
        x -> np.array

        return: np.array
        '''

        hessian = self.hessian_numerical
        n = len(hessian)
        h = np.empty((n,n))
        
        for i, line in enumerate(hessian):
            for j, f in enumerate(line):
                h[i,j] = f(x)
        return h

    def __points_optimize(self, x):
        p = []
        optimize = self.optimize_vars
        for i, point in enumerate(x):
            if optimize[i] is not True:
                continue
            p.append(point)

        return np.array(p)

    def __points_optimize_to_full_size(self, x, x_opt):
        if len(x) == len(x_opt):
            return x_opt

        n = self.dimension
        optimize = self.optimize_vars
        x_new = np.empty(n)
        for i, opt in enumerate(optimize):
            if opt is True:
                # Pop the first item from array
                x_new[i], x_opt = x_opt[0], x_opt[1:]
                continue
            x_new[i] = x[i]
        return x_new

    def converge_newtown(self,
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
                A Result object with informations about convergence process
        '''

        init_point = self.point
        init_value = self.function(init_point)
        converge = False
        iterations = 0
        p = np.copy(self.point)


        np.set_printoptions(linewidth=100)
        grad = self.__grad(p)
        norm_grad = linalg.norm(grad)
        print("\n Trying to converge: {}".format(p.tolist()))
        for _ in range(max_iterations):
            
            grad = self.__grad(p)
            norm_grad = linalg.norm(grad)
            if norm_grad < tolerance:
                converge = True
                break
            
            hess = self.__hess(p)
            det = linalg.det(hess)
            if det == 0: # Cannot get inv matrix when det is zero
                print("Matriz Singular")
                break

            hess_inv = linalg.inv(hess)
            p_optimize = self.__points_optimize(p)
            p_optimize = p_optimize - hess_inv @ grad
            p = self.__points_optimize_to_full_size(p, p_optimize)
            iterations += 1

        norm_grad = f"{norm_grad:.1e}"
        point = p.tolist()
        final_value = self.function(point)

        if final_value == self.FAIL_CONVERGE_VALUE:
            converge = False

        r = result.Result(
                converge = converge,
                init_point = init_point,
                final_point = point,
                init_value = init_value,
                final_value = final_value,
                iterations = iterations,
                gradient = norm_grad,
            )

        return r


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

        if final_value == self.FAIL_CONVERGE_VALUE:
            converge = False

        result = {
            'converge': converge,
            'point': point,
            'iterations': iterations,
            'init_value': init_value,
            'final_value': final_value
        }

        return result
