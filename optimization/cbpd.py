# CBPD
# Convergence Based in Partials Derivatives

import numpy as np

from optimization import base
from optimization import result

class FunctionCBPD(base.Function):
    '''
    Find criticial points of a function using numerical or analytic derivative

    ...

    Attributes
    ----------
    function: function
        Function that will be optimized

    point: float list
        Initial point of otimization

    derivatives: function list
        Partial derivatives of the main function


    Methods
    ----------
    converge_analytics():
        Try to convert a function using his analytics derivatives
    
    converge_numerical():
        Try to convert a function using his numerical derivatives
    '''

    def __check_converge(self, pts, derivatives, tolerance):
        '''
        Verify if every partial derivative in point is lower than tolerance value
        '''
        for derivative in derivatives:
            if abs(derivative(pts)) > tolerance:
                return False
        return True

    def __converge_step(self, f, old_point, old_points, new_points, new_point):
        '''
        Calculate the next point based in secant method
        '''
        try:
            v = old_point - f(old_points) * (( (old_point-new_point) / ( f(old_points) - f(new_points) ) ))
        except ZeroDivisionError:
            return old_point
        return v

    def __converge_method(self, gap, n, tolerance, numerical, norm):
        derivatives = self.derivatives_numerical if numerical else self.derivatives
        old = self._point
        for iteration in range(n):
            new = []
            for i, derivative in enumerate(derivatives):
                old_point = old[i]
                new_points = np.copy(old)
                new_points[i] += gap
                new.append(self.__converge_step(derivative, old_point, old, new_points, new_points[i]))
            
            old = np.copy(new)
            value = norm(derivatives, new)
            if value <= tolerance:
                return True, np.copy(new), iteration + 1

        return False, np.copy(old), n

    def converge_analytically(self, 
                              gap = base.Function.DEFAULT_GAP, 
                              max_iterations = base.Function.DEFAULT_MAX_ITERATIONS, 
                              tolerance = base.Function.DEFAULT_TOLERANCE,
                              norm = "euclidian"):
        '''
        Finds the closest convergence point based on the initial point using partial derivatives

            Parameters
            ----------
                gap: float, optional
                    Used to estimate the second derivative in actual point (default is 0.01)
                max_iterations: int, optional
                    Max iterations to try converge the function
                tolerance float, optional
                    How close gradient should be to converge
                norm
                    Norm that will be used to check convergence


            Returns
            ----------
                A Result object with informations about convergence process
        '''

        if not self.derivatives:
            raise ValueError("Partial derivatives not defined.")


        if norm == "euclidian":
            norm_func = self.euclidian_norm
        elif norm == "maximum":
            norm_func = self.maximum_norm
        else:
            raise ValueError("Invalid norm: {}".format(norm))


        init_point = self.point
        init_value = self.function(self.point)
        converge, point, iterations = self.__converge_method(gap, max_iterations, tolerance, False, norm_func)
        final_value = self.function(point)

        r = result.Result(
                converge          = converge,
                init_point        = init_point,
                final_point       = point,
                init_value        = init_value,
                final_value       = final_value,
                iterations        = iterations,
            )

        return r


    def converge_numerically(self, 
                            gap = base.Function.DEFAULT_GAP,
                            max_iterations = base.Function.DEFAULT_MAX_ITERATIONS,
                            tolerance = base.Function.DEFAULT_TOLERANCE,
                            norm = "euclidian"):
        '''
        Finds the closest convergence point based on the initial point using partial derivatives calculated numeracly

            Parameters
            ----------
                gap: float, optional
                    Used to estimate the second derivative in actual point (default is 0.01)
                max_iterations: int, optional
                    Max iterations to try converge the function
                tolerance float, optional
                    How close gradient should be to converge
                norm
                    Norm that will be used to check convergence

            Returns
            ----------
                A Result object with informations about convergence process
        '''

        if norm == "euclidian":
            norm_func = self.euclidian_norm
        elif norm == "maximum":
            norm_func = self.maximum_norm
        else:
            raise ValueError("Invalid norm: {}".format(norm))

        init_point = self.point
        init_value = self.function(self.point)
        converge, point, iterations = self.__converge_method(gap, max_iterations, tolerance, True, norm_func)
        point = point.tolist()
        final_value = self.function(point)

        r = result.Result(
                converge          = converge,
                init_point        = init_point,
                final_point       = point,
                init_value        = init_value,
                final_value       = final_value,
                iterations        = iterations,
            )

        return r
