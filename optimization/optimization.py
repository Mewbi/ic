import numpy as np

class Optimize:
    '''
    Find criticial points of a function using numerical or analytic derivate

    ...

    Attributes
    ----------
    function: function
        Function that will be optimized

    point: float list
        Initial point of otimization

    derivatives: function list
        Partial derivates of the main function


    Methods
    ----------
    converge_analytics():
        Try to convert a function using his analytics derivates
    
    converge_numerical():
        Try to convert a function using his numerical derivates
    '''

    DEFAULT_GAP_NUMERICAL = 1e-5
    DEFAULT_GAP = 0.1
    DEFAULT_MAX_ITERATIONS = 100000
    DEFAULT_TOLERANCE = 1e-5

    def __lambda_derivative(self, i):
        return lambda x: self.__derivative_numericly(x, i)

    def __init__(self, function, point):
        self.function = function
        self.point = point

        self.derivatives_numericly = []
        for i, _ in enumerate(point):
            self.derivatives_numericly.append(self.__lambda_derivative(i))

    @property
    def function(self):
        return self._function

    @function.setter
    def function(self, function):
        if not callable(function):
            raise ValueError("It's not a python callable function.")
        self._function = function

    @property
    def point(self):
        return self._point

    @point.setter
    def point(self, point):
        try:
            iter(point)
            self._point = point
        except:
            raise ValueError("Point is not iterable")

    @property
    def derivates(self):
        return self._derivates

    @derivates.setter
    def derivates(self, derivates):
        self._derivates = []
        if len(derivates) != len(self._point):
            raise ValueError("Derivatives has not the same dimesion of points")
        try:
            iter(derivates)
        except:
            raise ValueError("Derivates is not iterable")
        for derivate in derivates:
            if not callable(derivate):
                raise ValueError("Derivate {} is not a python callable function".format(derivate))
            self._derivates.append(derivate)

    def __check_converge(self, pts, tolerance):
        '''
        Verify if every partial derivative in point is lower than tolerance value
        '''
        for derivate in self.derivates:
            if abs(derivate(pts)) > tolerance:
                return False
        return True

    def __converge_step(self, f, old_point, old_points, new_points, new_point):
        v = old_point - f(old_points) * (( (old_point-new_point) / ( f(old_points) - f(new_points) ) ))
        return v

    def __converge_method(self, gap, n, tolerance, numerical):
        derivatives = self.derivatives_numericly if numerical else self.derivates
        old = self._point
        for _ in range(n):
            new = []
            for i, derivate in enumerate(derivatives):
                old_point = old[i]
                #old_points = np.array(old[:])
                new_points = old[:]
                new_points[i] += gap
                new.append(self.__converge_step(derivate, old_point, old, new_points, new_points[i]))
                
            old = new[:]
            if self.__check_converge(new, tolerance):
                return True, new[:]

        return False, []

    def converge_analytical(self, gap = DEFAULT_GAP, max_iterations = DEFAULT_MAX_ITERATIONS, tolerance = DEFAULT_TOLERANCE):
        '''
        Finds the closest convergence point based on the initial point using partial derivates

            Parameters
            ----------
                gap: float, optional
                    Used to estimate the second derivate in actual point (default is 0.01)
                max_iterations: int, optional
                    Max iterations to try converge the function
                tolerance float, optional
                    How close gradient should be to converge

            Returns
            ----------
                float list
                    Point where function has converged
        '''

        if not self.derivates:
            raise ValueError("Partial derivates not defined.")

        converge, point = self.__converge_method(gap, max_iterations, tolerance, False)

        if not converge:
            raise ValueError("Method didn't converge.")

        return point


    def __derivative_numericly(self, point, index):
        old = point[:]
        new = point[:]
        new[index] += self.DEFAULT_GAP_NUMERICAL
        f = self._function
        d = ( f(new) - f(old) ) / self.DEFAULT_GAP_NUMERICAL
        return d

    def converge_numerical(self, gap = DEFAULT_GAP, max_iterations = DEFAULT_MAX_ITERATIONS, tolerance = DEFAULT_TOLERANCE):
        '''
        Finds the closest convergence point based on the initial point using partial derivatives calculated numeracly

            Parameters
            ----------
                gap: float, optional
                    Used to estimate the second derivate in actual point (default is 0.01)
                max_iterations: int, optional
                    Max iterations to try converge the function
                tolerance float, optional
                    How close gradient should be to converge

            Returns
            ----------
                float list
                    Point where function has converged
        '''

        converge, point = self.__converge_method(gap, max_iterations, tolerance, True)

        if not converge:
            raise ValueError("Method didn't converge.")

        return point
