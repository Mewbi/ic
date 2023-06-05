import numpy as np

class Function:

    DEFAULT_GAP = 0.1
    DEFAULT_GAP_NUMERICAL = 1e-5
    DEFAULT_MAX_ITERATIONS = 100000
    DEFAULT_TOLERANCE = 1e-5

    def __call__(self, point):
        return self._function(point)

    def __lambda_derivative(self, i):
        return lambda x: self.__derivative_numerical(x, i)
    
    def __init__(self, function, point):
        self.function = function
        self.point = point

        self.derivatives_numerical = []
        for i, _ in enumerate(point):
            self.derivatives_numerical.append(self.__lambda_derivative(i))


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
    def derivatives(self):
        return self._derivatives

    @derivatives.setter
    def derivatives(self, derivatives):
        self._derivatives = []
        if len(derivatives) != len(self._point):
            raise ValueError("Derivatives has not the same dimesion of points")
        try:
            iter(derivatives)
        except:
            raise ValueError("Derivatives is not iterable")
        for derivative in derivatives:
            if not callable(derivative):
                raise ValueError("Derivative {} is not a python callable function".format(derivative))
            self._derivatives.append(derivative)
    
    def __derivative_numerical(self, point, index):
        '''
        Calculate the partial derivative numerically based in derivative definition
        '''
        old = np.copy(point)
        new = np.copy(point)
        new[index] += self.DEFAULT_GAP_NUMERICAL
        f = self._function
        d = ( f(new) - f(old) ) / self.DEFAULT_GAP_NUMERICAL
        return d