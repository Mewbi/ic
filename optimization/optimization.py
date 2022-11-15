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

    derivates: function list
        Partial derivates of the main function


    Methods
    ----------
    converge_analytics():
        Try to convert a function using his analytics derivates
    
    converge_numerical():
        Try to convert a function using his numerical derivates
    '''

    def __init__(self, function, point):
        self.function = function
        self.point = point

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
        try:
            iter(derivates)
        except:
            raise ValueError("Derivates is not iterable")
        for derivate in derivates:
            if not callable(derivate):
                raise ValueError("Derivate {} is not a python callable function".format(derivate))
            self._derivates.append(derivate)

    def __check_converge(self, pts, tolerance):
        for derivate in self.derivates:
            if abs(derivate(pts)) > tolerance:
                return False
        return True

    def __converge_step(self, f, old_point, old_points, new_points, new_point):
        dist = self.__distance(np.array(old_points)-np.array(new_points)) # Essa maneira não está funcionando, não sei pq

        v = old_point - f(old_points) * (( (old_point-new_point) / ( f(old_points) - f(new_points) ) ))
        print("Ponto encontrado: {}".format( v) )
        return v

    def __distance(self, points):
        s = 0
        for p in points:
            s += p**2
        return np.sqrt(s)

    def __converge_method(self, gap, n, tolerance):
        old = self._point
        for _ in range(n):
            new = []
            for i, derivate in enumerate(self.derivates):
                old_point = old[i]
                #old_points = np.array(old[:])
                new_points = old[:]
                new_points[i] += gap
                new.append(self.__converge_step(derivate, old_point, old, new_points, new_points[i]))
                
            old = new[:]
            print(new)
            if self.__check_converge(new, tolerance):
                print("Convergiu em {}".format(new))
                return True, new[:]

        return False, []

    def converge_analytical(self, gap = 0.01, max_iterations = 100000, tolerance = 1e-5):
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



        converge, point = self.__converge_method(gap, max_iterations, tolerance)

        if not converge:
            raise ValueError("Method didn't converge.")

        return point
