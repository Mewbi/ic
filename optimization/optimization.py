class Function:
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
            #print("Ponto: {} - Derivada {} - Função {}\nPonto: {} - Derivada {} - Função {}".format(old_points, f(old_points), self(old_points), new_points, f(new_points), self(new_points)))
            return old_point
        return v

    def __converge_method(self, gap, n, tolerance, numerical):
        derivatives = self.derivatives_numerical if numerical else self.derivatives
        old = self._point
        for iteration in range(n):
            new = []
            for i, derivative in enumerate(derivatives):
                old_point = old[i]
                new_points = old[:]
                new_points[i] += gap
                new.append(self.__converge_step(derivative, old_point, old, new_points, new_points[i]))
            
            #print(new)
            old = new[:]
            if self.__check_converge(new, derivatives, tolerance):
                return True, new[:], iteration + 1

        return False, [], 0

    def converge_analytically(self, gap = DEFAULT_GAP, max_iterations = DEFAULT_MAX_ITERATIONS, tolerance = DEFAULT_TOLERANCE):
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

            Returns
            ----------
                float list
                    Point where function has converged
        '''

        if not self.derivatives:
            raise ValueError("Partial derivatives not defined.")

        converge, point, iterations = self.__converge_method(gap, max_iterations, tolerance, False)

        if not converge:
            raise ValueError("Method didn't converge.")

        return point, iterations


    def __derivative_numerical(self, point, index):
        '''
        Calculate the partial derivative numerically based in derivative definition
        '''
        old = point[:]
        new = point[:]
        new[index] += self.DEFAULT_GAP_NUMERICAL
        f = self._function
        d = ( f(new) - f(old) ) / self.DEFAULT_GAP_NUMERICAL
        return d

    def converge_numerically(self, gap = DEFAULT_GAP, max_iterations = DEFAULT_MAX_ITERATIONS, tolerance = DEFAULT_TOLERANCE):
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

            Returns
            ----------
                float list
                    Point where function has converged
        '''

        converge, point, iterations = self.__converge_method(gap, max_iterations, tolerance, True)

        if not converge:
            raise ValueError("Method didn't converge.")

        return point, iterations
