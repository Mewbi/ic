import csv
import math
import numpy as np

class Function:

    DEFAULT_GAP = 0.1
    DEFAULT_GAP_NUMERICAL = 1e-5
    DEFAULT_MAX_ITERATIONS = 100000
    DEFAULT_TOLERANCE = 1e-5
    FAIL_CONVERGE_VALUE = 115.302739153 # When PES value 'explode'

    def __call__(self, point):
        return self._function(point)

    def __lambda_derivative(self, i):
        return lambda x: self.__derivative_numerical(x, i)
    
    def __init__(self, function, point, optimize_vars = []):
        self.function = function
        self.point = point
        self.dimension = len(point)
        self.optimize_vars = optimize_vars

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
    def dimension(self):
        return self._dimension
    
    @dimension.setter
    def dimension(self, dimension):
        self._dimension = dimension

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
    
    @property
    def optimize_vars(self):
        return self._optimize_vars

    @optimize_vars.setter
    def optimize_vars(self, optimize):
        if type(optimize) is not list:
            raise ValueError("Optimize vars must be a list of boolean")

        if len(optimize) == 0:
            optimize = [True for _ in range(self.dimension)]

        if len(optimize) != self.dimension:
            raise ValueError("Optimize vars must have the same length as the function dimension")

        self._optimize_vars = optimize

    def __derivative_numerical(self, point, index):
        '''
        Calculate the partial derivative numerically based in derivative definition
        '''
        old = np.copy(point)
        new = np.copy(point)
        old[index] -= self.DEFAULT_GAP_NUMERICAL
        new[index] += self.DEFAULT_GAP_NUMERICAL
        f = self._function
        d = ( f(new) - f(old) ) / ( 2 * self.DEFAULT_GAP_NUMERICAL )
        return d
    
    def __derivative_numerical_bkp(self, point, index):
        '''
        Calculate the partial derivative numerically based in derivative definition
        '''
        old = np.copy(point)
        new = np.copy(point)
        new[index] += self.DEFAULT_GAP_NUMERICAL
        f = self._function
        d = ( f(new) - f(old) ) / self.DEFAULT_GAP_NUMERICAL
        return d

class Result:
    '''
    Result of convergence process

    ...

    Attributes
    ----------
    converge: bool
        True if function has converged and False if not

    init_point: float list
        Initial point of optimization
    
    final_point: float list
        Point where function has converged or stopped

    iterations: int
        Number of iterations

    gradient: float
        Gradient value of function

    init_value: float
        Init value of function

    final_value: float 
        Final value of function
    '''


    def __init__(self, 
                 converge = False,
                 init_point = [],
                 final_point = [],
                 init_value = 0.0,
                 final_value = 0.0,
                 iterations = 0,
                 gradient = "0",
                 specie = "",
                 variable = "",
                 variation = 0,
                 ):
        self.converge = converge
        self.init_point = init_point
        self.final_point = final_point
        self.init_value = init_value
        self.final_value = final_value
        self.iterations = iterations
        self.gradient = gradient

        self.specie = specie
        self.variable = variable
        self.variation = variation

    @property
    def converge(self):
        return self._converge

    @converge.setter
    def converge(self, converge):
        if type(converge) is not bool:
            raise ValueError("Converge value is not a boolean")
        self._converge = converge

    @property
    def init_point(self):
        return self._init_point

    @init_point.setter
    def init_point(self, init_point):
        try:
            iter(init_point)
        except:
            raise ValueError("Init point is not iterable")
        self._init_point = init_point

    @property
    def final_point(self):
        return self._final_point

    @final_point.setter
    def final_point(self, final_point):
        try:
            iter(final_point)
        except:
            raise ValueError("Final point is not iterable")
        
        if len(final_point) != len(self.init_point):
            raise ValueError("Final point must have same dimension as init point")

        self._final_point = final_point

    @property
    def init_value(self):
        return self._init_value

    @init_value.setter
    def init_value(self, init_value):
        if type(init_value) is not float:
            raise ValueError("Init value is not a float number")
        self._init_value = init_value

    @property
    def final_value(self):
        return self._final_value

    @final_value.setter
    def final_value(self, final_value):
        if type(final_value) is not float:
            raise ValueError("Final value is not a float number")
        self._final_value = final_value

    @property
    def iterations(self):
        return self._iterations

    @iterations.setter
    def iterations(self, iterations):
        if type(iterations) is not int:
            raise ValueError("Iterations value is not a integer number")
        self._iterations = iterations

    @property
    def gradient(self):
        return self._gradient

    @gradient.setter
    def gradient(self, gradient):
        if type(gradient) is not str:
            raise ValueError("Gradient value must be a string value")
        self._gradient = gradient

    @property
    def specie(self):
        return self._specie

    @specie.setter
    def specie(self, specie):
        if type(specie) is not str:
            raise ValueError("Specie value must be a string value")
        self._specie = specie

    @property
    def variable(self):
        return self._variable

    @variable.setter
    def variable(self, variable):
        if type(variable) is not str:
            raise ValueError("Variable value must be a string value")
        self._variable = variable

    @property
    def variation(self):
        return self._variation

    @variation.setter
    def variation(self, variation):
        if type(variation) is not int and type(variation) is not float:
            raise ValueError("Variation value must be a number")
        self._variation = variation

class Results():
    '''
    Results of many optimizations process

    ...

    Attributes
    ----------
    results: Result list
        List of Result objects

    Methods
    ----------
    csv():
        Generate a CSV file with results values

    json():
        Generate a JSON file with results values
    '''

    FIELDS_ORDER = [
        "specie",
        "variable",
        "variation",
        "converge",
        "iterations",
        "gradient",
        "init_value",
        "final_value",
        "init_point",
        "final_point"
    ]

    VARIATION_TOLERANCE = 5e-1

    def __init__(self):
        self.results = None

    @property
    def results(self):
        return self._results

    @results.setter
    def results(self, results):
        if results == None:
            self._results = results
            return

        try:
            iter(results)
        except:
            raise ValueError("Results must be iterable")

        for result in results:
            if type(result) is not Result:
                raise ValueError("Values of results list must a a Result object")
        self._results = results

    def add_single_result(self, result):
        if type(result) is not Result:
            raise ValueError("Result must a be a Result object")

        if self.results == None:
            self.results = [result]
            return

        self.results.append(result)

    def add_multiple_results(self, results):
        try:
            iter(results)
        except:
            raise ValueError("Results must be iterable")

        for result in results:
            self.add_single_result(result)

    def normalize_final_points(self, refence_value=None):
        max_values     = np.copy([])
        min_values     = np.copy([])
        trunc_decimals = np.copy([])

        results = self.results
        for result in results:
            if result.converge == False:
                continue
            
            if refence_value is not None and abs(refence_value - result.final_value) > self.VARIATION_TOLERANCE:
                continue

            point = result.final_point

            if max_values.size == 0:
                max_values = np.copy(point)

            if min_values.size == 0:
                min_values = np.copy(point)

            for i, value in enumerate(point):
                if value > max_values[i]:
                    max_values[i] = value

                if value < min_values[i]:
                    min_values[i] = value

        for i, max_v in enumerate(max_values):
            min_v = min_values[i]
            diff = max_v - min_v

            if diff == 0:
                trunc = len(str(max_v).split('.')[1])
                trunc_decimals = np.append(trunc_decimals, trunc)
                continue

            trunc = 0
            while True:
                if math.trunc(diff * (10 ** trunc)) != 0:
                    break
                trunc += 1

            trunc_decimals = np.append(trunc_decimals, trunc - 1)


        for i, result in enumerate(results):
            if result.converge == False:
                continue

            if refence_value is not None and abs(refence_value - result.final_value) > self.VARIATION_TOLERANCE:
                continue

            point = result._final_point

            for i, value in enumerate(point):
                value = self.__truncate(value, trunc_decimals[i])
                point[i] = value

            result.final_point = point
            results[i] = result

        self.results = results


    def __truncate(self, number, digits):
        nb_decimals = len(str(number).split('.')[1])
        if nb_decimals <= digits:
            return number
        stepper = 10.0 ** digits
        return math.trunc(stepper * number) / stepper

    def csv(self, filename = 'result.csv'):
        with open(filename, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(self.FIELDS_ORDER)

            for result in self.results:
                row = []
                for field in self.FIELDS_ORDER:
                    value = getattr(result, field)
                    row.append(value)
                writer.writerow(row)
