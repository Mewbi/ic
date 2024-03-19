import csv
import math
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Dict
from matplotlib.animation import FuncAnimation

from numpy import linalg
from mpl_toolkits.mplot3d import Axes3D

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

    specie: string
        The specie of optimization

    variable: string
        Variable that was changed from optimal setting
    
    variation: number
        Variation of variable from optimal setting

    convergence_steps: float lists
        Each step os convergence process

    Methods
    ----------
    plot():
        Plot molecular configuration in final point
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
                 convergence_steps = [],
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
        self.convergence_steps = convergence_steps

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
        if type(variation) is not int and type(variation) is not float and not np.isreal(variation):
            raise ValueError("Variation value must be a number")
        self._variation = variation

    @property
    def convergence_steps(self):
        return self._convergence_steps

    @convergence_steps.setter
    def convergence_steps(self, steps):
        try:
            iter(steps)
        except:
            raise ValueError("Steps must be iterable")
        
        for step in steps:
            try:
                iter(step)
            except:
                raise ValueError("Each step must be iterable")

        self._convergence_steps = steps

    def __calculate_coords(self, p):
        if len(p) != 6:
            raise ValueError("Probably the result is not from the F+H20 reaction, there will be no plot")

        R1 = p[0]
        R2 = p[1]
        R3 = p[2]
        A1 = p[3]
        A2 = p[4]
        D1 = p[5]

        # Convert to radians
        A1_rad = np.radians(A1)
        A2_rad = np.radians(A2)
        D1_rad = np.radians(D1)

        # Elements
        # Colors: https://en.wikipedia.org/wiki/CPK_coloring
        H1 = Element(color="whitesmoke",  size=1, label="H1")
        H2 = Element(color="whitesmoke",  size=1, label="H2")
        O  = Element(color="red",        size=1, label="O")
        F  = Element(color="lightgreen", size=1, label="F")

        # Calculate points coordenates
        O.X = 0
        O.Y = 0
        O.Z = 0

        H1.X = 0
        H1.Y = R1
        H1.Z = 0
    
        H2.X = np.sin(A1_rad) * R2
        H2.Y = np.cos(A1_rad) * R2
        H2.Z = 0

        F.Y = (R1 - R3 * np.cos(A2_rad)) / R1
        # R3 ** 2 = F.X ** 2 + (F.Y - H1.Y) ** 2 + F.Z ** 2
        # np.cos(D1_rad) ** 2 = ( -F.X ** 2 ) / (F.X ** 2 + F.Z ** 2)
        m = np.array([
            [- np.sin(D1_rad) ** 2, np.cos(D1_rad) ** 2],
            [1, 1]
        ])

        m2 = np.array(
            [0, R3 ** 2 - (F.Y -R1) ** 2]
        )

        x = linalg.inv(m) @ m2

        F.X = np.sqrt(x[0])
        F.Z = np.sqrt(x[1])

        return H1, H2, O, F

    def __setup_plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_facecolor('darkgrey')
        ax.axis('off')

        return fig, ax

    def plot(self):

        p = self.final_point

        H1, H2, O, F = self.__calculate_coords(p)

        x = np.array([H1.X, H2.X, O.X, F.X])
        y = np.array([H1.Y, H2.Y, O.Y, F.Y])
        z = np.array([H1.Z, H2.Z, O.Z, F.Z])
        
        # Configure 3D plot
        fig, ax = self.__setup_plot()

        # Plot points
        ax.scatter(H1.X, H1.Y, H1.Z, c=H1.color, marker='o', s=H1.size)
        ax.scatter(H2.X, H2.Y, H2.Z, c=H2.color, marker='o', s=H2.size)
        ax.scatter(O.X,  O.Y,  O.Z,  c=O.color,  marker='o', s=O.size)
        ax.scatter(F.X,  F.Y,  F.Z,  c=F.color,  marker='o', s=F.size)

        ax.text(H1.X, H1.Y, H1.Z, H1.label, fontweight='bold')
        ax.text(H2.X, H2.Y, H2.Z, H2.label, fontweight='bold')
        ax.text(O.X,  O.Y,  O.Z,  O.label,  fontweight='bold')
        ax.text(F.X,  F.Y,  F.Z,  F.label,  fontweight='bold')

        # Plot lines between points
        ax.plot([H1.X, O.X], [H1.Y, O.Y], [H1.Z, O.Z], c='b')
        ax.plot([H2.X, O.X], [H2.Y, O.Y], [H2.Z, O.Z], c='b')
        ax.plot([H1.X, F.X], [H1.Y, F.Y], [H1.Z, F.Z], c='b')
        
        # Adjust scale
        ax.auto_scale_xyz([np.min(x), np.max(x)], [np.min(y), np.max(y)], [np.min(z), np.max(z)])

        # Plot
        plt.show()

    def save_steps(self, filename="convergence.gif"):
        if len(self.convergence_steps) == 0:
            raise ValueError("Empty convergence steps")

        fig, ax = self.__setup_plot()
        p = self.init_point
        H1, H2, O, F = self.__calculate_coords(p)

        # Dados dos pontos
        x = np.array([H1.X, H2.X, O.X, F.X])
        y = np.array([H1.Y, H2.Y, O.Y, F.Y])
        z = np.array([H1.Z, H2.Z, O.Z, F.Z])

        def init():
            ax.scatter(H1.X, H1.Y, H1.Z, c=H1.color, marker='o', s=H1.size)
            ax.scatter(H2.X, H2.Y, H2.Z, c=H2.color, marker='o', s=H2.size)
            ax.scatter(O.X,  O.Y,  O.Z,  c=O.color,  marker='o', s=O.size)
            ax.scatter(F.X,  F.Y,  F.Z,  c=F.color,  marker='o', s=F.size)

            ax.text(H1.X, H1.Y, H1.Z, H1.label, fontweight='bold')
            ax.text(H2.X, H2.Y, H2.Z, H2.label, fontweight='bold')
            ax.text(O.X,  O.Y,  O.Z,  O.label,  fontweight='bold')
            ax.text(F.X,  F.Y,  F.Z,  F.label,  fontweight='bold')

            ax.plot([H1.X, O.X], [H1.Y, O.Y], [H1.Z, O.Z], c='b')
            ax.plot([H2.X, O.X], [H2.Y, O.Y], [H2.Z, O.Z], c='b')
            ax.plot([H1.X, F.X], [H1.Y, F.Y], [H1.Z, F.Z], c='b')

            ax.auto_scale_xyz([np.min(x), np.max(x)], [np.min(y), np.max(y)], [np.min(z), np.max(z)])

        def update(step):
            nH1, nH2, nO, nF = self.__calculate_coords(step)
            H1.X = nH1.X
            H1.Y = nH1.Y
            H1.Z = nH1.Z

            H2.X = nH2.X
            H2.Y = nH2.Y
            H2.Z = nH2.Z

            O.X = nO.X
            O.Y = nO.Y
            O.Z = nO.Z

            F.X = nF.X
            F.Y = nF.Y
            F.Z = nF.Z
            ax.clear()
            init()

        animation = FuncAnimation(fig, update, frames=self.convergence_steps, init_func=init, blit=False)
        animation.save(filename, writer='pillow', fps=2)


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

    add_single_result():
        Add a single result to results list

    add_multiple_results():
        Add multiple results to results list

    normalize_final_points():
        Remove different significant figures
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
    def results(self) -> List[Result]:
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

    def add_single_result(self, result: Result):
        if type(result) is not Result:
            raise ValueError("Result must a be a Result object")

        if self.results == None:
            self.results = [result]
            return

        self.results.append(result)

    def add_multiple_results(self, results: List[Result]):
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

            for j, value in enumerate(point):
                if value > max_values[j]:
                    max_values[j] = value

                if value < min_values[j]:
                    min_values[j] = value

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

            for j, value in enumerate(point):
                value = self.__truncate(value, trunc_decimals[j])
                point[j] = value

            result.final_point = point
            results[i] = result

        self.results = results

    def get_results_metrics(self) -> Dict:
        '''
        Return result metrics

        Convergence is related to each stactionaty point

        Iterations is related only to convergence success
        Key means number of iterations necessary to converge
        Value means number of cases that it occoured
        '''
        metrics = {}

        # Convergence tax
        conv = {}
        for result in self.results:
            specie = result.specie
            if specie not in conv:
                conv[specie] = {
                        "conv_success": 0,
                        "conv_fail": 0,
                        "total": 0,
                        "conv_percent": 0,
                    }

            conv[specie]["total"] += 1
            if result.converge:
                conv[specie]["conv_success"] += 1
            else:
                conv[specie]["conv_fail"] += 1

        for specie, values in conv.items():
            values["conv_percent"] = (values["conv_success"] / values["total"]) * 100

        metrics["convergence"] = conv

        # Get total convergence
        total = 0
        for result in self.results:
            if result.converge:
                total += 1

        # Iterations cases
        iterations = {}
        for result in self.results:
            if not result.converge:
                continue

            it = result.iterations
            if it not in iterations:
                iterations[it] = {
                        "qtd": 0,
                        "percent": 0,
                    }

            iterations[it]["qtd"] += 1

        for it, data in iterations.items():
            data["percent"] = (data["qtd"] / total) * 100
            iterations[it] = data

        metrics["iterations"] = iterations
        return metrics


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

class Element:
    
    def __init__(self, color="red", size=1, label=""):
        self.X = 0
        self.Y = 0
        self.Z = 0
        self.color = color
        self.size = 20 * 4 * size
        self.label = label

    @property
    def X(self):
        return self._X

    @X.setter
    def X(self, value):
        self.__valid_point(value, "X")
        self._X = value

    @property
    def Y(self):
        return self._Y

    @Y.setter
    def Y(self, value):
        self.__valid_point(value, "Y")
        self._Y = value

    @property
    def Z(self):
        return self._Z

    @Z.setter
    def Z(self, value):
        self.__valid_point(value, "Z")
        self._Z = value

    def __valid_point(self, value, coord):
        if type(value) is not int and type(value) is not float and not np.isreal(value):
            raise ValueError("Value of {} must be a number. Receive: {}".format(coord, type(value)))

    @property
    def color(self):
        return self._color

    @color.setter
    def color(self, c):
        if type(c) is not str:
            raise ValueError("Color must be a string")

        c = c.strip()

        if len(c) == 0:
            raise ValueError("Color cannot be a empty string")

        self._color = c

    @property
    def label(self):
        return self._label

    @label.setter
    def label(self, l):
        if type(l) is not str:
            raise ValueError("Label must be a string")

        self._label = l
