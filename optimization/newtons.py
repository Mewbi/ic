import numpy as np
from autograd import grad 
from autograd import hessian

class Function:
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

    DEFAULT_MAX_ITERATIONS = 100000
    DEFAULT_TOLERANCE = 1e-5

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


    # Reference: https://jermwatt.github.io/machine_learning_refined/notes/4_Second_order_methods/4_4_Newtons.html
    def __newtons_method(self, g,max_its,w,**kwargs):
        # compute gradient module using autograd
        gradient = grad(g)
        hess = hessian(g)
    
        # set numericxal stability parameter / regularization parameter
        epsilon = 10**(-7)
        if 'epsilon' in kwargs:
            beta = kwargs['epsilon']

        # run the newtons method loop
        # TODO: Provavelmente vou apagar essas variáveis
        weight_history = [w]           # container for weight history
        cost_history = [g(w)]          # container for corresponding cost function history
        for k in range(max_its):
            # evaluate the gradient and hessian
            grad_eval = gradient(w)
            hess_eval = hess(w)
    
            # reshape hessian to square matrix for numpy linalg functionality
            hess_eval.shape = (int((np.size(hess_eval))**(0.5)),int((np.size(hess_eval))**(0.5)))
        
            # solve second order system system for weight update
            A = hess_eval + epsilon*np.eye(w.size)
            b = grad_eval
            w = np.linalg.solve(A,np.dot(A,w) - b)
        
            # record weight and cost
            weight_history.append(w)
            cost_history.append(g(w))
        return weight_history,cost_history
