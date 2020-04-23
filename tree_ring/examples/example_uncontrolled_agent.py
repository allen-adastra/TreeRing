from tree_ring.tree_ring import expand, generate_code_reps
import tree_ring.objects as tro
import sympy as sp
import numpy as np
import networkx as nx

class UncontrolledAgent(object):
    def __init__(self):
        # Declare SymPy variables.
        x = sp.Symbol("x")
        y = sp.Symbol("y")
        v = sp.Symbol("v")
        wv = sp.Symbol("wv")
        theta = sp.Symbol("theta")
        cos_theta = sp.Symbol("c")
        sin_theta = sp.Symbol("s")
        cos_wtheta = sp.Symbol("cw")
        sin_wtheta = sp.Symbol("sw")

        # Initialize state variables with discrete time dynamics.
        self._x = tro.StateVariable(x, x + v * cos_theta)
        self._yt = tro.StateVariable(y, y + v * sin_theta)
        self._vt = tro.StateVariable(v, v + wv)
        self._sin_thetat = tro.StateVariable(sin_theta, sin_theta * cos_wtheta + cos_theta * sin_wtheta)
        self._cos_thetat = tro.StateVariable(cos_theta, cos_theta * cos_wtheta - sin_theta * sin_wtheta)
        self._state_variables = [self._x, self._yt, self._vt, self._sin_thetat, self._cos_thetat]

        # Initialize disturbance variables.
        self._wvt = tro.DisturbanceVariable(wv)
        self._sin_wthetat = tro.DisturbanceVariable(sin_wtheta)
        self._cos_wthetat = tro.DisturbanceVariable(cos_wtheta)
        self._disturbance_variables = [self._sin_wthetat, self._cos_wthetat, self._wvt]

        # Specify dependence graph.
        self._dependence_graph = nx.Graph()
        self._dependence_graph.add_nodes_from(self._state_variables)
        self._dependence_graph.add_edges_from([(self._x, self._yt), (self._x, self._vt), (self._yt, self._vt),
                                                        (self._x, self._sin_thetat), (self._x, self._cos_thetat), (self._yt, self._sin_thetat),
                                                        (self._yt, self._cos_thetat), (self._cos_thetat, self._sin_thetat)])
        
    def test_reduced(self):
        # List of moments we want to propagate represented as dictionaries, in this case, the second moment of position.
        # Each dictionary maps the the variables in the moment to each power. Below, we have that:
        #   E[x^2] is represented by {self._x : 2}
        #   E[y^2] is represented by {self._yt : 2}
        #   E[xy] is represented by {self._x : 1, self._yt : 1}
        moments_to_propagate = [{self._x : 2}, {self._yt : 2}, {self._x : 1, self._yt : 1}]

        # Initialize the moment basis with the state variables.
        moment_basis = {var.to_basis_variable() for var in self._state_variables}
        for variable_power_map in moments_to_propagate:
            expand(variable_power_map, self._state_variables, self._disturbance_variables, self._dependence_graph, moment_basis, reduced_muf=True)
        #print_moment_basis(moment_basis)
        generate_code_reps(moment_basis, self._dependence_graph, self._state_variables, self._disturbance_variables)
        generate_underactuated_code(moment_basis)

    def test_unreduced(self):
        # List of moments we want to propagate represented as dictionaries, in this case, the second moment of position.
        # Each dictionary maps the the variables in the moment to each power. Below, we have that:
        #   E[x^2] is represented by {self._x : 2}
        #   E[y^2] is represented by {self._yt : 2}
        #   E[xy] is represented by {self._x : 1, self._yt : 1}
        moments_to_propagate = [{self._x : 2}, {self._yt : 2}, {self._x : 1, self._yt : 1}]

        # Initialize the moment basis using state variables.
        moment_basis = {var.to_basis_variable() for var in self._state_variables}
        for variable_power_map in moments_to_propagate:
            expand(variable_power_map, self._state_variables, self._disturbance_variables, self._dependence_graph, moment_basis, reduced_muf=False)
        print_moment_basis(moment_basis)
        # TODO: generate_code_rep does not work for the unreduced case.
        generate_code_reps(moment_basis, self._dependence_graph, self._state_variables, self._disturbance_variables)

def generate_underactuated_code(moment_basis):
    for i, basis_var in enumerate(moment_basis):
        print(str(basis_var.sympy_rep) + " = moment_state[t, " + str(i) + "]")
    print("\n")
    for i, basis_var in enumerate(moment_basis):
        print("prog.AddConstraint(moment_state[t + 1, " + str(i) + "] == " + str(basis_var.code_rep) + ")")
    print("\n")
    for i, basis_var in enumerate(moment_basis):
        print("state[" + str(i) + "] = " + str(basis_var.state_variable_rep()))

def print_moment_basis(moment_basis):
    print("Printing variables of the moment basis and their expansions.")
    for var in moment_basis:
        print(str(var.sympy_rep) + " expansion: ")
        print(sp.expand(var.update_relation))
        print("\n")

test = UncontrolledAgent()
test.test_reduced()