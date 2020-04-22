import sympy as sp
import numpy as np

class StateVariable(object):
    """
    A state variable of the polynomial system x_{t + 1} = f(x_t, w_t).
    """
    def __init__(self, sympy_rep, dynamics):
        """
        Args:
            sympy_rep (Instance of Sp.Symbol): Representation of this variable in SymPy.
            dynamics (SymPy expression): The dynamics of this state variable.
        """
        self._sympy_rep = sympy_rep
        self._update_relation = dynamics

    def __str__(self):
        return str(self._sympy_rep)

    def __repr__(self):
        return repr(self._sympy_rep)

    @property
    def update_relation(self):
        return self._update_relation
    
    @property
    def sympy_rep(self):
        return self._sympy_rep

class DisturbanceVariable(object):
    """
    A disturbance variable of the polynomial system x_{t + 1} = f(x_t, w_t).
    """
    def __init__(self, sympy_rep):
        """
        Args:
            sympy_rep (Instance of Sp.Symbol): Representation of this variable in SymPy.
        """
        self._sympy_rep = sympy_rep
    
    def __str__(self):
        return str(self._sympy_rep)

    def __repr__(self):
        return repr(self._sympy_rep)

    @property
    def sympy_rep(self):
        return self._sympy_rep

    def power(self, power):
        """
        Return a sympy variable representing this variable raised to a power in the notation:
            variable_name + power
        So for example:
            wv2
        Args:
            pow (int): [description]

        Returns:
            sp.Symbol: [description]
        """
        return sp.Symbol(self._sympy_rep.name + str(power))


class BasisVariable(object):
    def __init__(self, variable_power_mapping, update_relation):
        """
        Args:
            variable_power_mapping (Dictionary Mapping BaseVariable -> Natural Number):
                Every derived variable is of the form x^alpha where x is a vector and alpha is a multi-index.
                variable_power_mapping essentially maps x_i to alpha_i for all i s.t. alpha_i > 0.
            update_relation (SymPy expression): This variable in at time t + 1 as a polynomial in base variables at time t.
        """
        assert min(variable_power_mapping.values()) > 0
        self._variable_power_mapping = {var : power for var, power in variable_power_mapping.items() if power > 0}
        self._update_relation = update_relation
        self._sympy_rep = self.generate_sympy_rep()
        self._update_relation_code_rep = None

    def __hash__(self):
        return hash(self.sympy_rep)
    
    def __eq__(self, other_obj):
        return hash(self.sympy_rep) == hash(other_obj.sympy_rep)

    @property
    def sympy_rep(self):
        return self._sympy_rep

    @property
    def variable_power_mapping(self):
        return self._variable_power_mapping

    @property
    def update_relation(self):
        return self._update_relation

    @property
    def update_relation_code_rep(self):
        return self._code_rep
    
    def generate_sympy_rep(self):
        """
        Generate the sympy rep of the basis variable using the variables in self._variable_power_mapping.
        """
        string = ''
        for var, power in self._variable_power_mapping.items():
            string += var.sympy_rep.name + str(power) + "_"
        # Remove the underscore at the end.
        string = string.strip("_")
        return sp.Symbol(string)

    def equivalent_variable_power_mapping(self, variable_power_map):
        """
        Given another variable power mapping, check if it is equal to that of this variable.
        Args:
            variable_power_mapping (Dictionary Mapping BaseVariable -> Natural Number): variable power mapping to check.
        """
        reduced_vpm = {var : power for var, power in variable_power_map.items() if power > 0}
        if reduced_vpm == self._variable_power_mapping:
            return True