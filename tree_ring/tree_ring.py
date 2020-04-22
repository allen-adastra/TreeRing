import sympy as sp
import networkx as nx
import numpy as np
from tree_ring.objects import StateVariable, BasisVariable



def expand(variable_power_mapping, state_variables, disturbance_variables, dependence_graph, moment_basis, reduced_muf = True):
    """
    Args:
        variable_power_mapping (Dictionary): dictionary mapping base variables to powers
        state_variables (set of StateVariable): state variables of the system
        disturbance_variables (set of DisturbanceVariable): disturbance variables of the system
        dependence_graph (networkx graph): nodes are instances of StateVariable, edges represent dependence between two variables
        moment_basis (moment_basis): list of instances of BasisVariable
    """
    # Derive the new relation.
    reduced_vpm = {var : power for var, power in variable_power_mapping.items() if power != 0}
    variable_expansions = [var.update_relation**power for var, power in reduced_vpm.items()]

    # Create a new basis variable and add it to the set of basis variables.
    new_basis_variable = BasisVariable(reduced_vpm, np.prod(variable_expansions))
    moment_basis.add(new_basis_variable)

    # Express as a polynomial.
    state_variables_sympy = [var.sympy_rep for var in state_variables]
    poly = sp.poly(new_basis_variable.update_relation, state_variables_sympy)

    for multi_index in poly.monoms():
        # Iterate over multi-indicies for the monomials.
        # Construct a dictionary mapping base variables to their power in this monomial.
        new_vpm = {state_variables[i] : multi_index[i] for i in range(len(state_variables))}
        if reduced_muf == False:
            # We are finding a completion w.r.t. the un-reduced moment update form.
            update_relation_exists = any([d_var.equivalent_variable_power_mapping(new_vpm) for d_var in moment_basis])
            if update_relation_exists == False:
                expand(new_vpm, state_variables, disturbance_variables, dependence_graph, moment_basis, reduced_muf=reduced_muf)
        else:
            # We are finding a completion w.r.t. the reduced moment update form. 
            # Need to first find a factorization.
            # Find the subgraph induced by multi_index.
            variables_in_mono = {state_variables[i] for i, degree in enumerate(multi_index) if degree != 0}
            dependence_subgraph = dependence_graph.subgraph(variables_in_mono)

            # Find the connected components of the subgraph.
            connected_components = list(nx.connected_components(dependence_subgraph))

            # For each connected component there is a moment. Check if the moment is already in the moment basis.
            for comp in connected_components:
                component_var_power_map = {var : new_vpm[var] for var in comp}
                update_relation_exists = any([d_var.equivalent_variable_power_mapping(component_var_power_map) for d_var in moment_basis])
                if update_relation_exists == False:
                    expand(component_var_power_map, state_variables, disturbance_variables, dependence_graph, moment_basis, reduced_muf=reduced_muf)
    return moment_basis

def generate_code_reps(all_basis_variables, dependence_graph, state_variables, disturbance_variables):
    """

    Args:
        all_basis_variables ([type]): [description]
        dependence_graph ([type]): [description]
        state_variables ([type]): [description]
        disturbance_variables ([type]): [description]
    """

    # Generate a code_rep for every single basis variable.
    for basis_var in all_basis_variables:

        # Express the update in relation as a polynomial in the state and disturbance variables.
        state_variables_sympy = [var.sympy_rep for var in state_variables]
        disturbance_variables_sympy = [var.sympy_rep for var in disturbance_variables]
        update_relation = sp.poly(basis_var.update_relation, state_variables_sympy + disturbance_variables_sympy)
        
        # Indicies of the multi-indices corresponding to state and disturbance variables.
        state_var_idx = range(len(state_variables))
        dist_var_idx = range(len(state_variables), len(state_variables) + len(disturbance_variables))

        # List that will contain the terms in code rep.
        terms_in_code_rep = []

        # Iterate through the terms of the update relation.
        for multi_index in update_relation.monoms():

            # This term will be expressed in terms of the product of 
            # basis variables and disturbance variables expressed in code rep
            # in this list.
            term_code_rep = []

            # First include the relevant disturbance variables in code rep in this list.
            dist_vpm = {disturbance_variables[i - len(state_variables)] : multi_index[i] for i in dist_var_idx if multi_index[i] != 0}
            for dist_var, power in dist_vpm.items():
                term_code_rep.append(dist_var.power(power))

            # Find the connected components in the state dependence graph associated with this term.
            term_vpm = {state_variables[i] : multi_index[i] for i in state_var_idx if multi_index[i] != 0}
            term_vars = list(term_vpm.keys())
            dependence_subgraph = dependence_graph.subgraph(term_vars)
            connected_components = list(nx.connected_components(dependence_subgraph))

            # There should be a one-to-one mapping between connected components and
            # basis variables by virtue of how TreeRing operates. For each connected component,
            # we are essentially finding the equivalent basis variable and adding its sympy rep
            # to the list.
            for comp in connected_components:
                component_var_power_map = {var : term_vpm[var] for var in comp}
                equivalent_variables = [bvar for bvar in all_basis_variables if bvar.equivalent_variable_power_mapping(component_var_power_map)]
                assert len(equivalent_variables) == 1
                equivalent_basis_var = equivalent_variables[0]
                term_code_rep.append(equivalent_basis_var.sympy_rep)

            # If the list is not empty, we can generate the code rep now.
            if len(term_code_rep) > 0:
                basis_rep = np.prod(term_code_rep)
                terms_in_code_rep.append(basis_rep)
        basis_var.code_rep = sum(terms_in_code_rep)