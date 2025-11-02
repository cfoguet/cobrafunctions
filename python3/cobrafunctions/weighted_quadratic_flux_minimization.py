import time
import cobra
import pandas as pd
from cobra.util import solver as sutil
from optlang.symbolics import Zero as optlang_Zero, add #, Pow

#TODO Update qMTA to use this function

def add_quadratic_difference_minimization(model, target_fluxes={},target_fluxes_weight={},expanded_reaction_mapping_dict={},copy_model=True,verbose=False,):
    r"""
    ###Adapated from add_moma in cobrapy
    ###WARNING: Function will not work if there is only one target flux and its weight is 1
    Inputs are:
    model : cobra.Model
        The model to add Flux Difference minimization constraints and objective to. 
    target_fluxes : dict
        A dictionary of target fluxes where keys are fluxes ids and values are the target flux values. If flux to reaction mapping is not provided, flux ids are assumed to be the same as reaction ids 
    forward_reverse_net_flux_dict : dict
        A dictionary of weights for each target flux where keys are reaction ids and values are the weights. If flux to reaction mapping is not provided, flux ids are assumed to be the same as reaction ids 
    expanded_reaction_mapping_dict : dict
        A dictionary mapping flux ids to forward and reverse reactions in the model. If not provided, flux ids are assumed to be the same as reaction ids and all reactions are standard reactions (no forward/reverse split). The format is:
        {
           "flux_id1":{"forward_reactions":["reaction_id1","reaction_id2"],"reverse_reactions":["reaction_id3"]},
           "flux_id2":{"forward_reactions":["reaction_id4"],"reverse_reactions":[]},
           ...
        }       
    verbose : bool
        Whether to print verbose output or not. Default is False.             
    """    
    #CF added
    if isinstance(target_fluxes, pd.Series):
       target_fluxes = target_fluxes.to_dict()
    if(copy_model):
       model=model.copy()
    #  if forward_reverse_net_flux_dict is None assumed flux id is the same as reaction id and all reactions are standard reactions
    if expanded_reaction_mapping_dict is None or len(expanded_reaction_mapping_dict)==0: 
       expanded_reaction_mapping_dict={x:{"forward_reactions":[x],"reverse_reactions":[]} for x in target_fluxes.keys()}
    
    #Scale weights 
    scaling_factor=max([abs(target_fluxes_weight[x])for x in target_fluxes_weight])
    if(verbose):
       print("Scaling factor for weights: "+str(scaling_factor))
    #End CF added
    if "old_objective" in model.solver.variables:
        raise ValueError("The model is already adjusted for Flux Difference")
    # Fall back to default QP solver if current one has no QP capability
    if  sutil.interface_to_str(model.problem) not in sutil.qp_solvers:
        model.solver = sutil.choose_solver(model, qp=True)
    #if solution is None:
    #    solution = pfba(model)
    prob = model.problem
    v = prob.Variable("old_objective")
    c = prob.Constraint(
        model.solver.objective.expression - v,
        lb=0.0,
        ub=0.0,
        name="old_objective_constraint",
    )
    to_add = [v, c]
    model.objective = prob.Objective(optlang_Zero, direction="min", sloppy=True)
    obj_vars = []
    for flux_id in target_fluxes:
        target_flux_value = target_fluxes[flux_id]
        weight=target_fluxes_weight[flux_id]/scaling_factor
        
        forward_reactions=expanded_reaction_mapping_dict[flux_id]["forward_reactions"]
        reverse_reactions=expanded_reaction_mapping_dict[flux_id]["reverse_reactions"]

        #Build the net flux expression
        #As a default cobra splits every reactions into forward and reverse variables defined by r.flux_expression and r.forward_variable and r.reverse_variable and r.reverse_variable
        #We will make our own net flux expression to handle the case of forward/reverse reactions
        net_flux_expression=optlang_Zero
        ##Forward/Standard reaction
        for forward_rid in forward_reactions:
           r=model.reactions.get_by_id(forward_rid)
           if r.forward_variable.ub>0:
              net_flux_expression+=r.forward_variable
           if r.reverse_variable.ub>0:
              if len(reverse_reactions)>0:
                 raise Exception(forward_rid+":If using Forwards/Reverse reactions forward reaction should not have negative lower bound")
              net_flux_expression-=r.reverse_variable
        ##Reverse reaction if applicable
        for reverse_rid in reverse_reactions:
           r_rev=model.reactions.get_by_id(reverse_rid)
           if r_rev.forward_variable.ub>0:
              net_flux_expression-=r_rev.forward_variable
           if r_rev.reverse_variable.ub>0:
              raise Exception(reverse_rid+":Reverse reactions should not have reverse fluxes allowed")
              #net_flux_expression+=r_rev.reverse_variable
        if net_flux_expression is optlang_Zero:
           if(verbose):
              print("Skipping reaction "+flux_id+" as no flux is allowed in either direction")
           continue #Skip reactions with no flux allowed
        #Define net flux variable and constraint
        net_flux = prob.Variable("net_flux_" + flux_id)
        net_flux_const = prob.Constraint(
                net_flux_expression - net_flux,
                lb=0,
                ub=0,
                name="Cnet_flux" + flux_id,
            )
        to_add.extend([net_flux, net_flux_const])
        #Define quadratic expression
        # weight * (net_flux - target_flux_value)^2= weight * net_flux**2 - 2*weight*target_flux_value*net_flux + weight*target_flux_value**2
        #The last term is constant and can be ignored in the optimization        
        obj_vars.append(weight * net_flux**2 - 2*weight*target_flux_value*net_flux)
        if(verbose):
           print("######################## "+flux_id)
           print("Adding flux difference minimization for reaction: "+flux_id+" with target flux: "+str(target_flux_value)+" and weight: "+str(weight) )
           print("Net flux expression: "+str(net_flux_expression) )
           print("Objective term: "+str(weight * net_flux**2 - 2*weight*target_flux_value*net_flux) )
    model.add_cons_vars(to_add)
    model.objective = prob.Objective(add(obj_vars), direction="min", sloppy=True)
    #print(model.objective.expression)      
    return model 

"""
def add_quadratic_difference_minimization(model, target_fluxes={},target_fluxes_weight={},copy_model=True,):
    #CF added
    if isinstance(target_fluxes, pd.Series):
       target_fluxes = target_fluxes.to_dict()
    if(copy_model):model=model.copy()
    #Scale weights 
    scaling_factor=max([abs(target_fluxes_weight[x])for x in target_fluxes_weight])
    
    #End CF added
    if "old_objective" in model.solver.variables:
        raise ValueError("The model is already adjusted for Flux Difference")
    # Fall back to default QP solver if current one has no QP capability
    if  sutil.interface_to_str(model.problem) not in sutil.qp_solvers:
        model.solver = sutil.choose_solver(model, qp=True)
    #if solution is None:
    #    solution = pfba(model)
    prob = model.problem
    v = prob.Variable("old_objective")
    c = prob.Constraint(
        model.solver.objective.expression - v,
        lb=0.0,
        ub=0.0,
        name="old_objective_constraint",
    )
    to_add = [v, c]
    model.objective = prob.Objective(optlang_Zero, direction="min", sloppy=True)
    obj_vars = []
    for r_id in target_fluxes:
        r=model.reactions.get_by_id(r_id)
        flux = target_fluxes[r_id]
        weight=target_fluxes_weight[r_id]/scaling_factor
        #As a default cobra splits every reactions into forward and reverse reactions defined by r.flux_expression
        
        
        dist = prob.Variable("net_flux" + r.id)
        const = prob.Constraint(
                r.flux_expression - dist,
                lb=0,
                ub=0,
                name="net_flux" + r.id,
            )
        to_add.extend([dist, const])
        obj_vars.append(weight * dist**2 - 2*weight*flux*dist)
    model.add_cons_vars(to_add)
    model.objective = prob.Objective(add(obj_vars), direction="min", sloppy=True)
    #print(model.objective.expression)      
    return model 


"""

"""
Alternative formulation. Seemed to be slower than the previous one
def add_weighted_difference_minimization(model, target_fluxes={},target_fluxes_weight={}, linear=False,copy_model=True,):
    
    ###Adapated from cobrapy
    ###WARNING: Function will not work if there is only one target flux and its weight is 1
    ###Adapated from add_moma
    Add MOMA constraints and objective representing to the `model`.

    This adds variables and constraints for the minimization of metabolic
    adjustment (MOMA) to the model.

    Parameters
    ----------
    model : cobra.Model
        The model to add MOMA constraints and objective to.
    solution : cobra.Solution, optional
        A previous solution to use as a reference. If no solution is given,
        one will be computed using pFBA (default None).
    linear : bool, optional
        Whether to use the linear MOMA formulation or not (default True).

    Notes
    -----
    In the original MOMA [1]_ specification, one looks for the flux
    distribution of the deletion (v^d) closest to the fluxes without the
    deletion (v).
    In math this means:

    minimize: \sum_i (v^d_i - v_i)^2
    s.t.    : Sv^d = 0
              lb_i \le v^d_i \le ub_i

    Here, we use a variable transformation v^t := v^d_i - v_i. Substituting
    and using the fact that Sv = 0 gives:

    minimize: \sum_i (v^t_i)^2
    s.t.    : Sv^d = 0
              v^t = v^d_i - v_i
              lb_i \le v^d_i \le ub_i

    So, basically we just re-center the flux space at the old solution and
    then find the flux distribution closest to the new zero (center). This
    is the same strategy as used in cameo.

    In the case of linear MOMA [2]_, we instead minimize \sum_i abs(v^t_i).
    The linear MOMA is typically significantly faster. Also, quadratic MOMA
    tends to give flux distributions in which all fluxes deviate from the
    reference fluxes a little bit whereas linear MOMA tends to give flux
    distributions where the majority of fluxes are the same reference with
    few fluxes deviating a lot (typical effect of L2 norm vs L1 norm).

    The former objective function is saved in the optlang solver interface as
    ``"moma_old_objective"`` and this can be used to immediately extract the
    value of the former objective after MOMA optimization.

    See Also
    --------
    pfba : parsimonious FBA

    References
    ----------
    .. [1] Segrè, Daniel, Dennis Vitkup, and George M. Church. “Analysis of
           Optimality in Natural and Perturbed Metabolic Networks.”
           Proceedings of the National Academy of Sciences 99, no. 23
           (November 12, 2002): 15112. https://doi.org/10.1073/pnas.232349399.
    .. [2] Becker, Scott A, Adam M Feist, Monica L Mo, Gregory Hannum,
           Bernhard Ø Palsson, and Markus J Herrgard. “Quantitative
           Prediction of Cellular Metabolism with Constraint-Based Models:
           The COBRA Toolbox.” Nature Protocols 2 (March 29, 2007): 727.
    ""
    prstr="Creating Model to minimize difference to target fluxes: "
    if linear:
       prstr+="Lineal Mode"
    else:
       prstr+="Quadratic Mode"
    print(prstr)        
        
    #CF added
    if isinstance(target_fluxes, pd.Series):
       target_fluxes = target_fluxes.to_dict()
    if(copy_model):model=model.copy()
    #Scale weights 
    scaling_factor=max([abs(target_fluxes_weight[x])for x in target_fluxes_weight])
    
    #End CF added
    if "old_objective" in model.solver.variables:
        raise ValueError("The model is already adjusted for Flux Difference")
    # Fall back to default QP solver if current one has no QP capability
    if not linear and sutil.interface_to_str(model.problem) not in sutil.qp_solvers:
        model.solver = sutil.choose_solver(model, qp=True)
    #if solution is None:
    #    solution = pfba(model)
    prob = model.problem
    v = prob.Variable("old_objective")
    c = prob.Constraint(
        model.solver.objective.expression - v,
        lb=0.0,
        ub=0.0,
        name="old_objective_constraint",
    )
    to_add = [v, c]
    model.objective = prob.Objective(optlang_Zero, direction="min", sloppy=True)
    obj_vars = []
    linear_weights_dict = {}
    for r_id in target_fluxes:
        r=model.reactions.get_by_id(r_id)
        flux = target_fluxes[r_id]
        weight=target_fluxes_weight[r_id]/scaling_factor
        if linear:
            components = sutil.add_absolute_expression(
                model,
                r.flux_expression,
                name="target_flux_dist_" + r.id,
                difference=flux,
                add=False,
            )
            to_add.extend(components)
            obj_vars.append(components.variable)
            linear_weights_dict[components.variable] = weight
        else:
            dist = prob.Variable("target_flux_dist_" + r.id)
            const = prob.Constraint(
                r.flux_expression - dist,
                lb=flux,
                ub=flux,
                name="target_flux_constraint_" + r.id,
            )
            to_add.extend([dist, const])
            obj_vars.append(weight * dist**2)
    model.add_cons_vars(to_add)
    if linear:
        model.objective.set_linear_coefficients(linear_weights_dict)
    else:
        model.objective = prob.Objective(add(obj_vars), direction="min", sloppy=True)
    print(model.objective.expression)      
    return model 

"""
