import pandas as pd
from cobra.util import solver as sutil
from optlang.symbolics import Zero as optlang_Zero, add #, Pow

from cobra.exceptions import OptimizationError


try:
  from cplex import Cplex
except:
  pass  
        


import logging
logger = logging.getLogger(__name__)

#TODO Implement a direct cplex interface compatible with this as it might be faster

def add_quadratic_difference_minimization(model, target_fluxes={},target_fluxes_weight={},scaling_factor=None,expanded_reaction_mapping_dict={},copy_model=True,verbose=False,):
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
    if scaling_factor is None:
       scaling_factor = max(abs(target_fluxes_weight[x]) for x in target_fluxes_weight)
    if verbose:
        print(f"Scaling factor for weights: {scaling_factor}")
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
    #obj_vars = []
    optimized_fluxes_list=[]
    quad_terms = []
    lin_terms = []
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
        else:
          optimized_fluxes_list.append(flux_id)  #Keep track to the reactions that are actually optimized
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
        #obj_vars.append(weight * net_flux**2 - 2*weight*target_flux_value*net_flux)
        quad_terms.append(weight * net_flux**2)
        lin_terms.append(-2 * weight * target_flux_value * net_flux)
        if(verbose):
           print("######################## "+flux_id)
           print("Adding flux difference minimization for reaction: "+flux_id+" with target flux: "+str(target_flux_value)+" and weight: "+str(weight) )
           print("Net flux expression: "+str(net_flux_expression) )
           print("Objective term: "+str(weight * net_flux**2 - 2*weight*target_flux_value*net_flux) )
    model.add_cons_vars(to_add)
    #model.objective = prob.Objective(add(obj_vars), direction="min", sloppy=True)
    model.objective = prob.Objective( add(quad_terms) + add(lin_terms), direction="min", sloppy=True)
    #print(model.objective.expression)      
    return model,optimized_fluxes_list 


def update_quadratic_objective_coefficients(
    model,
    target_fluxes,
    target_fluxes_weight,
    scaling_factor=None,
    verbose=False
):
    if isinstance(target_fluxes, pd.Series):
        target_fluxes = target_fluxes.to_dict()
    if isinstance(target_fluxes_weight, pd.Series):
        target_fluxes_weight = target_fluxes_weight.to_dict()
    if scaling_factor is None:
       # Scale weights
       scaling_factor = max(abs(target_fluxes_weight[x]) for x in target_fluxes_weight)
    if verbose:
        print(f"Scaling factor for weights: {scaling_factor}")
    prob = model.problem
    quad_terms = []
    lin_terms = []
    for flux_id, target in target_fluxes.items():
        weight = target_fluxes_weight[flux_id] / scaling_factor
        v = model.solver.variables["net_flux_" + flux_id]
        quad_terms.append(weight * v * v)
        lin_terms.append(-2 * weight * target * v)
        if verbose:
            print(f"Updated flux {flux_id}: target={target:.4g}, weight={weight:.4g}")
    model.objective = prob.Objective(
        add(quad_terms) + add(lin_terms),
        direction="min",
        sloppy=True,
    )
    if verbose:
        print(f"Objective rebuilt for {len(target_fluxes)} fluxes")
    return model
    
def update_quadratic_objective_coefficients_cplex(
    model,
    target_fluxes,
    target_fluxes_weight,
    scaling_factor=None,
    verbose=False
):
    """
    Update quadratic objective coefficients using direct CPLEX interface.
    Orders of magnitude faster than optlang symbolic version.
    
    Parameters
    ----------
    model : cobra.Model
        Model with existing net_flux variables (must use CPLEX solver)
    target_fluxes : dict
        New target flux values
    target_fluxes_weight : dict
        New weights
    verbose : bool
        Print progress
        
    Returns
    -------
    cobra.Model
        Model with updated objective
    """
    # Verify CPLEX solver
    #try:
    #    from cplex import Cplex
    #except ImportError:
    #    raise ImportError("CPLEX not available")
    
    current_solver = sutil.interface_to_str(model.problem)
    if 'cplex' not in current_solver.lower():
        raise ValueError(
            f"Model is using {current_solver} solver, but this function requires CPLEX.\n"
            f"Use model.solver = 'cplex' before calling this function, or use the "
            f"non-CPLEX version: update_quadratic_objective_coefficients()"
        )
    
    if isinstance(target_fluxes, pd.Series):
        target_fluxes = target_fluxes.to_dict()
    
    if isinstance(target_fluxes_weight, pd.Series):
        target_fluxes_weight = target_fluxes_weight.to_dict()
    
    # Scale weights
    if scaling_factor is None:
       # Scale weights
       scaling_factor = max(abs(target_fluxes_weight[x]) for x in target_fluxes_weight)
    if verbose:
        print(f"Scaling factor for weights: {scaling_factor}")    
    #Reset Objective
    #model.objective = model.problem.Objective(optlang_Zero, direction='min')
    # Access CPLEX problem directly
    lp = model.solver.problem
    
    # Verify it's a CPLEX object
    if not isinstance(lp, Cplex):
        raise RuntimeError(f"Expected CPLEX problem but got {type(lp)}")
    
    # Build coefficient lists
    quadratic_updates = []
    linear_updates = []
    
    # Collect all net_flux variable names
    var_names = []
    flux_ids = []
    for flux_id in target_fluxes.keys():
        var_name = "net_flux_" + flux_id
        if var_name in model.solver.variables:
            var_names.append(var_name)
            flux_ids.append(flux_id)
        else:
            if verbose:
                print(f"Warning: Variable {var_name} not found, skipping")
    
    # Get all indices at once (batch operation)
    var_indices = lp.variables.get_indices(var_names)
    
    # Build update lists
    for var_index, flux_id in zip(var_indices, flux_ids):
        target_flux_value = target_fluxes[flux_id]
        weight = target_fluxes_weight[flux_id] / scaling_factor
        
        # Quadratic coefficient: weight
        quadratic_updates.append((var_index, var_index, weight*2)) #Weight needs to be multiplied by 2 because Cplex seems to divide the quadratic terms by 2. This is done by default int optlang
        
        # Linear coefficient: -2 * weight * target
        linear_updates.append((var_index, -2 * weight * target_flux_value))
        
        if verbose:# and len(quadratic_updates) <= 10:
            print(f"  {flux_id}: target={target_flux_value:.4g}, weight={weight:.4g}")
        #elif verbose and len(quadratic_updates) == 11:
        #    print(f"  ... (suppressing output for remaining {len(var_names)-10} fluxes)")
    
    # Apply updates using CPLEX batch operations
    if quadratic_updates:
        for i, j, coef in quadratic_updates:
            lp.objective.set_quadratic_coefficients(i, j, coef)
    
    if linear_updates:
        lp.objective.set_linear(linear_updates)
    
    if verbose:
        print(f"Objective rebuilt for {len(target_fluxes)} fluxes using CPLEX interface")
    
    return model


def qMTA_get_optimization_target_fluxes_and_weights(
    reference_fluxes,
    model=None,
    target_net_flux_fold_change=None,
    target_measured_flux=None,
    fold_change_weight=None,
    target_measured_flux_weight=None,
    fluxes_to_omit=None,
    target_fold_change_base_weight=0.5,
    target_measured_flux_base_weight=0.5,
    unchanged_reaction_base_weight=0.5,
    fold_change_is_Log2=True,
    reference_fluxes_are_net_fluxes=False,
    max_Log2FC_change=None, #max_fold_change
    scale_target_fold_change_weight_by_target_flux_difference=True, #normalize_by_scale_genes
    scale_measured_target_flux_weight_by_target_flux_difference=False, #normalize_by_scale_mets
    scale_unchanged_reactions_by_reference_flux=True, #normalize_by_scale_unchanged_reactions
    min_factor_for_scaling=1e-6, #min_flux4weight
    min_flux_to_apply_fold_change=1e-6, #min_flux_fold_change
    precision=None,
    expanded_reaction_mapping_dict=None,
    verbose=True
    ):
    """
    Get target net fluxes and weights for qMTA optimization.
    
    Parameters
    ----------
    reference_fluxes : dict or pd.Series
        Reference flux values
    model : cobra.Model, optional
        The metabolic model. Only required if expanded_reaction_mapping_dict or net fluxes are not provided. 
    target_net_flux_fold_change : dict or pd.Series, optional
        Fold changes to apply to reference fluxes
    target_measured_flux : dict or pd.Series, optional
        Direct target flux measurements
    fold_change_weight : dict or pd.Series, optional
        Weights for fold change targets
    target_measured_flux_weight : dict or pd.Series, optional
        Weights for measured flux targets
    fluxes_to_omit : list, optional
        Flux IDs to exclude from optimization
    target_fold_change_base_weight : float
        Base weight for fold change targets (default 0.5)
    target_measured_flux_base_weight : float
        Base weight for measured flux targets (default 0.5)
    unchanged_reaction_base_weight : float
        Base weight for unchanged reactions (default 0.5)
    fold_change_is_Log2 : bool
        Whether fold changes are in log2 scale (default True)
    reference_fluxes_are_net_fluxes : bool
        Whether reference fluxes are already net fluxes (default False)
    max_Log2FC_change : float
        Maximum allowed log2 fold change (default 99999999)
    scale_target_fold_change_weight_by_target_flux_difference : bool
        Scale fold change weights by flux difference (default True)
    scale_measured_target_flux_weight_by_target_flux_difference : bool
        Scale measured flux weights by flux difference (default False)
    scale_unchanged_reactions_by_reference_flux : bool
        Scale unchanged reaction weights by reference flux (default True)
    min_factor_for_scaling : float
        Minimum factor for weight scaling (default 1e-6)
    min_flux_to_apply_fold_change : float
        Minimum reference flux to apply fold change (default 1e-6)
    precision : int, optional
        Precision for rounding weights (applies to weights only, not target fluxes)
    expanded_reaction_mapping_dict : dict, optional
        Mapping of net fluxes to reactions
    verbose : bool
        Print detailed progress (default True)
        
    Returns
    -------
    tuple of (dict, dict)
        optimization_target_fluxes : Target flux values
        optimization_target_fluxes_weights : Weights for each target flux
    """
    # Handle mutable default arguments
    if target_net_flux_fold_change is None:
        target_net_flux_fold_change = {}
    if target_measured_flux is None:
        target_measured_flux = {}
    if fold_change_weight is None:
        fold_change_weight = {}
    if target_measured_flux_weight is None:
        target_measured_flux_weight = {}
    if fluxes_to_omit is None:
        fluxes_to_omit = []
    if precision is not None or verbose:
      from .cobra_functions import round_sig
    # Convert inputs to dicts if they are not already
    if isinstance(reference_fluxes, pd.Series):
        reference_fluxes = reference_fluxes.to_dict()
    
    if isinstance(target_net_flux_fold_change, pd.Series):
        target_net_flux_fold_change = target_net_flux_fold_change.to_dict()
    
    if isinstance(fold_change_weight, pd.Series):
        fold_change_weight = fold_change_weight.to_dict()
    
    if isinstance(target_measured_flux, pd.Series):
        target_measured_flux = target_measured_flux.to_dict()
    
    if isinstance(target_measured_flux_weight, pd.Series):
        target_measured_flux_weight = target_measured_flux_weight.to_dict()
    
    #Raise exception if reference fluxes is empty
    if len(reference_fluxes) == 0:
        raise Exception("Error: reference_fluxes is empty")
    # Auto-generate mapping if not provided
    if expanded_reaction_mapping_dict is None:
        print("Building default expanded_reaction_mapping_dict with reverse_reaction_pattern=_REV, isoenzyme_reaction_pattern=_EXP_\\d+ and patterns_to_omit=[^usage_prot_]")
        from .ec_base import get_ec_expanded_reaction_mapping
        expanded_reaction_mapping_dict, _ = get_ec_expanded_reaction_mapping(
            model, 
            reverse_reaction_pattern="_REV", 
            isoenzyme_reaction_pattern="_EXP_\\d+",
            patterns_to_omit=["^usage_prot_"],
            verbose=False
        )
    
    # Convert reference fluxes to net fluxes if needed
    if reference_fluxes_are_net_fluxes == False:
        from .ec_base import get_net_fluxes_from_ec_model
        
        print("Converting reference fluxes to net fluxes")
        net_fluxes = get_net_fluxes_from_ec_model(
            model,
            fluxes=reference_fluxes,
            output_flux_breakdown=False,
            ec_expanded_reaction_mapping_dict=expanded_reaction_mapping_dict
        )  
        reference_fluxes = net_fluxes['net_flux'].to_dict()
    
    # Convert log2 fold changes to linear fold changes
    if fold_change_is_Log2:
        target_net_flux_fold_change = {
            rid: pow(2, fc) for rid, fc in target_net_flux_fold_change.items()
        }
    
    # Get Max and Min fold changes thresholds
    if max_Log2FC_change is not None:
       max_fc = pow(2, max_Log2FC_change)
       min_fc = pow(2, -max_Log2FC_change)
    
    # Get fold change fluxes
    fold_change_fluxes = set(target_net_flux_fold_change.keys()) - set(fluxes_to_omit)
    
    # Get target net fluxes
    measured_target_fluxes = set(target_measured_flux.keys()) - set(fluxes_to_omit)
    
    # Raise an exception if there are overlapping reactions
    overlapping_reactions = fold_change_fluxes.intersection(measured_target_fluxes)
    if len(overlapping_reactions) > 0:
        raise Exception(
            "Error: The following reactions are present in both fold change and measured target fluxes: {}".format(
                ", ".join(overlapping_reactions)
            )
        )
    
    if len(fold_change_fluxes) == 0 and len(measured_target_fluxes) == 0:
        logger.warning("No fold change or measured target fluxes to process")
    
    if verbose:
        if len(measured_target_fluxes) > 0:
            print("Number of measured target fluxes to process: {}".format(len(measured_target_fluxes)))
        if len(fold_change_fluxes) > 0:
            print("Number of fold change fluxes to process: {}".format(len(fold_change_fluxes)))
    
    optimization_target_fluxes = {}
    optimization_target_fluxes_weights = {}
    fluxes_excluded_from_fold_change = set()
    counter_fold_change = 0
    counter_measured = 0
    counter_unchanged = 0
    
    # Build target fluxes and weights for fold change reactions
    for net_flux_id in fold_change_fluxes:
        vref = reference_fluxes.get(net_flux_id)
        if vref is None:
            if verbose:
                logger.warning("Reaction {} not found in reference fluxes, skipping".format(net_flux_id))
            continue
        
        if abs(vref) < min_flux_to_apply_fold_change:
            # Add them to excluded fluxes
            fluxes_excluded_from_fold_change.add(net_flux_id)
            if verbose:
                print("Reaction {} with reference flux {} below min_flux_to_apply_fold_change {}. Adding it to excluded from fold change reactions.".format(
                    net_flux_id, round_sig(vref, 4), min_flux_to_apply_fold_change
                ))
            continue
        
        fold_change = target_net_flux_fold_change[net_flux_id]
        if max_Log2FC_change is not None:
           fold_change = min(max_fc, max(min_fc, fold_change))
        target_flux = vref * fold_change
        optimization_target_fluxes[net_flux_id] = target_flux
        
        # Get weight
        weight = target_fold_change_base_weight * fold_change_weight.get(net_flux_id, 1)
        if scale_target_fold_change_weight_by_target_flux_difference:
            flux_diff_squared = max(pow(vref - target_flux, 2), min_factor_for_scaling)
            weight /= flux_diff_squared
        
        if precision is not None:
            weight = round_sig(weight, precision)
        
        optimization_target_fluxes_weights[net_flux_id] = weight
        counter_fold_change += 1
        
        if verbose:
            # Print vref, fold change, target flux, weight
            print("rid: {}, vref: {}, fold change: {}, target flux: {}, weight: {}".format(
                net_flux_id, round_sig(vref, 4), round_sig(fold_change, 4), 
                round_sig(target_flux, 4), weight
            ))
    
    # Build target fluxes and weights for measured target flux reactions
    for net_flux_id in measured_target_fluxes:
        vref = reference_fluxes.get(net_flux_id)
        if vref is None:
            if verbose:
                logger.warning("Reaction {} not found in reference fluxes, skipping".format(net_flux_id))
            continue
        
        # We will not filter by min_flux_to_apply_fold_change here, as these are measured fluxes         
        target_flux = target_measured_flux[net_flux_id]
        optimization_target_fluxes[net_flux_id] = target_flux
        
        # Get weight
        weight = target_measured_flux_base_weight * target_measured_flux_weight.get(net_flux_id, 1)
        if scale_measured_target_flux_weight_by_target_flux_difference:
            flux_diff_squared = max(pow(vref - target_flux, 2), min_factor_for_scaling)
            weight /= flux_diff_squared
        
        if precision is not None:
            weight = round_sig(weight, precision)
        
        optimization_target_fluxes_weights[net_flux_id] = weight
        counter_measured += 1
        
        if verbose:
            # Print vref, target flux, weight
            print("Measured reaction: rid: {}, vref: {}, target flux: {}, weight: {}".format(
                net_flux_id, round_sig(vref, 4), round_sig(target_flux, 4), weight
            ))
    
    # Get unchanged reactions
    unchanged_fluxes = (
        set(reference_fluxes.keys()) - 
        fold_change_fluxes - 
        measured_target_fluxes - 
        fluxes_excluded_from_fold_change - 
        set(fluxes_to_omit)
    )
    
    if verbose:
        print("Number of unchanged fluxes to process: {}".format(len(unchanged_fluxes)))
    # Build target fluxes and weights for unchanged reactions
    for net_flux_id in unchanged_fluxes:
        vref = reference_fluxes.get(net_flux_id)
        if vref is None:
            if verbose:
                logger.warning("Reaction {} not found in reference fluxes, skipping".format(net_flux_id))
            continue
        
        target_flux = vref
        optimization_target_fluxes[net_flux_id] = target_flux
        
        # Get weight
        weight = unchanged_reaction_base_weight
        if scale_unchanged_reactions_by_reference_flux:
            weight /= max(abs(vref), min_factor_for_scaling) 
        
        if precision is not None:
            weight = round_sig(weight, precision)
        
        optimization_target_fluxes_weights[net_flux_id] = weight
        counter_unchanged += 1
        
        if verbose:
            print("Unchanged reaction: rid: {}, vref: {}, weight: {}".format(
                net_flux_id, round_sig(vref, 4), weight
            ))
    
    if verbose:
        # Print how many passed filters
        print("Total fold change reactions after processing: {}".format(counter_fold_change))
        print("Total measured target flux reactions after processing: {}".format(counter_measured))
        print("Total unchanged reactions after processing: {}".format(counter_unchanged))
    
    return optimization_target_fluxes, optimization_target_fluxes_weights

    
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




def solve_qp_model(qp_model, params={}, cold_start=False):
    """
    Configure optlang / CPLEX solver settings on `qp_model` from a dict and
    then optimize, trying multiple QP methods in order until one succeeds.

    Any parameter whose key is absent from `params` is simply never set,
    i.e. it's left at whatever the solver's current/default value is.

    Parameters
    ----------
    qp_model : cobra.Model
        Model to configure and solve. Must already have `.solver` set.
    params : dict, optional
        Recognized keys (all optional):

        General (optlang) parameters
            solver_qp_methods              list[str] - QP methods to try,
                                            in order. If omitted, optimize()
                                            is called once without touching
                                            qp_method at all.
            solver_tolerances_feasibility  float
            solver_tolerances_optimality   float
            solver_verbosity               int
            solver_presolve                True / False / "auto"

        CPLEX-only parameters (only applied if the current solver is CPLEX;
        optlang has no bindings for these, so they're set on the raw
        `qp_model.solver.problem`):
            cplex_n_threads
            cplex_emphasis_numerical
            cplex_scaling
            cplex_network_netfind
            cplex_network_pricing
            cplex_network_tol_feasibility
            cplex_network_tol_optimality
            cplex_barrier_crossover
            cplex_barrier_ordering
            cplex_barrier_startalg
            cplex_barrier_convergetol
    cold_start : bool, default False
        CPLEX-only. If True, disables the advanced (warm) start basis
        (parameters.advance = 0), forcing a solve from scratch. If False,
        lets CPLEX reuse a previous basis (parameters.advance = 1).
        
    Example Params:
        {
        "solver_qp_methods": ["network", "barrier", "primal"],
        "solver_tolerances_feasibility": 1e-7,
        "solver_tolerances_optimality": 1e-6,
        "solver_verbosity": 2,
        "solver_presolve": False,
        "cplex_n_threads": 1,
        "cplex_emphasis_numerical": 0,
        "cplex_scaling": 0,
        "cplex_network_netfind": 1,
        "cplex_network_pricing": 0,
        "cplex_network_tol_feasibility": 1e-7,
        "cplex_network_tol_optimality": 1e-6,
        "cplex_barrier_crossover": 0,
        "cplex_barrier_ordering": 0,
        "cplex_barrier_startalg": 1,
        "cplex_barrier_convergetol": 1e-6,
    }


    Returns
    -------
    sol : solution object
    stats : str or Exception. "optimal" if a solution was found, otherwise the last status /
        OptimizationError encountered.

        
    """
    # ---- Validate that every key in params is a recognized option ----
    known_keys = {
        "solver_qp_methods",
        "solver_tolerances_feasibility",
        "solver_tolerances_optimality",
        "solver_verbosity",
        "solver_presolve",
        "cplex_n_threads",
        "cplex_emphasis_numerical",
        "cplex_scaling",
        "cplex_network_netfind",
        "cplex_network_pricing",
        "cplex_network_tol_feasibility",
        "cplex_network_tol_optimality",
        "cplex_barrier_crossover",
        "cplex_barrier_ordering",
        "cplex_barrier_startalg",
        "cplex_barrier_convergetol",
        "cplex_timelimit"
    }
    unknown_keys = set(params) - known_keys
    if unknown_keys:
        raise ValueError(
            "Unrecognized solver parameter(s): " + ", ".join(sorted(unknown_keys))
        )

    # ---- General optlang configuration ----
    if "solver_tolerances_feasibility" in params:
        qp_model.solver.configuration.tolerances.feasibility = params["solver_tolerances_feasibility"]
    if "solver_tolerances_optimality" in params:
        qp_model.solver.configuration.tolerances.optimality = params["solver_tolerances_optimality"]
    if "solver_verbosity" in params:
        qp_model.solver.configuration.verbosity = params["solver_verbosity"]
    if "solver_presolve" in params:
        qp_model.solver.configuration.presolve = params["solver_presolve"]

    # ---- Identify current solver ----
    current_solver = sutil.interface_to_str(qp_model.problem)
    print("Current Solver is " + current_solver)
    solver_is_cplex = "cplex" in current_solver.lower()

    if solver_is_cplex:
        if "cplex_n_threads" in params:
            qp_model.solver.problem.parameters.threads.set(params["cplex_n_threads"])
        if "cplex_emphasis_numerical" in params:
            qp_model.solver.problem.parameters.emphasis.numerical.set(params["cplex_emphasis_numerical"])
        if "cplex_scaling" in params:
            qp_model.solver.problem.parameters.read.scale.set(params["cplex_scaling"])
        if "cplex_network_netfind" in params:
            qp_model.solver.problem.parameters.network.netfind.set(params["cplex_network_netfind"])
        if "cplex_network_pricing" in params:
            qp_model.solver.problem.parameters.network.pricing.set(params["cplex_network_pricing"])
        if "cplex_network_tol_feasibility" in params:
            qp_model.solver.problem.parameters.network.tolerances.feasibility.set(params["cplex_network_tol_feasibility"])
        if "cplex_network_tol_optimality" in params:
            qp_model.solver.problem.parameters.network.tolerances.optimality.set(params["cplex_network_tol_optimality"])
        if "cplex_barrier_crossover" in params:
            qp_model.solver.problem.parameters.barrier.crossover.set(params["cplex_barrier_crossover"])
        if "cplex_barrier_ordering" in params:
            qp_model.solver.problem.parameters.barrier.ordering.set(params["cplex_barrier_ordering"])
        if "cplex_barrier_startalg" in params:
            qp_model.solver.problem.parameters.barrier.startalg.set(params["cplex_barrier_startalg"])
        if "cplex_barrier_convergetol" in params:
            qp_model.solver.problem.parameters.barrier.convergetol.set(params["cplex_barrier_convergetol"])
        if "cplex_timelimit" in params:
            qp_model.solver.problem.parameters.timelimit.set(params["cplex_timelimit"])
        # Cold start vs warm start (CPLEX "advance" basis parameter)
        qp_model.solver.problem.parameters.advance.set(0 if cold_start else 1)

    # ---- Try QP methods in order until an optimal solution is found ----
    # If the key is missing, don't touch qp_method at all — just optimize once.
    qp_methods = params["solver_qp_methods"] if "solver_qp_methods" in params else [None]

    status = None
    objective_value = None
    for solver_qp_method in qp_methods:
        if solver_qp_method is not None and qp_model.solver.configuration.qp_method != solver_qp_method:
            qp_model.solver.configuration.qp_method = solver_qp_method

        try:
            sol = qp_model.optimize()
            status = sol.status
            objective_value = sol.objective_value
        except OptimizationError as error:
            status = error
            objective_value =None
            sol=None

        print("\t Objective Value:" + str(status) + " " + str(objective_value) + " with " + str(solver_qp_method))

        if status == "optimal":
            break

    return sol, status

#Example Params
#    {
#        "solver_qp_methods": ["network", "barrier", "primal"],
#        "solver_tolerances_feasibility": 1e-7,
#        "solver_tolerances_optimality": 1e-6,
#        "solver_verbosity": 2,
#        "solver_presolve": False,
#        "cplex_n_threads": 1,
#        "cplex_emphasis_numerical": 0,
#        "cplex_scaling": 0,
#        "cplex_network_netfind": 1,
#        "cplex_network_pricing": 0,
#        "cplex_network_tol_feasibility": 1e-7,
#        "cplex_network_tol_optimality": 1e-6,
#        "cplex_barrier_crossover": 0,
#        "cplex_barrier_ordering": 0,
#        "cplex_barrier_startalg": 1,
#        "cplex_barrier_convergetol": 1e-6,
#        "cplex_timelimit":1000
#    }
