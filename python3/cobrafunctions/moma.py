import cobra
import pandas as pd
from cobra.util import solver as sutil
from optlang.symbolics import Zero, add



def simulate_reaction_ko_moma(model,reference_fluxes={},reactions_to_ko_dict={},ko_factor=0.5,ko_reference_flux={},quadratic=True):
   results_dict={}
   """for x in model.reactions:
       if x.id not in reference_fluxes:
          reference_fluxes[x.id]=0"""
   if isinstance(reference_fluxes, pd.Series):
       reference_fluxes = reference_fluxes.to_dict()
   if isinstance(ko_reference_flux, pd.Series):
       ko_reference_flux = ko_reference_flux.to_dict()
   constrained_model=add_moma(model,reference_fluxes=reference_fluxes,copy_model=True,linear= not quadratic)
   for nr,gene in enumerate(reactions_to_ko_dict):
        reaction_list=reactions_to_ko_dict[gene]
        with constrained_model:  #Outised this context changes are reverted automatically
          #print constrained_model.optimize().fluxes["BIOMASS"]#REMOVE ME
          for reaction_id in reaction_list:
              reaction=constrained_model.reactions.get_by_id(reaction_id)
              if reaction_id in ko_reference_flux:#To define a different reference flux
                 vref=ko_reference_flux[reaction_id] 
              else:    
                 vref=reference_fluxes[reaction_id]
              bound=vref*ko_factor#round_sig(vref*ko_factor,3) 
              reaction.lower_bound=max(reaction.lower_bound,-abs(bound))
              reaction.upper_bound=min(reaction.upper_bound,abs(bound))
          solution=constrained_model.optimize()
          results_dict[gene]={x.id:solution.fluxes[x.id] for x in constrained_model.reactions }
          output_str=str(nr)+"/"+str(len(reactions_to_ko_dict))
          print(output_str)
   return results_dict



def add_moma(model, reference_fluxes={}, linear=False,copy_model=True,):
    r"""
    ###Adapated from cobrapy 
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
    """
    prstr="Creating MOMA Model: "
    if linear:
       prstr+="Lineal Mode"
    else:
       prstr+="Quadratic Mode"
    print(prstr)        
    #CF added
    if isinstance(reference_fluxes, pd.Series):
       reference_fluxes = reference_fluxes.to_dict()
    if(copy_model):model=model.copy()
    #End CF added
    if "moma_old_objective" in model.solver.variables:
        raise ValueError("The model is already adjusted for MOMA.")
    # Fall back to default QP solver if current one has no QP capability
    if not linear and sutil.interface_to_str(model.problem) not in sutil.qp_solvers:
        model.solver = sutil.choose_solver(model, qp=True)
    #if solution is None:
    #    solution = pfba(model)
    prob = model.problem
    v = prob.Variable("moma_old_objective")
    c = prob.Constraint(
        model.solver.objective.expression - v,
        lb=0.0,
        ub=0.0,
        name="moma_old_objective_constraint",
    )
    to_add = [v, c]
    model.objective = prob.Objective(Zero, direction="min", sloppy=True)
    obj_vars = []
    for r in model.reactions:
        flux = reference_fluxes[r.id]
        if linear:
            components = sutil.add_absolute_expression(
                model,
                r.flux_expression,
                name="moma_dist_" + r.id,
                difference=flux,
                add=False,
            )
            to_add.extend(components)
            obj_vars.append(components.variable)
        else:
            dist = prob.Variable("moma_dist_" + r.id)
            const = prob.Constraint(
                r.flux_expression - dist,
                lb=flux,
                ub=flux,
                name="moma_constraint_" + r.id,
            )
            to_add.extend([dist, const])
            obj_vars.append(dist**2)
    model.add_cons_vars(to_add)
    if linear:
        model.objective.set_linear_coefficients({v: 1.0 for v in obj_vars})
    else:
        model.objective = prob.Objective(add(obj_vars), direction="min", sloppy=True)
    return model 



