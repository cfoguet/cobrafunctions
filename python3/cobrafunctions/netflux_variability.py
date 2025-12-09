"""Flux variability analysis for net fluxes defined by forward/reverse reaction mappings."""
""" Adapted from the variability cobrapy function"""

import logging
from typing import TYPE_CHECKING, List, Optional, Dict, Union
from warnings import warn

import numpy as np
import pandas as pd
from optlang.symbolics import Zero, add as optlang_add

from cobra.core import Configuration, get_solution
from cobra.util import ProcessPool
from cobra.util import solver as sutil
from .ec import get_ec_expanded_reaction_mapping


if TYPE_CHECKING:
    from cobra import Gene, Model, Reaction


if TYPE_CHECKING:
    from cobra import Gene, Model, Reaction


logger = logging.getLogger(__name__)
configuration = Configuration()


def _init_worker_net_flux(
    model: "Model", 
    sense: str,
    expanded_reaction_mapping_dict: Dict
) -> None:
    """Initialize a global model object for multiprocessing with net flux mapping.

    Parameters
    ----------
    model: cobra.Model
        The model to operate on.
    sense: {"max", "min"}
        Whether to maximise or minimise objective.
    expanded_reaction_mapping_dict: dict
        Mapping of flux IDs to forward/reverse reactions.
    """
    global _model
    global _expanded_mapping
    _model = model
    _model.solver.objective.direction = sense
    _expanded_mapping = expanded_reaction_mapping_dict


def _fva_step_net_flux(flux_id: str) -> tuple:
    """Take a step for calculating FVA on a net flux.

    Parameters
    ----------
    flux_id: str
        The ID of the flux (from expanded_reaction_mapping_dict keys).

    Returns
    -------
    tuple of (str, float)
        The flux ID with the net flux value.
    """
    global _model
    global _expanded_mapping
    
    forward_reactions = _expanded_mapping[flux_id]["forward_reactions"]
    reverse_reactions = _expanded_mapping[flux_id]["reverse_reactions"]
    
    # Build the net flux expression
    net_flux_expression = Zero
    
    # Forward/Standard reactions
    for forward_rid in forward_reactions:
        r = _model.reactions.get_by_id(forward_rid)
        if r.forward_variable.ub > 0:
            net_flux_expression += r.forward_variable
        if r.reverse_variable.ub > 0:
            if len(reverse_reactions) > 0:
                raise Exception(
                    f"{forward_rid}: If using forward/reverse reactions, "
                    "forward reaction should not have negative lower bound"
                )
            net_flux_expression -= r.reverse_variable
    
    # Reverse reactions if applicable
    for reverse_rid in reverse_reactions:
        r_rev = _model.reactions.get_by_id(reverse_rid)
        if r_rev.forward_variable.ub > 0:
            net_flux_expression -= r_rev.forward_variable
        if r_rev.reverse_variable.ub > 0:
            raise Exception(
                f"{reverse_rid}: Reverse reactions should not have reverse fluxes allowed"
            )
    
    if net_flux_expression is Zero:
        logger.warning(f"Skipping flux {flux_id} as no flux is allowed in either direction")
        return flux_id, float("nan")
    
    # Set objective to the net flux expression
    _model.solver.objective.set_linear_coefficients(
        {var: coef for var, coef in net_flux_expression.as_coefficients_dict().items()}
    )
    
    _model.slim_optimize()
    sutil.check_solver_status(_model.solver.status)
    
    value = _model.solver.objective.value
    
    # Handle infeasible case
    if value is None:
        value = float("nan")
        logger.warning(
            f"Could not get flux for {flux_id}, setting it to NaN. "
            "This is usually due to numerical instability."
        )
    
    # Reset coefficients
    _model.solver.objective.set_linear_coefficients(
        {var: 0 for var in net_flux_expression.as_coefficients_dict().keys()}
    )
    
    return flux_id, value


def flux_variability_analysis_net_flux(
    model: "Model",
    expanded_reaction_mapping_dict: Optional[Dict[str, Dict[str, List[str]]]] = None,
    flux_list: Optional[List[str]] = None,
    fraction_of_optimum: float = 1.0,
    processes: Optional[int] = None,
) -> pd.DataFrame:
    """Determine the minimum and maximum net flux value for each flux.

    Note precision can be incrased with: model.solver.configuration.tolerances.feasibility = 1e-9
    
    This is similar to standard FVA but operates on net fluxes defined by
    forward and reverse reaction mappings, useful for models with split
    reversible reactions or when you want to analyze net flux across
    multiple reactions.

    Parameters
    ----------
    model : cobra.Model
        The model for which to run the analysis. It will *not* be modified.
    expanded_reaction_mapping_dict : dict
        A dictionary mapping flux IDs to forward and reverse reactions in the model.
        Format:
        {
           "flux_id1": {
               "forward_reactions": ["reaction_id1", "reaction_id2"],
               "reverse_reactions": ["reaction_id3"]
           },
           "flux_id2": {
               "forward_reactions": ["reaction_id4"],
               "reverse_reactions": []
           },
        }
    flux_list : list of str, optional
        The flux IDs for which to obtain min/max net fluxes. If None will use
        all fluxes in expanded_reaction_mapping_dict (default None).
    fraction_of_optimum : float, optional
        Must be <= 1.0. Requires that the objective value is at least the
        fraction times maximum objective value (default 1.0).
    processes : int, optional
        The number of parallel processes to run (default None).

    Returns
    -------
    pandas.DataFrame
        A data frame with flux identifiers as the index and two columns:
        - maximum: indicating the highest possible net flux
        - minimum: indicating the lowest possible net flux

    Examples
    --------
    >>> # For a standard model with split reversible reactions
    >>> mapping = {
    ...     "PGI_net": {
    ...         "forward_reactions": ["PGI_forward"],
    ...         "reverse_reactions": ["PGI_reverse"]
    ...     },
    ...     "PFK": {
    ...         "forward_reactions": ["PFK"],
    ...         "reverse_reactions": []
    ...     }
    ... }
    >>> result = flux_variability_analysis_net_flux(
    ...     model,
    ...     expanded_reaction_mapping_dict=mapping,
    ...     fraction_of_optimum=0.9
    ... )
    """
    # Check that the model doesn't already contain reporter metabolites
    existing_reporter_mets = [met.id for met in model.metabolites if met.id.startswith("reporter_met_")]
    
    if existing_reporter_mets:
        raise ValueError(
            f"Model already contains reporter metabolites for net fluxes: {existing_reporter_mets}. "
            "Should not run flux_variability_analysis_net_flux on a model with reporters. "
            "Use the original model without reporters."
        )
    
    # Auto-generate mapping if not provided
    if expanded_reaction_mapping_dict is None:
        print("Building default expanded_reaction_mapping_dict with reverse_reaction_pattern=_REV, isoenzyme_reaction_pattern=_EXP_\\d+")
        expanded_reaction_mapping_dict, _ = get_ec_expanded_reaction_mapping(
            model, 
            reverse_reaction_pattern="_REV", 
            isoenzyme_reaction_pattern="_EXP_\\d+"
        )
    
    if flux_list is None:
        flux_list = list(expanded_reaction_mapping_dict.keys())
    else:
        # Validate that all requested fluxes are in the mapping
        missing = set(flux_list) - set(expanded_reaction_mapping_dict.keys())
        if missing:
            raise ValueError(
                f"The following flux IDs are not in expanded_reaction_mapping_dict: {missing}"
            )

    if processes is None:
        processes = configuration.processes

    num_fluxes = len(flux_list)
    processes = min(processes, num_fluxes)

    fva_result = pd.DataFrame(
        {
            "minimum": np.zeros(num_fluxes, dtype=float),
            "maximum": np.zeros(num_fluxes, dtype=float),
        },
        index=flux_list,
    )
    
    prob = model.problem
    with model:
        # Safety check before setting up FVA
        model.slim_optimize(
            error_value=None,
            message="There is no optimal solution for the chosen objective!",
        )
        
        # Add the previous objective as a constraint
        if model.solver.objective.direction == "max":
            fva_old_objective = prob.Variable(
                "fva_old_objective",
                lb=fraction_of_optimum * model.solver.objective.value,
            )
        else:
            fva_old_objective = prob.Variable(
                "fva_old_objective",
                ub=fraction_of_optimum * model.solver.objective.value,
            )
        fva_old_obj_constraint = prob.Constraint(
            model.solver.objective.expression - fva_old_objective,
            lb=0,
            ub=0,
            name="fva_old_objective_constraint",
        )
        model.add_cons_vars([fva_old_objective, fva_old_obj_constraint])

        model.objective = Zero
        
        for what in ("minimum", "maximum"):
            if processes > 1:
                chunk_size = len(flux_list) // processes
                with ProcessPool(
                    processes,
                    initializer=_init_worker_net_flux,
                    initargs=(model, what[:3], expanded_reaction_mapping_dict),
                ) as pool:
                    for flux_id, value in pool.imap_unordered(
                        _fva_step_net_flux, flux_list, chunksize=chunk_size
                    ):
                        fva_result.at[flux_id, what] = value
            else:
                _init_worker_net_flux(
                    model, what[:3], expanded_reaction_mapping_dict
                )
                for flux_id, value in map(_fva_step_net_flux, flux_list):
                    fva_result.at[flux_id, what] = value

    return fva_result[["minimum", "maximum"]]
