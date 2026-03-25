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
    
    # Collect variables and coefficients for objective
    coef_dict = {}
    
    # Forward/Standard reactions: forward_var gets +1, reverse_var gets -1
    for forward_rid in forward_reactions:
        r = _model.reactions.get_by_id(forward_rid)
        if r.forward_variable.ub > 0:
            coef_dict[r.forward_variable] = 1
        if r.reverse_variable.ub > 0:
            if len(reverse_reactions) > 0:
                raise Exception(
                    f"{forward_rid}: If using forward/reverse reactions, "
                    "forward reaction should not have negative lower bound"
                )
            coef_dict[r.reverse_variable] = -1
    
    # Reverse reactions: forward_var gets -1 (reverse the direction)
    for reverse_rid in reverse_reactions:
        r_rev = _model.reactions.get_by_id(reverse_rid)
        if r_rev.forward_variable.ub > 0:
            coef_dict[r_rev.forward_variable] = -1
        if r_rev.reverse_variable.ub > 0:
            raise Exception(
                f"{reverse_rid}: Reverse reactions should not have reverse fluxes allowed"
            )
    
    if not coef_dict:
        logger.warning(f"Skipping flux {flux_id} as no flux is allowed in either direction")
        return flux_id, float("nan")
    
    # Set objective coefficients directly
    _model.solver.objective.set_linear_coefficients(coef_dict)
    
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
        {var: 0 for var in coef_dict.keys()}
    )
    
    return flux_id, value


def flux_variability_analysis_net_flux(
    model: "Model",
    expanded_reaction_mapping_dict: Optional[Dict[str, Dict[str, List[str]]]] = None,
    flux_list: Optional[List[str]] = None,
    fraction_of_optimum: float = 1.0,
    processes: Optional[int] = None,
    solver_tolerance_feasibility: Optional[float] = None,
    solver_tolerance_optimality: Optional[float] = None,
    verbose: bool = True,
    progress_interval: Optional[int] = None,  # Print Update every N fluxes
) -> pd.DataFrame:
    """Determine the minimum and maximum net flux value for each flux.

    Note precision can be increased with: model.solver.configuration.tolerances.feasibility = 1e-9
    Setting processes to 1 might also increase reproducibility
    This is similar to standard FVA but operates on net fluxes defined by
    forward and reverse reaction mappings, useful for models with split
    reversible reactions or when you want to analyze net flux across
    multiple reactions.

    Parameters
    ----------
    model : cobra.Model
        The model for which to run the analysis. It will *not* be modified.
    expanded_reaction_mapping_dict : dict, optional
        A dictionary mapping flux IDs to forward and reverse reactions in the model.
    flux_list : list of str, optional
        The flux IDs for which to obtain min/max net fluxes. If None will use
        all fluxes in expanded_reaction_mapping_dict (default None).
    fraction_of_optimum : float, optional
        Must be <= 1.0. Requires that the objective value is at least the
        fraction times maximum objective value (default 1.0).
    processes : int, optional
        The number of parallel processes to run (default None).
    solver_tolerance_feasibility : float, optional
        If provided, temporarily sets the solver's feasibility tolerance.
    solver_tolerance_optimality : float, optional
        If provided, temporarily sets the solver's optimality tolerance.
    verbose : bool, optional
        Print progress information (default True).
    progress_interval : int, optional
        Print progress every N fluxes. If None, prints every 10% (default None).

    Returns
    -------
    pandas.DataFrame
        A data frame with flux identifiers as the index and two columns:
        - maximum: indicating the highest possible net flux
        - minimum: indicating the lowest possible net flux
    """
    import time
    
    # Check that the model doesn't already contain reporter metabolites
    #existing_reporter_mets = [met.id for met in model.metabolites if met.id.startswith("reporter_met_")]
    #
    #if existing_reporter_mets:
    #    raise ValueError(
    #        f"Model already contains reporter metabolites for net fluxes: {existing_reporter_mets}. "
    #        "Should not run flux_variability_analysis_net_flux on a model with reporters. "
    #        "Use the original model without reporters."
    #    )
    
    # Auto-generate mapping if not provided
    if expanded_reaction_mapping_dict is None:
        from .ec_base import get_ec_expanded_reaction_mapping
        if verbose:
            print("Building default expanded_reaction_mapping_dict with reverse_reaction_pattern=_REV, isoenzyme_reaction_pattern=_EXP_\\d+ and patterns_to_ommit=[^usage_prot_]")
        expanded_reaction_mapping_dict, _ = get_ec_expanded_reaction_mapping(
            model, 
            reverse_reaction_pattern="_REV", 
            isoenzyme_reaction_pattern="_EXP_\\d+",
            patterns_to_ommit=["^usage_prot_"],
            verbose=False
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
    
    if verbose:
        print(f"Calculating FVA for {len(flux_list)} net fluxes...")
    
    if processes is None:
        processes = configuration.processes

    num_fluxes = len(flux_list)
    processes = min(processes, num_fluxes)
    
    # Determine progress update interval
    if progress_interval is None:
        # Default: update every 10%
        progress_interval = max(1, num_fluxes // 10)
    
    if verbose:
        print(f"Progress updates every {progress_interval} fluxes")

    fva_result = pd.DataFrame(
        {
            "minimum": np.zeros(num_fluxes, dtype=float),
            "maximum": np.zeros(num_fluxes, dtype=float),
        },
        index=flux_list,
    )
    
    # Adjust solver tolerances if specified
    old_tolerance = None
    old_optimality_tolerance = None
    
    if solver_tolerance_feasibility is not None:
        if verbose:
            print(f"Setting solver feasibility tolerance to {solver_tolerance_feasibility} for FVA")
        old_tolerance = model.solver.configuration.tolerances.feasibility
        model.solver.configuration.tolerances.feasibility = solver_tolerance_feasibility
    
    if solver_tolerance_optimality is not None:
        if verbose:
            print(f"Setting solver optimality tolerance to {solver_tolerance_optimality} for FVA")
        old_optimality_tolerance = model.solver.configuration.tolerances.optimality
        model.solver.configuration.tolerances.optimality = solver_tolerance_optimality
    
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
            if verbose:
                print(f"\nCalculating {what} fluxes...")
            
            start_time = time.time()
            completed = 0
            last_update = 0
            
            if processes > 1:
                # Better chunk size for load balancing
                chunk_size = min(len(flux_list) // processes,1000)

                
                if verbose:
                    print(f"  Using {processes} processes with chunk_size={chunk_size}")
                
                with ProcessPool(
                    processes,
                    initializer=_init_worker_net_flux,
                    initargs=(model, what[:3], expanded_reaction_mapping_dict),
                ) as pool:
                    for flux_id, value in pool.imap_unordered(
                        _fva_step_net_flux, flux_list, chunksize=chunk_size
                    ):
                        fva_result.at[flux_id, what] = value
                        completed += 1
                        
                        # Print progress every progress_interval fluxes
                        if verbose and (completed - last_update >= progress_interval or completed == num_fluxes):
                            elapsed = time.time() - start_time
                            rate = completed / elapsed if elapsed > 0 else 0
                            percent = 100 * completed / num_fluxes
                            eta = (num_fluxes - completed) / rate if rate > 0 else 0
                            
                            print(f"  {completed}/{num_fluxes} ({percent:.1f}%) | "
                                  f"{rate:.1f} flux/s | "
                                  f"Elapsed: {elapsed:.0f}s | "
                                  f"ETA: {eta:.0f}s")
                            last_update = completed
            else:
                # Single process mode
                _init_worker_net_flux(
                    model, what[:3], expanded_reaction_mapping_dict
                )
                
                for flux_id, value in map(_fva_step_net_flux, flux_list):
                    fva_result.at[flux_id, what] = value
                    completed += 1
                    
                    # Print progress every progress_interval fluxes
                    if verbose and (completed - last_update >= progress_interval or completed == num_fluxes):
                        elapsed = time.time() - start_time
                        rate = completed / elapsed if elapsed > 0 else 0
                        percent = 100 * completed / num_fluxes
                        eta = (num_fluxes - completed) / rate if rate > 0 else 0
                        
                        print(f"  {completed}/{num_fluxes} ({percent:.1f}%) | "
                              f"{rate:.1f} flux/s | "
                              f"Elapsed: {elapsed:.0f}s | "
                              f"ETA: {eta:.0f}s")
                        last_update = completed
            
            if verbose:
                elapsed = time.time() - start_time
                print(f"  Completed in {elapsed:.1f}s ({num_fluxes/elapsed:.1f} flux/s)")
    
    # Restore solver tolerances
    if old_tolerance is not None:
        model.solver.configuration.tolerances.feasibility = old_tolerance
    if old_optimality_tolerance is not None:
        model.solver.configuration.tolerances.optimality = old_optimality_tolerance
    
    return fva_result[["minimum", "maximum"]]
