"""
reaction_ratio.py
-----------------
Enforce a flux ratio (or ratio range) between two reactions using a
metabolite-based approach that is fully compatible with any cobrapy solver.

#THIS ONLY WORKS WITH IRREVERSIBLE REACTIOSN 

Architecture
------------
For each call, four model objects are created/managed:

rxn1   -> rep1
rxn2   -> rep2 
Ratio reaction 1
upper * rep_1 + rep2 ->
Ratio  reaction 2
lower * rep_1 + rep2 ->

lower  <=  v_1 / v_2  <=  upper

If Only one value is provided only one reaction will be created and ratio will be a single value 

When lower == upper (exact ratio), both reactions carry the same
stoichiometry and the constraint collapses to v_1 / v_2 = ratio exactly.


Example and test 

    test_model = cobra.Model("test_ratio")
    A = cobra.Metabolite("A", compartment="c",name="A")
    B = cobra.Metabolite("B", compartment="c",name="B")
    C = cobra.Metabolite("C", compartment="c",name="C")
 
    EX_A = cobra.Reaction("EX_A", lower_bound=1, upper_bound=10) #Force it to 1
    EX_A.add_metabolites({A: 1.0})
 
    rxn_1 = cobra.Reaction("rxn_1", lower_bound=0, upper_bound=100)
    rxn_1.add_metabolites({A: -1.0, B: 1.0})
 
    rxn_2 = cobra.Reaction("rxn_2", lower_bound=0, upper_bound=100)
    rxn_2.add_metabolites({A: -1.0, C: 1.0})
 
    EX_B = cobra.Reaction("EX_B", lower_bound=0, upper_bound=100)
    EX_B.add_metabolites({B: -1.0})
    
 
    EX_C = cobra.Reaction("EX_C", lower_bound=0, upper_bound=100)
    EX_C.add_metabolites({C: -1.0})
 
    test_model.add_reactions([EX_A, rxn_1, rxn_2, EX_B, EX_C])


    set_reaction_ratio(
    model= test_model,
    rxn_id_1= "rxn_1",
    rxn_id_2= "rxn_2",
    ratio= [2,15],
    constraint_name= "TestRatio",
)

    for reaction in test_model.reactions:
        print(reaction.id,get_equation(test_model,reaction.id))

        from cobra.sampling import sample
    n_samples=1000
    tol=1e-6
    rxn_id_1="rxn_1"
    rxn_id_2="rxn_2"
    samples = sample(test_model, n_samples, processes=1)
    v1, v2  = samples[rxn_id_1], samples[rxn_id_2]
    valid   = v2.abs() > tol
    ratios  = v1[valid] / v2[valid]
    print(min(ratios))
    print(max(ratios))  


"""

from __future__ import annotations

import math
from typing import Union

import cobra


# --------------------------------------------------------------------------- #
#  Public API                                                                  #
# --------------------------------------------------------------------------- #

def set_reaction_ratio(
    model: cobra.Model,
    rxn_id_1: str,
    rxn_id_2: str,
    ratio: Union[float, tuple[float, float]],
    constraint_name: str,
    ratio_reaction_bound: float= 1000,
    verbose: bool = True

) -> None:
    """
    Enforce  lower <= flux(rxn_id_1) / flux(rxn_id_2) <= upper  by adding a
    pair of reporter metabolites and two ratio reactions to the model.

    If the ratio reactions identified by ``constraint_name`` already exist
    in the model, their stoichiometry is **updated in-place** to the new
    ratio.  No duplicate objects are created.

    Parameters
    ----------
    model : cobra.Model
        The metabolic model to modify (edited in-place).
    rxn_id_1 : str
        Reaction ID of the *numerator* reaction.
    rxn_id_2 : str
        Reaction ID of the *denominator* reaction.
    ratio : float  or  (float, float)
        * **float** – enforce an *exact* ratio: flux(rxn_1) / flux(rxn_2) = ratio.
        * **(lower, upper)** – enforce a ratio *range*:
          lower =< flux(rxn_1) / flux(rxn_2) =< upper.
    constraint_name : str
        Base name used as a unique identifier (reporter label).  All model
        objects created by this function are prefixed with this string so
        they can be found and updated later.

    Raises
    ------
    KeyError
        If either reaction ID is not found in the model.
    ValueError
        If any ratio value is non-finite, or if lower > upper.
    
    """
    # ------------------------------------------------------------------ #
    #  1. Parse and validate ratio                                         #
    # ------------------------------------------------------------------ #
    if isinstance(ratio, (int, float)):
        lower = upper = float(ratio)
        flag_ratio_single_value=True
    else:
        lower, upper = float(ratio[0]), float(ratio[1])
        flag_ratio_single_value=False

    for val in (lower, upper):
        if not math.isfinite(val):
            raise ValueError(f"Ratio values must be finite; got {val!r}.")
    if lower > upper:
        raise ValueError(f"lower ({lower}) must be =< upper ({upper}).")

    # ------------------------------------------------------------------ #
    #  2. Fetch source reactions (raises KeyError if absent)              #
    # ------------------------------------------------------------------ #
    rxn_1 = model.reactions.get_by_id(rxn_id_1)
    rxn_2 = model.reactions.get_by_id(rxn_id_2)

    #Raise exception if either can be reversible
    if (rxn_1.lower_bound<0 or rxn_2.lower_bound<0) and flag_ratio_single_value==False:
      #Technically would also work as along as both were only negative
      raise Exception("This function only works for irreversible reactions. Please make sure both reactions have a lower bound of 0 or higher.")

    # ------------------------------------------------------------------ #
    #  3. Derive object IDs from the constraint name                      #
    # ------------------------------------------------------------------ #
    #Names if ratio is is single value
    single_ratio_rxn_id=constraint_name    
    #Names if ratio is an interval 
    upper_rxn_id = f"{constraint_name}_upper"
    lower_rxn_id = f"{constraint_name}_lower"

    rep1_id      = f"{constraint_name}__{rxn_id_1}__reporter"
    rep2_id      = f"{constraint_name}__{rxn_id_2}__reporter"
    #Sanity check that we are not tying to add double ratio to a model that has single ratio
    if (flag_ratio_single_value and upper_rxn_id in model.reactions) or (flag_ratio_single_value==False and single_ratio_rxn_id in model.reactions):
       #Remove reactions and start from scratch
       remove_reaction_ratio(model, constraint_name)  
    if flag_ratio_single_value:
       upper_rxn_id= lower_rxn_id = single_ratio_rxn_id
       
    # ------------------------------------------------------------------ #
    #  4. Add or update                                                    #
    # ------------------------------------------------------------------ #
    if upper_rxn_id in model.reactions and lower_rxn_id in model.reactions: 
        _update_ratio(
            model,
            rep1_id, rep2_id,
            upper_rxn_id, lower_rxn_id,
            lower, upper,
            constraint_name,verbose=verbose
        )
    else:
        _add_ratio(
            model,
            rxn_1, rxn_2,
            rep1_id, rep2_id,
            upper_rxn_id, lower_rxn_id,
            lower, upper,
            constraint_name,
             ratio_reaction_bound=ratio_reaction_bound,verbose=verbose
        )


# --------------------------------------------------------------------------- #
#  Internal helpers                                                            #
# --------------------------------------------------------------------------- #

def _add_ratio(
    model: cobra.Model,
    rxn_1: cobra.Reaction,
    rxn_2: cobra.Reaction,
    rep1_id: str,
    rep2_id: str,
    upper_rxn_id: str,
    lower_rxn_id: str,
    lower: float,
    upper: float,
    constraint_name: str,
    ratio_reaction_bound: float,
    verbose: bool = True
) -> None:
    """Create reporter metabolites and ratio reactions from scratch."""

    # ---- Reporter metabolites ----------------------------------------- #
    # Reuse existing metabolites if present (partial setup edge case).
    if rep1_id in model.metabolites:
        rep1 = model.metabolites.get_by_id(rep1_id)
    else:
        rep1 = cobra.Metabolite(
            rep1_id,
            name=f"{constraint_name} reporter for {rxn_1.id}",
            compartment="ratio_reporter",
        )
        rxn_1.add_metabolites({rep1: 1.0})  # rxn_1 produces rep1

    if rep2_id in model.metabolites:
        rep2 = model.metabolites.get_by_id(rep2_id)
    else:
        rep2 = cobra.Metabolite(
            rep2_id,
            name=f"{constraint_name} reporter for {rxn_2.id}",
            compartment="ratio_reporter",
        )
        rxn_2.add_metabolites({rep2: 1.0})  # rxn_2 produces rep2

    # ---- Upper ratio reaction ----------------------------------------- #
    # Stoichiometry: -upper * rep1  -1 * rep2  →  (nothing)
    # When only this reaction is active: v_1 / v_2 = upper.
    upper_rxn = cobra.Reaction(
        upper_rxn_id,
        name=f"{constraint_name} upper ratio ({upper})",
        lower_bound=0.0,
        upper_bound=ratio_reaction_bound#cobra.Configuration().upper_bound,
    )
    upper_rxn.add_metabolites({rep1: -upper, rep2: -1.0})
    model.add_reactions([upper_rxn])
    # ---- Lower ratio reaction ----------------------------------------- #
    # Stoichiometry: -lower * rep1  -1 * rep2  →  (nothing)
    # When only this reaction is active: v_1 / v_2 = lower.
    if(upper_rxn_id!=lower_rxn_id): #Only add if they are two reactions
        lower_rxn = cobra.Reaction(
            lower_rxn_id,
            name=f"{constraint_name} lower ratio ({lower})",
            lower_bound=0.0,
            upper_bound=ratio_reaction_bound##cobra.Configuration().upper_bound,
        )
        lower_rxn.add_metabolites({rep1: -lower, rep2: -1.0})
        model.add_reactions([lower_rxn])
    else:
        #Rename the single reaction
        upper_rxn.name=f"{constraint_name} ratio ({upper})"
        #Ratio can be reversible
        upper_rxn.bounds=(-1*ratio_reaction_bound,ratio_reaction_bound)
        
    if verbose:
       print(
           f"Added '{constraint_name}': "
           f"{lower} <= flux({rxn_1.id}) / flux({rxn_2.id}) <= {upper}"
       )


def _update_ratio(
    model: cobra.Model,
    rep1_id: str,
    rep2_id: str,
    upper_rxn_id: str,
    lower_rxn_id: str,
    lower: float,
    upper: float,
    constraint_name: str,
    verbose: bool = True
) -> None:
    """Update stoichiometry of existing ratio reactions to new bounds."""

    rep1      = model.metabolites.get_by_id(rep1_id)
    rep2      = model.metabolites.get_by_id(rep2_id)
    upper_rxn = model.reactions.get_by_id(upper_rxn_id)
    lower_rxn = model.reactions.get_by_id(lower_rxn_id)

    # combine=False replaces existing coefficients for these metabolites.
    # rep2 coefficient stays -1; only the rep1 coefficient may change.
    upper_rxn.add_metabolites({rep1: -upper, rep2: -1.0}, combine=False)
    upper_rxn.name = f"{constraint_name} upper ratio ({upper})"
    if upper_rxn_id!=lower_rxn_id:
        lower_rxn.add_metabolites({rep1: -lower, rep2: -1.0}, combine=False)
        lower_rxn.name = f"{constraint_name} lower ratio ({lower})"

    if verbose:
       print(
           f"Updated '{constraint_name}': "
           f"{lower} <= flux(rxn_1) / flux(rxn_2) <= {upper}"
       )


# --------------------------------------------------------------------------- #
#  Utility: remove a ratio constraint entirely                                 #
# --------------------------------------------------------------------------- #

def remove_reaction_ratio(model: cobra.Model, constraint_name: str, verbose: bool = True) -> None:
    """
    Remove all model objects (reporter metabolites + ratio reactions) that
    were created by :func:`enforce_reaction_ratio` for *constraint_name*.

    The stoichiometric coefficients for the reporter metabolites are also
    removed from the original source reactions.

    Parameters
    ----------
    model : cobra.Model
    constraint_name : str
        The same base name passed to ``enforce_reaction_ratio``.
    """
    ratio_rxn_ids = [f"{constraint_name}_upper", f"{constraint_name}_lower",constraint_name]
    #reporter_met_ids = [f"{constraint_name}__{rxn_id_1}__reporter", f"{constraint_name}__reporter_lower",f"{constraint_name}__reporter"]

    # Remove ratio reactions first (they reference the reporter metabolites).
    rxns_to_remove = [
        model.reactions.get_by_id(r)
        for r in ratio_rxn_ids
        if r in model.reactions
    ]
    #We can get the metabolites to remove from the ratio reactions
    reporter_met_ids=set()
    for rxn in rxns_to_remove:
        for met in rxn.metabolites:
            reporter_met_ids.add(met.id)
    
    model.remove_reactions(rxns_to_remove, remove_orphans=False)

    # Remove reporter metabolites from source reactions, then from the model.
    for met_id in reporter_met_ids:
        if met_id not in model.metabolites:
            continue
        met = model.metabolites.get_by_id(met_id)
        for rxn in list(met.reactions):
            rxn.add_metabolites({met: 0.0}, combine=False)  # zero coefficient → dropped
        model.remove_metabolites([met])

    if verbose:
       print(f"[remove_reaction_ratio] Removed constraint '{constraint_name}'.")
