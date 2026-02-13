import pandas as pd

def prune_mapping_dict_based_on_net_flux_direction(
    expanded_reaction_mapping_dict,
    reference_net_fluxes
):
    """
    Filter reaction mappings to retain only the active direction based on net flux values.
    
    This function processes an expanded reaction mapping dictionary and keeps only the 
    forward or reverse reactions based on the sign of the corresponding net flux. Reactions 
    with zero flux are excluded from the output.
    
    Parameters
    ----------
    expanded_reaction_mapping_dict : dict
        Dictionary mapping net flux IDs to their forward and reverse reactions.
        Structure: {net_flux_id: {"forward_reactions": [...], "reverse_reactions": [...]}}
    reference_net_fluxes : dict
        Dictionary mapping net flux IDs to their flux values.
        Positive values indicate forward direction, negative indicate reverse.
    
    Returns
    -------
    dict
        Updated mapping dictionary containing only the active direction for each reaction.
        Structure: {net_flux_id: {"forward_reactions": [...]}} for positive flux
                  {net_flux_id: {"reverse_reactions": [...]}} for negative flux
    
    Raises
    ------
    KeyError
        If a net_flux_id from expanded_reaction_mapping_dict is not found in 
        reference_net_fluxes. This is intentional to prevent using mismatched flux data.
    
    Notes
    -----
    - Reactions with net_flux_value == 0 are intentionally excluded from the output
    - The KeyError for missing net_flux_ids is intentional to ensure data consistency
    
    Examples
    --------
    >>> mapping = {"R1": {"forward_reactions": ["R1_FWD"], "reverse_reactions": ["R1_REV"]}}
    >>> fluxes = {"R1": 5.0}
    >>> result = prune_mapping_dict_based_on_net_flux_direction(mapping, fluxes)
    >>> result
    {'R1': {'forward_reactions': ['R1_FWD']}}
    """
    updated_expanded_reaction_mapping_dict = {}
    
    for net_flux_id, reaction_data in expanded_reaction_mapping_dict.items():
        # Will raise KeyError if net_flux_id is missing - this is intentional
        # to prevent using wrong/mismatched flux data
        net_flux_value = reference_net_fluxes[net_flux_id]
        
        # Skip reactions with zero flux
        if net_flux_value == 0:
            continue
            
        # Keep only the active direction based on flux sign
        if net_flux_value > 0:
            updated_expanded_reaction_mapping_dict[net_flux_id] = {
                "forward_reactions": reaction_data["forward_reactions"]
            }
        else:  # net_flux_value < 0
            updated_expanded_reaction_mapping_dict[net_flux_id] = {
                "reverse_reactions": reaction_data["reverse_reactions"]
            }
    
    return updated_expanded_reaction_mapping_dict


def get_enzyme_usage_bounds_and_components_from_gene_expression(
    model,
    gene_expression_dict,
    enzyme_kcat_scaling_factor_dict,
    gene_expression_to_enzyme_factor=1,
    reactions_to_omit=[],
    proteins_to_omit=[],
    verbose=False
):
    """
    Calculate enzyme usage bounds based on gene expression data.
    
    Function is a slight variation of from set_enzyme_usage_bounds_from_gene_expression
    
    This function converts gene expression levels to enzyme usage bounds by applying
    scaling factors. It processes enzyme-catalyzed reactions in the model and computes
    both the scaled protein levels and their usage bounds.
    
    Parameters
    ----------
    model : cobra.Model
        Genome-scale metabolic model containing enzyme-catalyzed reactions
    gene_expression_dict : dict
        Gene expression data mapping gene IDs to expression values
    enzyme_kcat_scaling_factor_dict : dict
        kcat scaling factors for each enzyme used in the EC model
    gene_expression_to_enzyme_factor : float, optional
        Conversion factor from gene expression to enzyme level (default: 1)
    reactions_to_omit : list, optional
        List of reaction IDs to exclude from calculations (default: [])
    proteins_to_omit : list, optional
        List of protein IDs to exclude from calculations (default: [])
    verbose : bool, optional
        If True, print detailed processing information (default: False)
    
    Returns
    -------
    tuple
        - scaled_protein_level_dict : dict
            Scaled protein levels {enzyme_id: scaled_level}
        - scaled_protein_usage_bounds_dict : dict
            Usage bounds for protein reactions {reaction_id: bound}
        - scaled_protein_level_components_dict : dict
            Components used to calculate protein levels {enzyme_id: components}
        - scaled_protein_usage_bounds_dict_components : dict
            Components used to calculate usage bounds {reaction_id: components}
        - missing_genes : list
            Genes not found in gene_expression_dict or with negative expression
    
    Notes
    -----
    - Enzyme usage reactions are expected to have exactly one gene
    - Protein usage bounds are negative as reactions are structured as: prot_A <--
    - Missing or negative gene expression values result in genes being skipped
    - Components dict contains the factors used in calculations for later reuse
    
    Raises
    ------
    Exception
        If an enzyme usage reaction has more than one gene associated with it
    """
    scaled_protein_level_dict = {}
    scaled_protein_usage_bounds_dict = {}
    scaled_protein_level_components_dict = {}
    scaled_protein_usage_bounds_dict_components = {}
    missing_genes = []
    
    # Expand proteins_to_omit with proteins from reactions_to_omit
    proteins_to_omit = proteins_to_omit.copy()
    for rid in reactions_to_omit:
        if rid in model.reactions:
            reaction = model.reactions.get_by_id(rid)
            reaction_proteins = [
                x.id.replace("prot_", "") 
                for x in reaction.metabolites 
                if x.id.startswith("prot_")
            ]
            proteins_to_omit += reaction_proteins
        else:
            if verbose:
                print(f"Reaction {rid} to omit not in model")
    
    if verbose:
        print(f"Processing enzyme usage reactions (omitting {len(proteins_to_omit)} proteins)...")
    
    # Process each enzyme usage reaction
    processed_count = 0
    for reaction in model.reactions.query("usage_prot_"):
        enzyme = reaction.id.replace("usage_prot_", "")
        
        if enzyme in proteins_to_omit:
            if verbose:
                print(f"  Skipping enzyme {enzyme} (in omit list)")
            continue
        
        # Validate that reaction has exactly one gene
        genes = list(reaction.genes)
        if len(genes) != 1:
            raise Exception(f"Wrong number of genes ({len(genes)}) in {reaction.id}")
        
        gene = genes[0].id
        gene_expression = gene_expression_dict.get(gene, None)
        
        # Skip genes with no expression data or negative values
        if gene_expression is None or gene_expression < 0:
            missing_genes.append(gene)
            if verbose:
                print(f"  No expression for gene {gene} (enzyme {enzyme})")
            continue
        
        # Calculate scaled enzyme usage
        max_enzyme_usage = gene_expression_to_enzyme_factor * gene_expression
        enzyme_kcat_factor = enzyme_kcat_scaling_factor_dict.get(enzyme, 1)
        scaled_max_enzyme_usage = max_enzyme_usage * enzyme_kcat_factor
        
        scaled_protein_level_dict[enzyme] = scaled_max_enzyme_usage
        # Usage bound is negative due to reaction structure: prot_A <--
        scaled_protein_usage_bounds_dict[reaction.id] = scaled_max_enzyme_usage * -1
        
        # Store components for later reuse
        aggregate_factor = -1 * enzyme_kcat_factor * gene_expression_to_enzyme_factor
        scaled_protein_level_components_dict[enzyme] = {
            "gene": gene,
            "gene_expression_to_enzyme_factor": gene_expression_to_enzyme_factor,
            "enzyme_kcat_scaling_factor": enzyme_kcat_factor,
            "aggregate_gene_to_prot_usage_factor": aggregate_factor
        }
        scaled_protein_usage_bounds_dict_components[reaction.id] = (
            scaled_protein_level_components_dict[enzyme].copy()
        )
        processed_count += 1
    
    if verbose:
        print(f"  Processed {processed_count} enzymes successfully")
        print(f"  Missing genes: {len(missing_genes)}")
    
    return (
        scaled_protein_level_dict,
        scaled_protein_usage_bounds_dict,
        scaled_protein_level_components_dict,
        scaled_protein_usage_bounds_dict_components,
        missing_genes
    )


def get_net_reaction_capacities(
    model,
    gene_expression_data,
    enzyme_kcat_scaling_factor,
    expanded_reaction_mapping_dict=None,
    mapping_dict_includes_both_directions=None,
    gene_expression_to_enzyme_factor=1,
    gene_expression_fold_change=None,
    fold_change_is_Log2=False,
    reactions_to_omit=[],
    proteins_to_omit=[],
    net_reaction_capacity_components_dict=None,
    protein_usage_capacities_components_dict=None,
    save_components=False,
    verbose=False
):
    """
    Calculate maximum and minimum net capacity for each reaction based on gene expression.
    
    This function computes reaction capacities by converting gene expression levels to
    enzyme activities, accounting for enzyme kinetics (kcat) and stoichiometry. It can
    operate in two modes: (1) calculate from scratch, or (2) reuse pre-calculated 
    component dictionaries for faster repeated calculations.
    
    Parameters
    ----------
    model : cobra.Model
        Genome-scale metabolic model with enzyme-catalyzed reactions
    gene_expression_data : dict or pd.Series
        Gene expression values mapping gene IDs to expression levels
    enzyme_kcat_scaling_factor : dict or pd.Series
        kcat scaling factors for each enzyme (from EC model construction)
    expanded_reaction_mapping_dict : dict, optional
        Mapping of net reactions to their expanded component reactions.
        If None, will auto-generate with default patterns (default: None)
    mapping_dict_includes_both_directions : bool or None, optional
        Whether mapping dict contains both forward and reverse reactions.
        If None, auto-detected from mapping structure (default: None)
    gene_expression_to_enzyme_factor : float, optional
        Conversion factor from gene expression to enzyme concentration (default: 1)
    gene_expression_fold_change : dict or pd.Series, optional
        Fold changes to apply to gene expression values (default: None)
    fold_change_is_Log2 : bool, optional
        If True, fold changes are in log2 scale and will be converted (default: False)
    reactions_to_omit : list, optional
        Reaction IDs to exclude from calculations (default: [])
    proteins_to_omit : list, optional
        Protein IDs to exclude from calculations (default: [])
    net_reaction_capacity_components_dict : dict, optional
        Pre-calculated component dictionary for faster reuse mode (default: None)
    protein_usage_capacities_components_dict : dict, optional
        Pre-calculated protein usage components for reuse mode (default: None)
    save_components : bool, optional
        If True, return component dictionaries for later reuse (default: False)
    verbose : bool, optional
        If True, print detailed processing information (default: False)
    
    Returns
    -------
    tuple
        - net_reaction_capacity_dict : dict
            Reaction capacities. Structure depends on mapping_dict_includes_both_directions:
            - If True: {net_reaction_id: {"forward": value, "reverse": value}}
            - If False: {net_reaction_id: value}
        - net_reaction_capacity_components_dict : dict or None
            Components for capacity calculations (None if save_components=False)
        - protein_usage_capacities_dict : dict
            Protein usage capacities {reaction_id: capacity}
        - protein_usage_capacities_components_dict : dict or None
            Components for protein usage (None if save_components=False)
    
    Raises
    ------
    KeyError
        In reuse mode, if a gene that was present when components were created
        is now missing from gene_expression_data. This is intentional to prevent
        silent errors from data mismatches.
    
    Notes
    -----
    **Calculation Modes:**
    
    - Fresh calculation: When component dicts are None, calculates from scratch
    - Reuse mode: When component dicts are provided, uses pre-calculated factors
      for much faster computation. In this mode, missing genes will raise KeyError.
    
    **Fold Changes:**
    
    - Applied ONLY to genes present in gene_expression_fold_change dict
    - If log2 scale, converted to linear: fold_change = 2^(log2_fold_change)
    - Only modifies genes that exist in both fold_change and gene_expression_data
    - Does NOT use .get() - intentionally only processes explicitly listed genes
    
    **Enzyme Complexes:**
    
    - For multi-subunit complexes, the minimum capacity across subunits is used
    - This represents the rate-limiting subunit in the complex
    
    **Output Structure:**
    
    - Single value mode (mapping_dict_includes_both_directions=False):
      Returns {reaction_id: capacity_value}
      Better for creating pandas DataFrames with one column per sample
    - Dual value mode (mapping_dict_includes_both_directions=True):
      Returns {reaction_id: {"forward": value, "reverse": value}}
      Useful when you need both directions for analysis
    
    **Performance Optimization:**
    
    For analyzing multiple samples with the same model:
    1. First call: Set save_components=True
    2. Subsequent calls: Pass the returned component dicts for 10-100x speedup
    
    Examples
    --------
    >>> # Fresh calculation with component saving
    >>> capacities, comps, protein_usage, prot_comps = get_net_reaction_capacities(
    ...     model, gene_expr, kcat_factors, save_components=True, verbose=True
    ... )
    
    >>> # Reuse mode for second sample (much faster)
    >>> capacities2, _, protein_usage2, _ = get_net_reaction_capacities(
    ...     model, gene_expr2, kcat_factors,
    ...     net_reaction_capacity_components_dict=comps,
    ...     protein_usage_capacities_components_dict=prot_comps,
    ...     verbose=True
    ... )
    
    >>> # With fold changes (only genes in fold_change dict are modified)
    >>> capacities3, _, _, _ = get_net_reaction_capacities(
    ...     model, gene_expr, kcat_factors,
    ...     gene_expression_fold_change=fold_changes,
    ...     fold_change_is_Log2=True
    ... )
    """
    if verbose:
        print("=" * 70)
        print("NET REACTION CAPACITY CALCULATION")
        print("=" * 70)
    
    # Convert pandas Series to dictionaries if needed
    if isinstance(enzyme_kcat_scaling_factor, pd.Series):
        enzyme_kcat_scaling_factor = enzyme_kcat_scaling_factor.to_dict()
        if verbose:
            print(f"Converted enzyme_kcat_scaling_factor from Series to dict")
    
    if isinstance(gene_expression_data, pd.Series):
        gene_expression_data = gene_expression_data.to_dict()
        if verbose:
            print(f"Converted gene_expression_data from Series to dict")
    
    if isinstance(gene_expression_fold_change, pd.Series):
        gene_expression_fold_change = gene_expression_fold_change.to_dict()
        if verbose:
            print(f"Converted gene_expression_fold_change from Series to dict")
    
    if verbose:
        print(f"\nInput data summary:")
        print(f"  Genes in expression data: {len(gene_expression_data)}")
        print(f"  Enzyme kcat scaling factors: {len(enzyme_kcat_scaling_factor)}")
    
    # Apply fold changes to gene expression data
    if gene_expression_fold_change is not None:
        if verbose:
            print(f"\nApplying fold changes:")
            print(f"  Fold changes provided for: {len(gene_expression_fold_change)} genes")
        
        # Convert log2 fold changes to linear scale if needed
        if fold_change_is_Log2:
            if verbose:
                print(f"  Converting log2 fold changes to linear scale")
            gene_expression_fold_change = {
                gene: pow(2, fc) for gene, fc in gene_expression_fold_change.items()
            }
        
        # Apply fold changes only to genes in fold_change dict
        # INTENTIONAL: Do NOT use .get() - only modify genes explicitly in fold_change dict
        gene_expression_data = gene_expression_data.copy()
        genes_modified = 0
        genes_not_in_data = 0
        
        for gene in gene_expression_fold_change:
            if gene in gene_expression_data:
                original_value = gene_expression_data[gene]
                gene_expression_data[gene] *= gene_expression_fold_change[gene]
                genes_modified += 1
                if verbose and genes_modified <= 5:  # Show first 5 examples
                    print(f"    {gene}: {original_value:.2f} -> "
                          f"{gene_expression_data[gene]:.2f} "
                          f"(FC={gene_expression_fold_change[gene]:.2f})")
            else:
                genes_not_in_data += 1
        
        if verbose:
            print(f"  Applied fold changes to: {genes_modified} genes")
            if genes_not_in_data > 0:
                print(f"  Genes in fold_change but not in expression data: "
                      f"{genes_not_in_data}")
    
    # Auto-generate mapping dict if not provided
    if expanded_reaction_mapping_dict is None:
        if verbose:
            print("\nNo mapping dict provided - generating default mapping...")
        print("Building default expanded_reaction_mapping_dict with "
              "reverse_reaction_pattern=_REV, isoenzyme_reaction_pattern=_EXP_\\d+ "
              "and patterns_to_ommit=[^usage_prot_]")
        
        from .ec_base import get_ec_expanded_reaction_mapping
        expanded_reaction_mapping_dict, _ = get_ec_expanded_reaction_mapping(
            model,
            reverse_reaction_pattern="_REV",
            isoenzyme_reaction_pattern="_EXP_\\d+",
            patterns_to_ommit=["^usage_prot_"],
            verbose=False
        )
        mapping_dict_includes_both_directions = True
        if verbose:
            print(f"  Generated mapping for {len(expanded_reaction_mapping_dict)} reactions")
    
    # Auto-detect if mapping includes both directions
    if mapping_dict_includes_both_directions is None:
        n_directions = max([
            len(expanded_reaction_mapping_dict[x]) 
            for x in expanded_reaction_mapping_dict
        ])
        mapping_dict_includes_both_directions = (n_directions > 1)
        
        if verbose:
            print(f"\nAuto-detected mapping structure:")
            print(f"  Maximum directions per reaction: {n_directions}")
            print(f"  Includes both directions: {mapping_dict_includes_both_directions}")
    
    if verbose:
        if mapping_dict_includes_both_directions:
            print(f"  Output format: {{reaction_id: {{'forward': value, 'reverse': value}}}}")
        else:
            print(f"  Output format: {{reaction_id: value}}")
    
    net_reaction_capacity_dict = {}
    
    # === FRESH CALCULATION MODE ===
    if (net_reaction_capacity_components_dict is None or 
        protein_usage_capacities_components_dict is None):
        
        if verbose:
            print("\n" + "=" * 70)
            print("MODE: FRESH CALCULATION")
            print("=" * 70)
        
        # Get scaled gene expression for proteins and their components
        (scaled_protein_level_dict,
         protein_usage_capacities_dict,
         scaled_protein_level_components_dict,
         protein_usage_capacities_components_dict,
         missing_genes) = get_enzyme_usage_bounds_and_components_from_gene_expression(
            model=model,
            gene_expression_dict=gene_expression_data,
            enzyme_kcat_scaling_factor_dict=enzyme_kcat_scaling_factor,
            gene_expression_to_enzyme_factor=gene_expression_to_enzyme_factor,
            reactions_to_omit=reactions_to_omit,
            proteins_to_omit=proteins_to_omit,
            verbose=verbose
        )
        
        if verbose:
            print(f"\nEnzyme processing results:")
            print(f"  Proteins with expression data: {len(scaled_protein_level_dict)}")
            print(f"  Protein usage reactions: {len(protein_usage_capacities_dict)}")
            print(f"  Missing genes: {len(missing_genes)}")
            if len(missing_genes) > 0 and len(missing_genes) <= 10:
                print(f"    Missing: {', '.join(missing_genes)}")
            elif len(missing_genes) > 10:
                print(f"    First 10 missing: {', '.join(missing_genes[:10])}...")
        
        if verbose:
            print(f"\nCalculating net reaction capacities...")
            print(f"  Processing {len(expanded_reaction_mapping_dict)} net reactions")
        
        net_reaction_capacity_components_dict = {}
        reactions_processed = 0
        reactions_with_proteins = 0
        
        # Calculate capacities for each net reaction
        for net_reaction_id in expanded_reaction_mapping_dict:
            for reaction_type in expanded_reaction_mapping_dict[net_reaction_id]:
                net_reaction_involves_proteins = False
                net_reaction_capacity = 0
                net_reaction_capacity_components_local = {}
                reaction_ids = expanded_reaction_mapping_dict[net_reaction_id][reaction_type]
                
                # Process each component reaction
                for reaction_id in reaction_ids:
                    reaction_object = model.reactions.get_by_id(reaction_id)
                    reaction_protein_activities = []
                    reaction_protein_activities_components = []
                    
                    # Process each protein in the reaction
                    for metabolite in reaction_object.metabolites:
                        if metabolite.id.startswith("prot_"):
                            protein_id = metabolite.id.replace("prot_", "")
                            protein_coefficient = reaction_object.metabolites[metabolite]
                            scaled_expression_value = scaled_protein_level_dict.get(protein_id, None)
                            
                            # Skip if protein data is missing (gene not in expression data)
                            if scaled_expression_value is None:
                                continue
                            
                            # Calculate capacity contributed by this protein
                            # If coefficient is -0.1, scaled_expression_value*10 can be consumed
                            capacity_contributed_by_protein = (
                                -1 * scaled_expression_value / protein_coefficient
                            )
                            reaction_protein_activities.append(capacity_contributed_by_protein)
                            
                            # Save components for later reuse if requested
                            if save_components:
                                if protein_id in scaled_protein_level_components_dict:
                                    protein_component_src = (
                                        scaled_protein_level_components_dict[protein_id]
                                    )
                                    protein_component_dict = {
                                        "gene": protein_component_src["gene"],
                                        "gene_expression_to_enzyme_factor": 
                                            protein_component_src["gene_expression_to_enzyme_factor"],
                                        "enzyme_kcat_scaling_factor": 
                                            protein_component_src["enzyme_kcat_scaling_factor"],
                                        "protein_coefficient": protein_coefficient,
                                        # Final constant: gene expression -> capacity
                                        "aggregate_factor_gene_to_capacity": (
                                            protein_component_src["aggregate_gene_to_prot_usage_factor"] / 
                                            protein_coefficient
                                        )
                                    }
                                    reaction_protein_activities_components.append(
                                        protein_component_dict
                                    )
                    
                    # For enzyme complexes, use minimum capacity (rate-limiting subunit)
                    if len(reaction_protein_activities) > 0:
                        min_capacity = min(reaction_protein_activities)
                        net_reaction_capacity += min_capacity
                        net_reaction_capacity_components_local[reaction_id] = (
                            reaction_protein_activities_components
                        )
                        net_reaction_involves_proteins = True
                
                # Store capacity if reaction involves proteins
                if net_reaction_involves_proteins:
                    if net_reaction_id not in net_reaction_capacity_dict:
                        if mapping_dict_includes_both_directions:
                            net_reaction_capacity_dict[net_reaction_id] = {}
                        net_reaction_capacity_components_dict[net_reaction_id] = {}
                    
                    direction_key = reaction_type.replace("_reactions", "")
                    if mapping_dict_includes_both_directions:
                        net_reaction_capacity_dict[net_reaction_id][direction_key] = (
                            net_reaction_capacity
                        )
                    else:
                        net_reaction_capacity_dict[net_reaction_id] = net_reaction_capacity
                    
                    # Always store components with direction information
                    net_reaction_capacity_components_dict[net_reaction_id][direction_key] = (
                        net_reaction_capacity_components_local
                    )
                    reactions_with_proteins += 1
                
                reactions_processed += 1
        
        if verbose:
            print(f"  Processed: {reactions_processed} reaction directions")
            print(f"  Reactions with proteins: {reactions_with_proteins}")
            print(f"  Net reactions in output: {len(net_reaction_capacity_dict)}")
    
    # === REUSE MODE ===
    else:
        if verbose:
            print("\n" + "=" * 70)
            print("MODE: REUSE (using pre-calculated components)")
            print("=" * 70)
            print(f"  Net reactions in components dict: "
                  f"{len(net_reaction_capacity_components_dict)}")
            print(f"  Protein usage components: "
                  f"{len(protein_usage_capacities_components_dict)}")
        
        reactions_processed = 0
        genes_accessed = set()
        
        try:
            # Calculate capacities using pre-computed factors
            for net_reaction_id in net_reaction_capacity_components_dict:
                for reaction_type in net_reaction_capacity_components_dict[net_reaction_id]:
                    net_reaction_capacity = 0
                    
                    for reaction_id in net_reaction_capacity_components_dict[net_reaction_id][reaction_type]:
                        enzyme_subunits = (
                            net_reaction_capacity_components_dict[net_reaction_id][reaction_type][reaction_id]
                        )
                        
                        # Calculate capacity for each subunit
                        # INTENTIONAL: Will raise KeyError if gene is missing
                        # This prevents silent errors from data mismatches
                        subunit_capacities = []
                        for subunit in enzyme_subunits:
                            gene = subunit["gene"]
                            # Direct indexing - will raise KeyError if gene is missing
                            gene_expr = gene_expression_data[gene]
                            capacity = subunit["aggregate_factor_gene_to_capacity"] * gene_expr
                            subunit_capacities.append(capacity)
                            genes_accessed.add(gene)
                        
                        # Use minimum capacity for enzyme complexes (rate-limiting subunit)
                        reaction_protein_activity = min(subunit_capacities)
                        net_reaction_capacity += reaction_protein_activity
                    
                    # Store results
                    if net_reaction_id not in net_reaction_capacity_dict and mapping_dict_includes_both_directions:
                        net_reaction_capacity_dict[net_reaction_id] = {}
                    
                    if mapping_dict_includes_both_directions:
                        net_reaction_capacity_dict[net_reaction_id][reaction_type] = net_reaction_capacity
                    else:
                        net_reaction_capacity_dict[net_reaction_id] = net_reaction_capacity
                    
                    reactions_processed += 1
        
        except KeyError as e:
            raise KeyError(
                f"Gene '{e.args[0]}' was present when component dictionary was created "
                f"but is now missing from gene_expression_data. This indicates a data "
                f"mismatch. Please ensure you're using the correct gene expression "
                f"dataset or recreate the component dictionaries with save_components=True."
            ) from e
        
        # Calculate protein usage capacities using pre-computed factors
        try:
            protein_usage_capacities_dict = {}
            for rid in protein_usage_capacities_components_dict:
                gene = protein_usage_capacities_components_dict[rid]["gene"]
                # Direct indexing - will raise KeyError if gene is missing
                gene_expr = gene_expression_data[gene]
                aggregate_factor = protein_usage_capacities_components_dict[rid]["aggregate_gene_to_prot_usage_factor"]
                protein_usage_capacities_dict[rid] = aggregate_factor * gene_expr
                genes_accessed.add(gene)
        
        except KeyError as e:
            raise KeyError(
                f"Gene '{e.args[0]}' was present when component dictionary was created "
                f"but is now missing from gene_expression_data. This indicates a data "
                f"mismatch. Please ensure you're using the correct gene expression "
                f"dataset or recreate the component dictionaries with save_components=True."
            ) from e
        
        if verbose:
            print(f"\nReuse mode results:")
            print(f"  Reaction capacities calculated: {reactions_processed}")
            print(f"  Unique genes accessed: {len(genes_accessed)}")
            print(f"  Protein usage capacities: {len(protein_usage_capacities_dict)}")
            print(f"  Net reactions in output: {len(net_reaction_capacity_dict)}")
        
        # Don't return components unless we generated them
        save_components = False
    
    # Clean up component dicts if not saving
    if save_components == False:
        net_reaction_capacity_components_dict = None
        protein_usage_capacities_components_dict = None
        if verbose:
            print("\nComponent dictionaries NOT returned (save_components=False)")
    else:
        if verbose:
            print("\nComponent dictionaries returned for reuse (save_components=True)")
            print(f"  To reuse: pass these dicts to net_reaction_capacity_components_dict "
                  f"and protein_usage_capacities_components_dict parameters")
    
    if verbose:
        print("=" * 70)
        print("CALCULATION COMPLETE")
        print("=" * 70)
    
    return (
        net_reaction_capacity_dict,
        net_reaction_capacity_components_dict,
        protein_usage_capacities_dict,
        protein_usage_capacities_components_dict
    )