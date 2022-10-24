from .smiles import canonicalise


def canonicalise_reaction_smiles(reaction_smiles: str) -> str:
    """
    Canonicalise the consituent SMILES of a reaction SMILES string.

    Parameters
    ----------
    reaction_smiles: str

    Returns
    -------
    canonicalised_rxn_smiles: str
    """

    reactants, products = reaction_smiles_split(reaction_smiles)

    canonical_reacts = [canonicalise(r) for r in reactants]
    canonical_prods = [canonicalise(p) for p in products]

    lhs = ".".join(canonical_reacts)
    rhs = ".".join(canonical_prods)

    canonicalised_rxn_smiles = f"{lhs}>>{rhs}"

    return canonicalised_rxn_smiles


def reaction_smiles_split(reaction_smiles: str) -> tuple[list[str], list[str]]:
    """
    Parameters
    ----------
    reaction_smiles: str


    Returns
    -------
    reactants, products: tuple(list[str])
    """

    split_rxn = reaction_smiles.split(">>")

    reactants = split_rxn[0].split(".")
    products = split_rxn[1].split(".")

    reactants.sort()
    products.sort()

    return reactants, products
