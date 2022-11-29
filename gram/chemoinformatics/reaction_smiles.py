from .core import smiles as sm
from .core import reaction_smiles as r_sm


def reaction_from_smiles(reaction_smiles: str) -> r_sm.ChemicalReaction:
    """
    For converting reaction SMILES strings into RDKit AllChem.ChemicalReaction
    objects.

    Parameters
    ----------
    smiles: str
        Reaction SMILES.

    Returns
    -------
    rdkit_reaction: r_sm.ChemicalReaction
        Converted reaction.
    """

    return r_sm.reaction_from_smiles(reaction_smiles)


def reaction_to_smiles(reaction: r_sm.ChemicalReaction) -> str:
    """
    Convert a ChemicalReaction to a reaction SMILES string.

    Parameters
    ----------
    reaction: r_sm.ChemicalReaction

    Returns
    -------
    str
    """

    return r_sm.reaction_to_smiles(reaction)


def split_reaction_smiles(reaction_smiles: str) -> tuple[list[str], list[str]]:
    """
    Split a reaction SMILES string into its components.

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

    reactants, products = split_reaction_smiles(reaction_smiles)

    canonical_reacts = [sm.canonicalise(r) for r in reactants]
    canonical_prods = [sm.canonicalise(p) for p in products]

    lhs = ".".join(canonical_reacts)
    rhs = ".".join(canonical_prods)

    canonicalised_rxn_smiles = f"{lhs}>>{rhs}"

    return canonicalised_rxn_smiles
