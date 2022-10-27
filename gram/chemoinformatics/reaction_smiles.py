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
    """ """

    return r_sm.reaction_to_smiles(reaction)
