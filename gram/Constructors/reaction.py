from gram.Classes import Reaction


def reaction_from_string(reaction_smiles: str) -> Reaction:
    """
    Create an Reaction object from a string.

    Parameters
    ----------
    reaction_smiles: str

    Returns
    -------
    Reaction
    """

    return Reaction(reaction_smiles)
