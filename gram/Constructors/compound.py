from gram.Classes import Compound


def compound_from_string(smiles: str) -> Compound:
    """
    Create a Compound object from a string.

    Parameters
    ----------
    smiles: str

    Returns
    -------
    Compound
    """

    return Compound(smiles)
