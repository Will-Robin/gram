from gram.Classes import Substructure


def substructure_from_string(smarts: str) -> Substructure:
    """
    Create an Substructure object from a string.

    Parameters
    ----------
    smarts: str

    Returns
    -------
    Substructure
    """

    return Substructure(smarts)
