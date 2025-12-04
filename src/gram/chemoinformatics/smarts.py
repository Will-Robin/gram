from .core import smarts as core_smarts


def mol_from_smarts(smarts: str) -> core_smarts.Mol:
    """
    Create a molecule object from a SMARTS string.

    Parameters
    ----------
    smarts: str

    Returns
    -------
    Mol
    """

    return core_smarts.mol_from_smarts(smarts)


def mol_to_smarts(substructure: core_smarts.Mol) -> str:
    """
    Convert a Mol object to a SMARTS string.

    Parameters
    ----------
    smarts: core_smarts.Mol

    Returns
    -------
    str
    """

    return core_smarts.mol_to_smarts(substructure)
