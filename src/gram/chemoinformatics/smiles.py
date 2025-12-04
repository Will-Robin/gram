from .core import smiles as core_sm


def mol_from_smiles(smiles: str) -> core_sm.Mol:
    """
    Create a molecule object from a SMILES string.

    Parameters
    ----------
    smiles: str

    Returns
    -------
    core_sm.Mol
    """

    return core_sm.mol_from_smiles(smiles)


def canonicalise(smiles: str) -> str:
    """
    Canonicalise a SMILES string.

    Parameters
    ----------
    smiles: str

    Returns
    -------
    str
    """

    return core_sm.canonicalise(smiles)
