from rdkit import Chem
from rdkit.Chem.rdchem import Mol


def mol_from_smarts(smarts: str) -> Mol:
    """
    Create a Mol from a SMARTS string.

    Parameters
    ----------
    smarts: str

    Returns
    -------
    substructure: Mol
    """

    substructure: Mol = Chem.MolFromSmarts(smarts)

    return substructure


def mol_to_smarts(substructure: Mol) -> str:
    """
    Convert a Mol to a SMARTS string.

    Parameters
    ----------
    substructure: Mol

    Returns
    -------
    smarts: str
    """

    smarts: str = Chem.MolToSmarts(substructure)

    return smarts
