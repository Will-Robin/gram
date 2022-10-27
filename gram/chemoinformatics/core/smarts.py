from rdkit import Chem
from rdkit.Chem.rdchem import Mol

def mol_from_smarts(smarts: str) -> Mol:
    """ """

    substructure: Mol = Chem.MolFromSmarts(smarts)

    return substructure


def mol_to_smarts(substructure: Mol) -> str:
    """ """

    smarts: str = Chem.MolToSmarts(substructure)

    return smarts
