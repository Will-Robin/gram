from rdkit import Chem
from rdkit.Chem.rdchem import Mol as Mol


def mol_from_smiles(smiles: str) -> Mol:
    """
    Create
    """

    mol = Chem.MolFromSmiles(smiles)

    return mol
