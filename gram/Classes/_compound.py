from gram.chemoinformatics.smiles import mol_from_smiles
from gram.chemoinformatics.smiles import Mol


class Compound:
    """
    A class to store compound information.
    """

    def __init__(self, smiles: str):
        """
        Parameters
        ----------
        smiles: str
            SMILES corresponding to compound.

        Attributes
        ----------
        mol: smiles.Mol
            RDKit representation of the compound to allow for chemoinformatics
            operations.
        product_of: list[str]
            List of reaction tokens for which the compound is a product.
        reactant_in: list[str]
            List of reaction tokens for which the compound is a reactant.
        reactive_substructures: list[str]
            List of keys to contextual reactive substructures which either
            engage in, or are products of, reactions.
        """

        self.smiles: str = smiles
        self.mol: Mol = mol_from_smiles(smiles)
        self.product_of: list[str] = []
        self.reactant_in: list[str] = []
        self.reactive_substructures: list[str] = []
