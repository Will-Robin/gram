from gram.chemoinformatics import smarts


class Substructure:
    """
    Class to store substructures.
    """

    def __init__(self, SMARTS: str):
        """
        Parameters
        ----------
        SMARTS: str
            SMARTS corresponding to substructure.

        Attributes
        ----------
        mol: smarts.Mol
        smarts: str
        matching_compounds: list[str]
            Contextual list of compound tokens which contain this substructure.
        reaction_participations: list[str]
            List of reaction tokens for which the substructure is a reactant.
        producing_reactions: list[str]
            List of reaction tokens of which the substructure is a product.
        """

        self.mol: smarts.Mol = smarts.mol_from_smarts(SMARTS)

        self.smarts: str = SMARTS
        if smarts.mol_to_smarts(self.mol) is not None:
            self.smarts = smarts.mol_to_smarts(self.mol)

        self.matching_compounds: list[str] = []
        self.reaction_participations: list[str] = []
        self.producing_reactions: list[str] = []
