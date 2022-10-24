from gram.chemoinformatics import reaction_smarts


class ReactionTemplate:
    """
    Class for reaction templates.
    """

    def __init__(
        self,
        name: str,
        reaction_SMARTS: str,
        reactant_substructs: list[str],
        product_substructs: list[str],
    ):
        """
        Parameters
        ----------
        name: str
            name for reaction.
        reaction_SMARTS: str
            reaction string.
        reactant_substructs: list[str]
            reacting substructures
        product_substructs: list[str]
            substructures to which reactant substructures are converted.

        Attributes
        ----------
        name: str
        reaction: rdkit.Chem.AllChem.ChemicalReaction
        reaction_smarts: str
        reactant_substructures: list[str]
        product_substructures: list[str]
        """

        self.name = name
        if reaction_SMARTS != "":
            self.reaction = reaction_smarts.reaction_from_smarts(reaction_SMARTS)
            self.reaction_smarts = reaction_smarts.reaction_to_smarts(self.reaction)
        else:
            self.reaction = None
            self.reaction_smarts = ""

        self.reactant_substructures = reactant_substructs
        self.product_substructures = product_substructs
