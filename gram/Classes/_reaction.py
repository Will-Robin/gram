from gram.chemoinformatics.reaction_smiles import reaction_from_smiles
from gram.chemoinformatics.core import reaction_smiles
from gram.chemoinformatics.core.reaction_smiles import reaction_smiles_split
from gram.chemoinformatics.core.reaction_smiles import canonicalise_reaction_smiles

from ._reactionTemplate import ReactionTemplate


class Reaction:
    """
    A class representing chemical reactions
    """

    def __init__(
        self,
        rxn_smiles: str,
        reaction_template=None,
        info: dict[str, str] = dict(),
    ):
        """
        Parameters
        ----------
        rxn_smiles: str
            Reaction object for the reaction. Please provide a valid reaction
            SMILES string with canonicalised SMILES reactants and products.
        reaction_template: ReactionTemplate object
            Reaction template which created the reaction.
        info: dict
            Dictionary of information (e.g. database entries)

        Attributes
        ----------
        reaction_smiles: str
            Reaction SMILES string
        reaction_template: ReactionTemplate or None
            Reaction template for the reaction. Defaults to None
        data: dict
            Dictionary of information (e.g. database entries)
        reactants: list[str]
            List of reactant SMILES tokens.
        products: list[str]
            List of product SMILES tokens.
        """

        canonical_rxn_smiles: str = canonicalise_reaction_smiles(rxn_smiles)

        self.reaction_smiles = canonical_rxn_smiles
        self.reaction = reaction_from_smiles(rxn_smiles)

        if reaction_template is None:
            self.reaction_template = ReactionTemplate("none", "", [], [])
        else:
            self.reaction_template = reaction_template

        self.data = info

        self.reactants, self.products = reaction_smiles_split(canonical_rxn_smiles)
