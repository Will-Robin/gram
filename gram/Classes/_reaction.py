from gram.chemoinformatics.reaction_smiles import reaction_from_smiles
from gram.chemoinformatics.core import reaction_smiles
from gram.chemoinformatics.core.reaction_smiles import reaction_smiles_split
from gram.chemoinformatics.core.reaction_smiles import canonicalise_reaction_smiles

from ._reactionRule import ReactionRule


class Reaction:
    """
    A class representing chemical reactions
    """

    def __init__(
        self,
        rxn_smiles: str,
        reaction_rule=None,
        info: dict[str, str] = dict(),
    ):
        """
        Parameters
        ----------
        rxn_smiles: str
            Reaction object for the reaction. Please provide a valid reaction
            SMILES string with canonicalised SMILES reactants and products.
        reaction_rule: ReactionRule object
            Reaction rule which created the reaction.
        info: dict
            Dictionary of information (e.g. database entries)

        Attributes
        ----------
        reaction_smiles: str
            Reaction SMILES string
        reaction_rule: ReactionRule or None
            Reaction rule for the reaction. Defaults to None
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

        if reaction_rule is None:
            self.reaction_rule = ReactionRule("none", "", [], [])
        else:
            self.reaction_rule = reaction_rule

        self.data = info

        self.reactants, self.products = reaction_smiles_split(canonical_rxn_smiles)
