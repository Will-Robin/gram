from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import ChemicalReaction


def reaction_from_smarts(smarts: str) -> ChemicalReaction:
    """ """

    reaction: ChemicalReaction = AllChem.ReactionFromSmarts(smarts)

    return reaction


def reaction_to_smarts(reaction: ChemicalReaction) -> str:
    """ """

    smarts: str = AllChem.ReactionToSmarts(reaction)

    return smarts
