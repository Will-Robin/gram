from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import ChemicalReaction


def reaction_from_smarts(smarts: str) -> ChemicalReaction:
    """
    Create a reaction from a SMARTS string.

    Parameters
    ----------
    smarts: str

    Returns
    -------
    reaction: ChemicalReaction
    """

    reaction: ChemicalReaction = AllChem.ReactionFromSmarts(smarts)

    return reaction


def reaction_to_smarts(reaction: ChemicalReaction) -> str:
    """
    Convert a ChemicalReaction to a reaction SMARTS string.

    Parameters
    ----------
    reaction: ChemicalReaction

    Returns
    -------
    smarts: str
    """

    smarts: str = AllChem.ReactionToSmarts(reaction)

    return smarts
