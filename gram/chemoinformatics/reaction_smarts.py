from .core import reaction_smarts as r_sm


def reaction_from_smarts(smarts: str) -> r_sm.ChemicalReaction:
    """
    Create a ChemicalReaction object from a SMARTS string.

    Parameters
    ----------
    smarts: str

    Returns
    -------
    r_sm.ChemicalReaction
    """

    return r_sm.reaction_from_smarts(smarts)


def reaction_to_smarts(reaction: r_sm.ChemicalReaction) -> str:
    """
    Convert a ChemicalReaction to a reaction SMARTS string.

    Parameters
    ----------
    reaction: r_sm.ChemicalReaction

    Returns
    -------
    str
    """

    return r_sm.reaction_to_smarts(reaction)
