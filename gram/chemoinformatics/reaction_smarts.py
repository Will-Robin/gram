from .core import reaction_smarts as r_sm


def reaction_from_smarts(smarts: str) -> r_sm.ChemicalReaction:
    """ """

    return r_sm.reaction_from_smarts(smarts)


def reaction_to_smarts(reaction: r_sm.ChemicalReaction) -> str:
    """ """

    return r_sm.reaction_to_smarts(reaction)
