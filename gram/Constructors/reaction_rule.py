from gram.Classes import ReactionRule


def reaction_rule_from_string(reaction_smarts: str) -> ReactionRule:
    """
    Create an ReactionRule object from a string.

    Parameters
    ----------
    reaction_smarts: str

    Returns
    -------
    ReactionRule
    """

    return ReactionRule("", reaction_smarts, [], [])
