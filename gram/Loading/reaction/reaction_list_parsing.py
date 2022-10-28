"""
For loading reaction information fron text files.
"""
from gram.Classes import ReactionRule


def load_reaction_rules_from_file(
    fname: str, delimiter: str = "\t"
) -> dict[str, ReactionRule]:
    """
    Reads reaction rules from a .csv file.

    Assumed that the file is structure as follows:
    header\n
    name\treactant SMARTS\tproduct SMARTS\tReaction SMARTS
    etc.

    (ignores anything beyond 4th column)

    Parameters
    ----------
    fname: str
        file containing reaction rules and substructures.
    delimiter: str or pathlib Path
        Column delimiter for the file

    Returns
    -------
    reaction_rules: dict
        Dictionary of reaction rules.
        {reaction class name: ReactionRule}
    """

    lines = []
    with open(fname, "r") as file:
        for c, line in enumerate(file):
            lines = file.readlines()

    reaction_rules = {}
    for c, line in enumerate(lines):
        if c == 0:
            pass
        else:
            ins = line.strip("\n").split(delimiter)
            reaction_rules[ins[0]] = ReactionRule(
                ins[0], ins[3], ins[1].split("."), ins[2].split(".")
            )
    return reaction_rules
