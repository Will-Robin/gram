from pathlib import Path
from gram.Classes import ReactionRule


def load_reaction_rules(filename: str) -> dict[str, ReactionRule]:
    """ """

    text = Path(filename).read_text()

    lines = [line for line in text.split("\n") if line != ""]

    reaction_rules = {}
    for l in lines[1:]:
        entries = [x for x in l.split("\t") if x != ""]

        name = entries[0]
        reactant_substrs = entries[1].split(".")
        product_substrs = entries[2].split(".")
        reaction_smarts = entries[3]

        reaction_rules[name] = ReactionRule(
            name, reaction_smarts, reactant_substrs, product_substrs
        )

    return reaction_rules
