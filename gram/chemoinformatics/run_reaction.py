from gram.Classes import Compound
from gram.Classes import Reaction
from gram.Classes import ReactionRule

from .core import run_reaction as core_run


def run_reaction(
    reactant_compounds: list[Compound],
    reaction_rule: ReactionRule,
    sanitize: bool = True,
) -> list[Reaction]:
    """
    Run a chemical reaction process.

    Parameters
    ----------
    reactant_compounds: list[Compound]
    reaction_rule: ReactionRule

    Returns
    -------
    output: list[Reaction]
    """

    reactants: list[core_run.Mol] = [r.mol for r in reactant_compounds]

    reactions = core_run.run_rdkit_reaction(
        reactants, reaction_rule.reaction, sanitize=sanitize
    )

    output: list[Reaction] = []
    for r in reactions:
        output.append(Reaction(r, reaction_rule=reaction_rule))

    return output
