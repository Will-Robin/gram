from gram.Classes import Compound
from gram.Classes import Reaction
from gram.Classes import ReactionTemplate

from .core import run_reaction as core_run


def run_reaction(
    reactant_compounds: list[Compound], reaction_template: ReactionTemplate
) -> list[Reaction]:

    reactants: list[core_run.Mol] = [r.mol for r in reactant_compounds]

    reactions = core_run.run_rdkit_reaction(reactants, reaction_template.reaction)

    output: list[Reaction] = []
    for r in reactions:
        output.append(Reaction(r, reaction_template=reaction_template))

    return output
