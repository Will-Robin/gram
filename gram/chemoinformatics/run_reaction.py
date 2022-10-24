from rdkit.Chem.rdchem import Mol

from gram.Classes import Compound
from gram.Classes import Reaction
from gram.Classes import ReactionTemplate

from .core.run_reaction import run_rdkit_reaction


def run_reaction(
    reactant_compounds: list[Compound], reaction_template: ReactionTemplate
) -> list[Reaction]:

    reactants: list[Mol] = [r.mol for r in reactant_compounds]

    reactions = run_rdkit_reaction(reactants, reaction_template.reaction)

    output: list[Reaction] = []
    for r in reactions:
        output.append(Reaction(r, reaction_template=reaction_template))

    return output
