import itertools

from gram.Classes import Network
from gram.Classes import Compound
from gram.Classes import ReactionTemplate
from gram.chemoinformatics.run_reaction import run_reaction

from .filters import get_reactive_compounds
from .inspect import check_reaction_input


def apply_specific_reaction_to_network(
    network: Network, reagents: list[Compound], reaction_template: ReactionTemplate
) -> None:
    """
    Extend the network using a single reagent set and reaction template.

    Parameters
    ----------
    network: Classes.Network
        Network to be extrapolated from. Modified in place.
    reagents: list[Classes.Compound]
        Reagents to be applied to the network.
    reaction_template: Classes.ReactionTemplate
        Reaction template to be used on the network.

    Returns
    -------
    None
    """

    rxn_substructs = reaction_template.reactant_substructures
    compounds_in_network = list(network.compounds.values())

    reactants = get_reactive_compounds(compounds_in_network, rxn_substructs)

    for reactant in reactants:
        insert = [reactant] + reagents

        input_valid = check_reaction_input(insert, rxn_substructs)

        if input_valid:
            resulting_reactions = run_reaction(insert, reaction_template)
            network.add_reactions(resulting_reactions)
        else:
            pass


def apply_reaction_to_network(
    network: Network, reaction_template: ReactionTemplate
) -> None:
    """
    Extend the network using any members of the network which can interact
    according to the reaction template provided.

    Parameters
    ----------
    network: Classes.Network
        Network to be extrapolated from. Modified in place.
    reagents: list[Compound]
        Reagents to be applied to the network.
    reaction_template: Classes.Reaction_Template
        Reaction template to be used on the network.
    secondary_substructure: Classes.Substructure
        Substructure of the second reaction component.

    Returns
    -------
    None
    """

    current_compounds = list(network.compounds.values())
    reactants = []
    for substruct in reaction_template.reactant_substructures:
        matches = get_reactive_compounds(current_compounds, [substruct])
        reactants.append(matches)

    # Build reactant combinations
    inputs = [list(prod) for prod in itertools.product(*reactants)]
    for inp in inputs:
        resulting_reactions = run_reaction(inp, reaction_template)
        network.add_reactions(resulting_reactions)
