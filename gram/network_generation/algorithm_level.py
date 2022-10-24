import itertools

from gram.Classes import Network
from gram.Classes import Compound
from gram.Classes import Substructure
from gram.Classes import ReactionTemplate
from gram.chemoinformatics.run_reaction import run_reaction

from .inspect import get_reactive_compounds
from .inspect import check_reaction_input


def extend_network_specific(
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

    template_substructures = reaction_template.reactant_substructures
    compounds_in_network = list(network.compounds.values())

    reactive_substrs = [Substructure(x) for x in template_substructures]

    reactants = get_reactive_compounds(compounds_in_network, reactive_substrs)

    for reactant in reactants:
        insert = [reactant] + reagents

        input_valid = check_reaction_input(insert, reactive_substrs)

        if input_valid:
            resulting_reactions = run_reaction(insert, reaction_template)
            network.add_reactions(resulting_reactions)
        else:
            pass


def extend_network_self(network: Network, reaction_template: ReactionTemplate) -> None:
    """
    Extend the network using any members of the network which can interact
    according to the reaction template provided.

    Parameters
    ----------
    network: Classes.Network
        Network to be extrapolated from. Modified in place.
    reagents: list[Classes.Compound]
        Reagents to be applied to the network.
    reaction_template: Classes.Reaction_Template
        Reaction template to be used on the network.
    secondary_substructure: Classes.Substructure
        Substructure of the second reaction component.

    Returns
    -------
    None
    """
    compounds_in_network = list(network.compounds.values())

    substructures = [Substructure(x) for x in reaction_template.reactant_substructures]

    reactants1 = get_reactive_compounds(compounds_in_network, [substructures[0]])

    reactants2 = get_reactive_compounds(compounds_in_network, [substructures[1]])

    for reactant_1 in reactants1:
        for reactant_2 in reactants2:

            insert = [reactant_1] + [reactant_2]
            resulting_reactions = run_reaction(insert, reaction_template)
            network.add_reactions(resulting_reactions)


def extend_network_task(network: Network, reaction_template: ReactionTemplate) -> None:
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
        substructure = Substructure(substruct)
        matches = get_reactive_compounds(current_compounds, [substructure])
        reactants.append(matches)

    # Build reactant combinations
    inputs = [list(prod) for prod in itertools.product(*reactants)]
    for inp in inputs:
        resulting_reactions = run_reaction(inp, reaction_template)
        network.add_reactions(resulting_reactions)
