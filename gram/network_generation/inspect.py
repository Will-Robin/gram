"""
These functions check if criteria are met for a reaction.
"""
from gram.Classes import Network
from gram.Classes import Compound
from gram.Classes import Substructure
from gram.Classes import ReactionRule


def check_reaction_input(
    reactant_list: list[Compound], reactive_substructs: list[Substructure]
) -> bool:
    """
    Checks for valid reaction input by checking the reactants
    contain the supplied substructures in the same order.

    Returns False if reaction inputs are incompatible with
    the reactive substructures provided.

    Parameters
    ----------
    reactant_list: Compound

    reaction_rule: ReactionRule

    Returns
    -------
    bool
        Whether the supplied list of compounds is compatible with the
        reaction substructure order supplied.
    """

    test_list = [
        r.mol.HasSubstructMatch(s) for r, s in zip(reactant_list, reactive_substructs)
    ]

    return all(test_list)


def check_reaction_occurence(
    compound: Compound, network: Network, reaction_rule: ReactionRule
) -> bool:
    """
    Check if a reaction rule has already been applied to a compound in a
    reaction network.

    Parameters
    ----------
    compound: Compound
    network: Network
    reaction_rule: ReactionRule

    Returns
    -------
    reaction_performed: bool
        Whether reaction type has been applied to the compound or not.
    """

    # Reaction to check for
    reaction_name = reaction_rule.name

    # Get the tokens of reactions in which the compound is a reactant
    used_reactions = compound.reactant_in
    # Get the names of the reaction rules for each reaction
    used_classes = [network.get_reaction_name(r) for r in used_reactions]

    return reaction_name in used_classes
