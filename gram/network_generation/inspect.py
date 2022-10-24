from itertools import compress

from gram.Classes import Network
from gram.Classes import Compound
from gram.Classes import Reaction
from gram.Classes import Substructure
from gram.Classes import ReactionTemplate


def get_reactive_compounds(
    species_list: list[Compound], substructures: list[Substructure]
) -> list[Compound]:
    """
    Find species with the given subtructure.
    Returns a list of matching compounds

    Parameters
    ----------
    species_list: list[Classes.Compound]
        List of Compound objects from which
        reactive molecules are extracted.
    substructure: list[Classes.Substructure]
        Reactive substructures.

    Returns
    -------
    matches: list[Classes.Compound]
        List of Compound objects which contain the substructure.
    """
    matches = []
    # testing each molecule in the list of species for reaction group matches
    for mol in species_list:
        for substructure in substructures:
            if mol.mol.HasSubstructMatch(substructure.mol):
                matches.append(mol)
    return matches


def remove_invalid_reactions(
    reactions: list[Reaction], invalid_substructures: list[Substructure]
):
    """
    Removes reactions with products that contain invalid substructures.

    Parameters
    ----------
    reactions: list[Classes.Reaction]
        List of Reactions
    invalid_substructures: list[Classes.Substructure]
        list of invalid Substructures.

    Returns
    -------
    reactions: list[Classes.Reaction]
        list of reactions with those that produce
        invalid substructures removed.
    """

    sortlist = []
    for reaction in reactions:
        tag = True
        for exc in invalid_substructures:

            products = reaction.reaction.GetProducts()

            substruct_matches = [p.HasSubstructMatch(exc.mol) for p in products]

            if any(substruct_matches):
                tag = False

        sortlist.append(tag)

    # new list with compounds tagged as false removed
    reactions = list(compress(reactions, sortlist))

    return reactions


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
    reactant_list: Classes.Compound

    reaction_template: Classes.ReactionTemplate

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
    compound: Compound, network: Network, reaction_template: ReactionTemplate
) -> bool:
    """
    Check if a reaction template has already been applied to a compound in a
    reaction network.

    Parameters
    ----------
    compound: NorthNet.Classes.Compound
    network: NorthNet.Classes.Network
    reaction_template: NorthNet.Classes.ReactionTemplate

    Returns
    -------
    reaction_performed: bool
        Whether reaction type has been applied to the compound or not.
    """

    reaction_name = reaction_template.name

    used_reactions = compound.reactant_in
    used_classes = [network.get_reaction_name(r) for r in used_reactions]
    reaction_performed = reaction_name in used_classes

    return reaction_performed
