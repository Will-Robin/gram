from itertools import compress

from gram.Classes import Network
from gram.Classes import Reaction
from gram.Classes import Compound
from gram.Classes import Substructure


def get_reactive_compounds(
    species_list: list[Compound], substructures: list[Substructure]
) -> list[Compound]:
    """
    Find species with the given subtructure.
    Returns a list of matching compounds

    Parameters
    ----------
    species_list: list[Compound]
        List of Compound objects from which
        reactive molecules are extracted.
    substructure: list[ Substructure]
        Reactive substructures.

    Returns
    -------
    matches: list[Compound]
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
    reactions: list[Reaction]
        List of Reactions
    invalid_substructures: list[Substructure]
        list of invalid Substructures.

    Returns
    -------
    reactions: list[Reaction]
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


def remove_reactions_by_product_substruct(network: Network, substruct: Substructure):
    """
    Remove the reactions in a nework if any of their products contain a defined
    substructure.

    Parameters
    ----------
    network: Network
        Network to check and modify.
    substruct: Substructure
        Substructure to check for in products.

    Returns
    -------
    None
    """

    remove_reactions = []

    for reaction in network.reactions:

        reaction_object = network.reactions[reaction]
        products = reaction_object.products

        for product in products:
            mol = network.compounds[product].mol

            if mol.HasSubstructMatch(substruct.mol):
                remove_reactions.append(reaction)

    network.remove_reactions(remove_reactions)
