"""
Filters are used to determine if a reaction process can occur.
"""
from itertools import compress

from gram.Classes import Network
from gram.Classes import Reaction
from gram.Classes import Compound
from gram.Classes import Substructure

from gram.chemoinformatics import substructure_match as substr


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
            if substr.has_substructure_match(mol, substructure):
                matches.append(mol)
    return matches


def find_invalid_reactions(
    reactions: list[Reaction], substructure: Substructure, network: Network
) -> list[Reaction]:
    """
    Finds reactions with products that contain invalid substructures.

    Parameters
    ----------
    reactions: list[Reaction]
        List of Reactions
    substructure: Substructure
        list of invalid Substructures.
    network: Network
        Network containing reaction data.

    Returns
    -------
    reactions: list[Reaction]
        list of reactions that produce invalid substructures.
    """

    invalid_reactions = []

    for reaction in reactions:

        products = [network.compounds[x] for x in reaction.products]

        for product in products:
            if substr.has_substructure_match(product, substructure):
                invalid_reactions.append(reaction)

    return invalid_reactions
