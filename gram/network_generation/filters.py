from gram.Classes import Network
from gram.Classes import Substructure


def remove_reactions_by_product_substruct(network: Network, substruct: Substructure):
    """
    Remove the reactions in a nework if any of their products contain a defined
    substructure.

    Parameters
    ----------
    network: Classes.Network
        Network to check and modify.
    substruct: Classes.Substructure
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
