from gram.Classes import Network
from gram.Classes import Reaction


def load_network_from_reaction_list(
    reaction_list: list[str], name: str = "", description: str = ""
) -> Network:
    """
    Create a Network from a list of reactions

    Parameters
    ----------
    reaction_list: list[str]
        list of reaction SMILES strings
        Format: e.g. C=O.OC=C(O)CO>>O=C([C@@H](CO)O)CO

    Returns
    -------
    network: Network
    """

    rxns = []
    for reaction in reaction_list:
        rxns.append(Reaction(reaction))

    network = Network(rxns, name, description)

    return network
