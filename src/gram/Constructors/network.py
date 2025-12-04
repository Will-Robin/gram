from gram.Classes import Network
from gram.Classes import Reaction


def network_from_string(net_string: str) -> Network:
    """
    Create an Network object from a string.

    Parameters
    ----------
    net_string: str

    Returns
    -------
    Network
    """
    return Network([], net_string, "")


def network_from_reaction_smiles(
    reaction_smiles: list[str], name: str = "", description: str = ""
) -> Network:
    """
    Create a network from a list of reaction SMILES strings.

    Parameters
    ----------
    reaction_smiles: list[str]
    name: str
    description: str

    Returns
    -------
    Network
    """

    reaction_list = []
    for r in reaction_smiles:
        reaction_list.append(Reaction(r))

    return Network(reaction_list, name, description)
