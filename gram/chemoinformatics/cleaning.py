from .core import cleanup as clean


def cleanup(reactions: list[str]) -> list[str]:
    """
    Clean up a list of reaction SMILES.

    Parameters
    ----------
    reactions: list[str]

    Returns
    -------
    list[str]

    """
    return clean.cleanup(reactions)
