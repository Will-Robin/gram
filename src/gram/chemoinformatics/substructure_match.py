from gram.Classes import Compound
from gram.Classes import Substructure
from .core import substructure_match


def has_substructure_match(compound: Compound, substructure: Substructure) -> bool:
    """
    Check if a Compound has a substructure match.

    Parameters
    ----------
    compound: Compound
    substructure: Substructure

    Returns
    -------
    bool
    """

    result = substructure_match.has_substructure_match(compound.mol, substructure.mol)

    return result


def get_substructure_matches(
    compound: Compound, substructure: Substructure
) -> tuple[tuple[int]]:
    """
    Get substructure matches for a Compound.

    Parameters
    ----------
    compound: Compound
    substructure: Substructure

    Returns
    -------
    tuple[tuple[int]]
    """

    matches = substructure_match.get_substructure_matches(
        compound.mol, substructure.mol
    )

    return matches
