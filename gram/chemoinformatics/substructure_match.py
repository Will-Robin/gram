from gram.Classes import Compound
from gram.Classes import Substructure
from .core import substructure_match as substr


def has_substructure_match(compound: Compound, substructure: Substructure) -> bool:
    """
    Check if a Compound has a substructure match.
    """

    result = substr.has_substructure_match(compound.mol, substructure.mol)

    return result


def get_substructure_matches(compound: Compound, substructure: Substructure) -> tuple[tuple[int]]:
    """
    Get substructure matches for a Compound.
    """

    matches = substr.get_substructure_matches(compound.mol, substructure.mol)



    return matches
