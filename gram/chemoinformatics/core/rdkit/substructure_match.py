from rdkit.Chem.rdchem import Mol


def has_substructure_match(molecule: Mol, substructure: Mol) -> bool:
    """
    Check if a molecule has a substructure match.
    """

    result: bool = molecule.HasSubstructMatch(substructure)

    return result


def get_substructure_matches(molecule: Mol, substructure: Mol) -> tuple[tuple[int]]:
    """
    Get substructure matches for a molecule.
    """

    matches: tuple[tuple[int]] = molecule.GetSubstructMatches(substructure)

    return matches
