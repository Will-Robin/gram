from rdkit.Chem.rdchem import Mol


def has_substructure_match(molecule: Mol, substructure: Mol) -> bool:
    """
    Check if a molecule has a substructure match.
    """

    result: bool = molecule.HasSubstructMatch(substructure)

    return result


def get_substructure_matches(molecule: Mol, substructure: Mol) -> list[Mol]:
    """
    Get substructure matches for a molecule.
    """

    matches: list[Mol] = molecule.GetSubstructMatches(substructure)

    return matches
