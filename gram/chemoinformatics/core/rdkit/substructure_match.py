from rdkit.Chem.rdchem import Mol


def has_substructure_match(molecule: Mol, substructure: Mol) -> bool:
    """
    Check if a molecule has a substructure match.

    Parameters
    ----------
    molecule: Mol
    substructure: Mol

    Returns
    -------
    result: bool
    """

    result: bool = molecule.HasSubstructMatch(substructure)

    return result


def get_substructure_matches(molecule: Mol, substructure: Mol) -> tuple[tuple[int]]:
    """
    Get substructure matches for a molecule.

    Parameters
    ----------
    molecule: Mol
    substructure: Mol

    Returns
    -------
    matches: tuple[tuple[int]]
    """

    matches: tuple[tuple[int]] = molecule.GetSubstructMatches(substructure)

    return matches
