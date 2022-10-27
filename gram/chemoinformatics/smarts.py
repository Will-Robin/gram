from .core import smarts as core_smarts


def mol_from_smarts(smarts: str) -> core_smarts.Mol:
    """ """

    return core_smarts.mol_from_smarts(smarts)


def mol_to_smarts(substructure: core_smarts.Mol) -> str:
    """ """

    return core_smarts.mol_to_smarts(substructure)
