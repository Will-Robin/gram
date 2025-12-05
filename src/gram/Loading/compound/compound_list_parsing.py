"""
For loading compound objects from a text file.
"""

from gram.Classes import Compound


def load_compounds_from_file(
    fname: str, name_col: int = 0, SMILES_col: int = 1, delimiter: str = ","
) -> dict[str, Compound]:
    """
    Reads compounds from a .csv file.

    Parameters
    ----------
    fname: str
        Path to the file containing nformation.
    name_col: int
        Column in the file which will give the keys for the output dict
    SMILES_col: int
        Column containing the compound SMILES

    Returns
    -------
    reagents: dict
        Dictionary containing extracted compounds.
        {SMILES string: Compound object}
    """

    reagents = {}
    with open(fname, "r") as file:
        for c, line in enumerate(file):
            if c == 0:
                pass
            else:
                ins = line.strip("\n").split(delimiter)
                name = ins[name_col]

                reagents[name] = Compound(ins[SMILES_col])

    return reagents
