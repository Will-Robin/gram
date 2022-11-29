from rdkit import Chem


def cleanup(reactions: list[str]) -> list[str]:
    """
    Used to clear up the output of rdkit reaction function, parsing multiple
    reaction outcomes. Not perfect.

    Parameters
    ----------
    reactions: list[str]
        list of reaction SMILES strings

    Returns
    -------
    reactions_out: list[str]
        list of cleaned reaction SMILES strings
    """

    reactions = [r.replace("[H+][O-]", "O") for r in reactions]
    reactions = [r.replace("[O-][H+]", "O") for r in reactions]
    reactions = [r.replace("/", "") for r in reactions]
    reactions = [r.replace(r"\\", "") for r in reactions]

    reactions_out = []
    for _, reaction in enumerate(reactions):
        # Split the reaction up for inspection
        LHS = reaction.split(">>")[0].split(".")
        RHS = reaction.split(">>")[1].split(".")
        ins = ""
        for l_elem in LHS:
            mol = Chem.MolFromSmiles(l_elem)
            smiles = Chem.MolToSmiles(mol, canonical=True)
            ins += smiles + "."

        ins = ins.strip(".")
        ins += ">>"

        for r_elem in RHS:
            mol = Chem.MolFromSmiles(r_elem)
            smiles = Chem.MolToSmiles(mol, canonical=True)
            ins += smiles + "."

        reactions_out.append(ins.strip("."))

    return reactions_out
