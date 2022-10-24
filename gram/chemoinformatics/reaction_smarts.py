from rdkit.Chem import AllChem


def reaction_from_smarts(smarts: str) -> AllChem.ChemicalReaction:
    """ """

    reaction: AllChem.ChemicalReaction = AllChem.ReactionFromSmarts(smarts)

    return reaction


def reaction_to_smarts(reaction: AllChem.ChemicalReaction) -> str:
    """ """

    smarts: str = AllChem.ReactionToSmarts(reaction)

    return smarts
