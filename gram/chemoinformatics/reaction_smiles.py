from rdkit import Chem
from rdkit.Chem import AllChem
from .core.reaction_smiles import reaction_smiles_split


def reaction_from_smiles(reaction_smiles: str) -> AllChem.ChemicalReaction:
    """
    For converting reaction SMILES strings into RDKit AllChem.ChemicalReaction
    objects.

    Parameters
    ----------
    smiles: str
        Reaction SMILES.

    Returns
    -------

    rdkit_reaction: AllChem.ChemicalReaction
        Converted reaction.
    """

    reactants, products = reaction_smiles_split(reaction_smiles)

    reactants_as_mol = [Chem.MolFromSmiles(r) for r in reactants]
    products_as_mol = [Chem.MolFromSmiles(p) for p in products]

    rdkit_reaction = AllChem.ChemicalReaction()

    [rdkit_reaction.AddReactantTemplate(r) for r in reactants_as_mol]
    [rdkit_reaction.AddProductTemplate(p) for p in products_as_mol]

    return rdkit_reaction


def reaction_to_smiles(reaction: AllChem.ChemicalReaction) -> str:
    """ """

    reaction_smiles: str = AllChem.ReactionToSmiles(reaction)

    return reaction_smiles
