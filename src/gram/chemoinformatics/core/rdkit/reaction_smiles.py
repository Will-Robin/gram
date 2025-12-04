from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import ChemicalReaction

from .smiles import canonicalise


def reaction_from_smiles(reaction_smiles: str) -> ChemicalReaction:
    """
    For converting reaction SMILES strings into RDKit ChemicalReaction
    objects.

    Parameters
    ----------
    smiles: str
        Reaction SMILES.

    Returns
    -------

    rdkit_reaction: ChemicalReaction
        Converted reaction.
    """

    split_rxn = reaction_smiles.split(">>")

    reactants = split_rxn[0].split(".")
    products = split_rxn[1].split(".")

    reactants.sort()
    products.sort()

    reactants_as_mol = [Chem.MolFromSmiles(r) for r in reactants]
    products_as_mol = [Chem.MolFromSmiles(p) for p in products]

    rdkit_reaction = ChemicalReaction()

    [rdkit_reaction.AddReactantTemplate(r) for r in reactants_as_mol]
    [rdkit_reaction.AddProductTemplate(p) for p in products_as_mol]

    return rdkit_reaction


def reaction_to_smiles(reaction: ChemicalReaction) -> str:
    """
    Convert a ChemicalReaction to a reaction SMILES string.

    Parameters
    ----------
    reaction: ChemicalReaction

    Returns
    -------
    reaction_smiles: str
    """

    reaction_smiles: str = AllChem.ReactionToSmiles(reaction)

    return reaction_smiles
