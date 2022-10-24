from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol

from .smiles import incorrect_chiral_H_solve


def run_rdkit_reaction(
    reactant_compounds: list[Mol], reaction_template: AllChem.ChemicalReaction
) -> list[str]:
    """
    Performs a chemical reaction.

    Mapping can be added:
        Set atom map numbers of reactants before performing the reaction.
        Atom mappings can be removed:
        rdkit.Chem.rdChemReactions.RemoveMappingNumbersFromReactions(
                                                 (ChemicalReaction)reaction
                                                 ) → None

    Parameters
    ----------
    reactant_compounds: list[Mol]
        tuple of molecules which take part in the reaction.
    reaction_template: list[reaction]

    Returns
    -------
    reactions: list[str]
        A list of reaction SMILES.
    """

    reactions = []

    # Run the reaction to give a list of products sets
    product_sets = reaction_template.RunReactants(reactant_compounds)

    for product_set in product_sets:
        for product in product_set:
            # The products are checked for chiral information which has been
            # transferred to achiral carbons, which is removed.
            product = incorrect_chiral_H_solve(product)
            Chem.SanitizeMol(product)

        rxn = AllChem.ChemicalReaction()  # Create an empty chemical reaction
        [rxn.AddReactantTemplate(r) for r in reactant_compounds]
        [rxn.AddProductTemplate(p) for p in product_set]

        reaction_SMILES = AllChem.ReactionToSmiles(rxn)

        reactions.append(reaction_SMILES)

    return reactions
