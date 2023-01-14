"""
An example which generates the formose reaction network similar to that used in
Robinson, Daines, van Duppen, de Jong, Huck, *Nature Chemistry*, 2022, **14**,
623-631. DOI: 10.1038/s41557-022-00956-7

This script demonstrates a lot of the functionality in `gram`, from loading
reaction rules from text, to creating an applying transformations to a reaction
network.
"""
import yaml
from pathlib import Path

from gram.Classes import Network
from gram.Classes import Compound
from gram.Classes import Substructure
from gram.Classes import ReactionRule

from gram import network_generation as n_gen
from gram.chemoinformatics import substructure_match as substr

from loading import load_reaction_rules


def generate_epimers(
    network: Network,
    deprotonation_rules: list[ReactionRule] = [],
    protonation_rules: list[ReactionRule] = [],
) -> None:
    """
    A function to generate all of the epimers of a sugar: over-iterates
    multiple times through a series of protonation/deprotonation reactions.

    Parameters
    ----------
    network: Network

    deprotonation_rules: list[ReactionRule]

    protonation_rules: list[ReactionRule]

    Returns
    -------
    None
    """

    hydroxide = Compound("[OH-]")
    water = Compound("O")

    i = 0
    reaction_number = len(network.reactions)
    while i < 0:
        for d_rule in deprotonation_rules:
            n_gen.apply_specific_reaction_to_network(network, [hydroxide], d_rule)

        for p_rule in protonation_rules:
            n_gen.apply_specific_reaction_to_network(network, [water], p_rule)

        new_reaction_number = len(network.reactions)

        i = new_reaction_number - reaction_number

        reaction_number = new_reaction_number


# Load in variables
info = yaml.load(Path("params.yaml").read_text(), Loader=yaml.FullLoader)

# Create ReactionRule objects from reaction SMARTS strings.
reactions = load_reaction_rules(info["reaction-smarts-file"])

# A function for counting carbons.
C_patt = Substructure("[C]")
count_carbons = lambda x: substr.get_substructure_matches(x, C_patt)

# Information for reaction network generation
network_name = info["network-name"]
description = info["network-description"]
iterations = info["iterations"]
initiator_smiles = info["initiator-smiles"]

# Initialise reaction network
reaction_network = Network([], network_name, description)
initiator_species = [Compound(x) for x in initiator_smiles]
reaction_network.add_compounds(initiator_species)

# Reaction templates to use (by name)
reaction_pattern = info["reaction-rules"]
deprotonation_rules = [x for x in reaction_pattern if "deprotonation" in x]
protonation_rules = [x for x in reaction_pattern if "protonation" in x]

# Expansion operation
for _ in range(iterations):

    # Apply reaction rules to compatible compounds in the network.
    for reaction_type in reaction_pattern:
        n_gen.apply_reaction_to_network(reaction_network, reactions[reaction_type])

    # Migrate carbonyls
    generate_epimers(
        reaction_network,
        deprotonation_rules=[reactions[d] for d in deprotonation_rules],
        protonation_rules=[reactions[p] for p in protonation_rules],
    )

    # Removal of compounds > C6
    remove_compounds = [
        reaction_network.compounds[c]
        for c in reaction_network.compounds
        if len(count_carbons(reaction_network.compounds[c])) > 6
    ]
    reaction_network.remove_compounds(remove_compounds)

compound_number = len(reaction_network.compounds)
reaction_number = len(reaction_network.reactions)

print(f"Generated {compound_number} compounds and {reaction_number} reactions.")
