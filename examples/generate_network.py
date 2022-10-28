import yaml

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


# Load in some variables
with open("params.yaml", "r") as file:
    text = file.read()

info = yaml.load(text, Loader=yaml.FullLoader)

"""Get reaction components"""
reaction_SMARTS_file = info["reaction-smarts-file"]

reactions = load_reaction_rules(reaction_SMARTS_file)

C_patt = Substructure("[C]")
count_carbons = lambda x: substr.get_substructure_matches(x, C_patt)

"""Name"""
network_name = info["network-name"]
description = info["network-description"]

"""Boundary conditions"""
# iterations  overshoot for C6, but do so to get
# all reaction paths and compounds possible up
# to C6 compounds
iterations = info["iterations"]
start_smiles = info["initiator-smiles"]

initiator_species = [Compound(x) for x in start_smiles]
reaction_network = Network([], network_name, description)

reaction_network.add_compounds(initiator_species)

"""Reactivity constructor"""
reaction_pattern = info["reaction-rules"]
deprotonation_rules = [x for x in reaction_pattern if "deprotonation" in x]
protonation_rules = [x for x in reaction_pattern if "protonation" in x]

"""Expansion operation"""
x = 0
while x < iterations:
    for task in reaction_pattern:
        n_gen.apply_reaction_to_network(reaction_network, reactions[task])
    generate_epimers(
        reaction_network,
        deprotonation_rules=[reactions[d] for d in deprotonation_rules],
        protonation_rules=[reactions[p] for p in protonation_rules],
    )

    """Removing compounds > C6"""
    # In effect, this is equivalent to setting all chain-growing reaction
    # rules to not occur for C6 compounds
    # i.e. [$(C(O)=CO)!$(C(O)=C(O)C(O)C(O)C(O)CO)], etc.
    remove_compounds = [
        reaction_network.compounds[c]
        for c in reaction_network.compounds
        if len(count_carbons(reaction_network.compounds[c])) > 6
    ]
    reaction_network.remove_compounds(remove_compounds)

    x += 1

compound_number = len(reaction_network.compounds)
reaction_number = len(reaction_network.reactions)

print(f"Generated {compound_number} compounds and {reaction_number} reactions.")

# save the reactions and their reaction rules
header = [
    "ReactionSMILES",
    "reaction name",
    "reactionSMARTS",
    "reacting substructures",
    "product substructures",
]

reaction_text = "\t".join(header) + "\n"
for reaction in reaction_network.reactions:
    rxn = reaction_network.reactions[reaction]
    rule = rxn.reaction_rule
    reaction_text += reaction + "\t"
    reaction_text += rule.name + "\t"

    reaction_text += rule.reaction_smarts + "\t"
    reactant_substrs = [r.smarts for r in rule.reactant_substructures]
    product_substrs = [p.smarts for p in rule.product_substructures]
    reaction_text += ".".join(reactant_substrs) + "\t"
    reaction_text += ".".join(product_substrs) + "\t"

    reaction_text += "\n"


filename = "reactions_and_rules.tsv"
with open(filename, "w") as file:
    file.write(reaction_text)
