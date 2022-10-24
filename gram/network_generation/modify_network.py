from rdkit import Chem
from gram.Classes import InputProcess
from gram.Classes import OutputProcess
from gram.Classes import Network
from gram.Classes import Substructure
from gram.Classes import Reaction


def add_flow_inputs(network: Network, inputs: list[str]) -> None:
    """
    Add flow inputs into network.

    Parameters
    ----------
    network: NorthNet.Classes.Network
        Network to add flow terms to.
    inputs: list[str]
        List of input string tokens to add.

    Returns
    -------
    None
    """

    for inp in inputs:
        i_obj = InputProcess(f"{inp}_#0>>{inp}")
        network.add_input_process(i_obj)

    for compound in network.compounds:
        o_obj = OutputProcess(f"{compound}>>Sample")
        network.add_output_process(o_obj)


def skip_compound(network: Network, substructure: Substructure) -> None:
    """
    TODo: refactor

    Finds compounds in the network which have the specified substructure and
    removes them from the network (including the reactions for which they are
    reactants and products). New reactions between the reactant compounds which
    created the removed compounds and those which are products of reactions of
    the removed compounds are created to fill the gaps.

    Parameters
    ----------
    network: Network
        Network to modify

    Returns
    -------
    None
    """

    reaction_removal_list = []
    compounds_removal_list = []
    new_reaction_list = []

    for compound in network.compounds:

        if network.compounds[compound].mol.HasSubstructMatch(substructure.mol):
            compounds_removal_list.append(compound)
            # Find reactions connected to the compound to remove
            # and their reactants and products

            # Cycle through the incoming reactions, add them into a list for
            # removal.
            for i in network.compounds[compound].product_of:
                reaction_removal_list.append(i)

                # Find the reactants for this reaction and any
                # other products produced with the removed compound
                in_rs = network.reactions[i].reactants
                outs = []

                for product in network.reactions[i].products:
                    if product != compound:
                        outs.append(product)

                # Go through the outgoing reactions add them into a list for
                # removal
                for reaction in network.compounds[compound].reactant_in:
                    reaction_removal_list.append(reaction)

                    # Find the products for this reaction and any other
                    # reactants which are required with the removed compound.
                    out_ps = network.reactions[reaction].products
                    ins = []

                    for product in network.reactions[reaction].reactants:
                        if product != compound:
                            ins.append(product)

                    # Process the tokens for the new reaction
                    new_reactants = ins + in_rs
                    new_products = outs + out_ps

                    new_reactants.sort()
                    new_products.sort()

                    # Create a new reaction between the incoming reaction's
                    # reactants and the outgoing reaction's products, including
                    # the other compounds required as reactants or products of
                    # both these reactions (which were not tagged for removal
                    LHS = ".".join(new_reactants)
                    RHS = ".".join(new_products)
                    re_str = f"{LHS}>>{RHS}"

                    rdkit_reaction = Chem.ReactionFromSmiles(re_str)
                    new_reaction = Reaction(rdkit_reaction)
                    new_reaction_list.append(new_reaction)

    network.add_reactions(new_reaction_list)

    # Remove reactions from the network and any compounds which are exclusively
    # produced by these reactions (if not already picked up in the compounds_removal_list.
    for reaction in set(reaction_removal_list):
        del network.reactions[reaction]
        for compound in network.compounds:
            for in_reaction in network.compounds[compound].product_of:
                if in_reaction == reaction:
                    network.compounds[compound].product_of.remove(reaction)
            for out_reaction in network.compounds[reaction].reactant_in:
                if out_reaction == reaction:
                    network.compounds[compound].reactant_in.remove(reaction)

    for compound in set(compounds_removal_list):
        del network.compounds[compound]
