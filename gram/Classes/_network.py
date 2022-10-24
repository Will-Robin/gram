from ._input import Input
from ._output import Output
from ._reaction import Reaction
from ._compound import Compound
from ._inputProcess import InputProcess
from ._outputProcess import OutputProcess


class Network:
    """
    An object which stored compounds and reactions and the connections
    between them.
    """

    def __init__(self, reactions: list[Reaction], name: str, description: str):
        """
        The Network object is initialised with a list of Reaction objects. If
        the list is empty, then the network is initialised as an empty network.

        Parameters
        ----------
        reactions: list[Reaction]
            List of reactions to create the network.
        name: str
            A label name for the network.
        description: str
            A description for the network.

        Attributes
        ----------
        name: str
            A name for the network.
        description: str
            A description of the network.
        reactions: dict
            A dictionary containing the Reaction objects keyed
            by their reaction SMILES: {reactionSMILES: Reaction}
        compounds: dict
            A dictionary containing the Compound objects keyed
            by their SMILES: {SMILES: Compound}
        inputs: dict
            A dictionary containing the Input objects keyed
            by their tokens: {token: Input}
        outputs: dict
            A dictionary containing the Output objects keyed
            by their tokens: {token: NetworkOutput}
        input_processes: dict
            A dictionary containing the InputProcess objects keyed
            by their tokens: {token: ReactionInput}
        output_processes: dict
            A dictionary containing the OutputProcess objects keyed
            by their tokens: {token: ReactionOutput}
        """

        self.name: str = name

        self.description: str = description

        self.reactions: dict[str, Reaction] = dict()

        self.compounds: dict[str, Compound] = dict()

        self.inputs: dict[str, Input] = dict()

        self.outputs: dict[str, Output] = dict()

        self.input_processes: dict[str, InputProcess] = dict()

        self.output_processes: dict[str, OutputProcess] = dict()

        if len(reactions) == 0:
            pass
        else:
            self.add_reactions(reactions)

    def add_compound(self, compound: Compound):
        """
        Add a compound to the network.

        Parameters
        ----------
        compound: Compound

        Returns
        -------
        None
        """

        if compound.smiles not in self.compounds:
            self.compounds[compound.smiles] = compound

    def add_compounds(self, compounds: list[Compound]):
        """
        Add list of Compound objects to the Network

        Parameters
        ----------
        compounds: list[Compound]
            compounds to be added

        Returns
        -------
        None
        """

        for compound in compounds:
            self.add_compound(compound)

    def remove_compound(self, compound: Compound):
        """
        Remove a compound and the reactions in which it is involved
        from the network.

        Parameters
        ----------
        compound: Compound
            compound to be removed

        Returns
        -------
        None
        """

        remove_reactions = []
        remove_reactions.extend(self.compounds[compound.smiles].product_of)
        remove_reactions.extend(self.compounds[compound.smiles].reactant_in)

        remove_reactions = list(set(remove_reactions))

        self.remove_reactions(remove_reactions)

        del self.compounds[compound.smiles]

    def remove_compounds(self, compounds: list[Compound]):
        """
        Remove list of compounds and the reactions in which they are involved
        from the network.

        Parameters
        ----------
        compounds: list[Compound]
            compounds to be removed

        Returns
        -------
        None
        """

        remove_reactions = []
        for compound in compounds:
            remove_reactions.extend(self.compounds[compound.smiles].product_of)
            remove_reactions.extend(self.compounds[compound.smiles].reactant_in)

        remove_reactions = list(set(remove_reactions))

        self.remove_reactions(remove_reactions)

        for compound in compounds:
            del self.compounds[compound.smiles]

    def add_reaction(self, reaction: Reaction):
        """
        Adds a reaction and its associated reactants and products into the
        Network.

        Parameters
        ----------
        reaction: Reaction
            reaction to be added

        Returns
        -------
        None
        """

        if reaction.reaction_smiles not in self.reactions:
            reaction_smiles = reaction.reaction_smiles

            self.reactions[reaction_smiles] = reaction

            for reactant in reaction.reactants:
                if reactant not in self.compounds:
                    new_compound = Compound(reactant)
                    self.add_compound(new_compound)
                    reactant = new_compound.smiles

                if reaction_smiles not in self.compounds[reactant].reactant_in:
                    self.compounds[reactant].reactant_in.append(reaction_smiles)

            for product in reaction.products:
                if product not in self.compounds:
                    new_compound = Compound(product)
                    self.add_compound(new_compound)
                    product = new_compound.smiles

                if reaction_smiles not in self.compounds[product].product_of:
                    self.compounds[product].product_of.append(reaction_smiles)

    def add_reactions(self, reactions: list[Reaction]):
        """
        Use the standardised strings information in the Reaction objects
        to build them into the Network.

        Parameters
        ----------
        reactions: list[Reaction]
            reactions to be added to the Network

        Returns
        -------
        None
        """

        for reaction in reactions:
            self.add_reaction(reaction)

    def remove_reactions(self, remove_reactions: list[str]):
        """
        Remove a list of reactions from the Network.

        Parameters
        ----------
        remove_reactions: list[str]
            reactions to be removed (either in reaction SMILES format or as
            Reaction).

        Returns
        -------
        None
        """

        for reaction in remove_reactions:

            for reactant in set(self.reactions[reaction].reactants):
                self.compounds[reactant].reactant_in.remove(reaction)
            for product in set(self.reactions[reaction].products):
                self.compounds[product].product_of.remove(reaction)

            del self.reactions[reaction]

    def add_input_process(self, input_addition: InputProcess):
        """
        Add a InputProcess to the Network

        Parameters
        ----------
        input_addition: InputProcess
            Input to be added.

        Returns
        -------
        None
        """

        if input_addition.token not in self.input_processes:
            for compound in input_addition.input_compound:
                if compound in self.compounds:
                    self.compounds[compound].product_of.append(input_addition.token)

                    self.input_processes[input_addition.token] = input_addition

                    if input_addition.input_id not in self.inputs:
                        self.inputs[input_addition.input_id] = Input(
                            input_addition.input_id
                        )

                    self.inputs[input_addition.input_id].out.append(
                        input_addition.token
                    )

    def add_input_processes(self, inputs: list[InputProcess]):
        """
        For adding a list of InputProcess objects to the out

        Parameters
        ----------
        inputs: list[InputProcess]
            Inputs to be added to the network.

        Returns
        -------
        None
        """

        for i in inputs:
            self.add_input_process(i)

    def add_output_process(self, output: OutputProcess):
        """
        Add a OutputProcess to the Network

        Parameters
        ----------
        output OutputProcess
            Output to be added

        Returns
        -------
        None
        """

        if output.output_compound in self.compounds:

            self.output_processes[output.token] = output

            self.compounds[output.output_compound].reactant_in.append(output.token)

            if output.output_id not in self.outputs:
                self.outputs[output.output_id] = Output(output.output_id)

            self.outputs[output.output_id].product_of.append(output.token)

        else:
            # The output compound is not in the network, so cannot be an output
            pass

    def add_output_processes(self, outputs: list[OutputProcess]):
        """
        For adding a list of OutputProcess objects to the network.

        Parameters
        ----------
        outputs: list[OutputProcess]
            Outputs to be added to the network

        Returns
        -------
        None
        """

        for out in outputs:
            self.add_output_process(out)

    def get_reaction(self, reaction: str):
        """
        Get a reaction using a key.

        Parameters
        ----------
        reaction: str
            Key in self.reactions

        Returns
        -------
        Reaction or None
        """

        if reaction in self.reactions:
            return self.reactions[reaction]

        print("Reaction not found in Network")

        return None

    def get_reactants(self, reaction: str):
        """
        Conveniently get the reactants of a reaction

        Parameters
        ----------
        reaction: str
            Key in self.reactions

        Returns
        -------
        None
        """

        reaction_entry = self.get_reaction(reaction)
        if reaction_entry is None:
            print("Reaction not found in Network")
            return None

        return reaction_entry.reactants

    def get_products(self, reaction: str):
        """
        Conveniently get the products of a reaction

        Parameters
        ----------
        reaction: str
            Key in self.reactions

        Returns
        -------
        None
        """

        reaction_entry = self.get_reaction(reaction)
        if reaction_entry is None:
            print("Reaction not found in Network")
            return None

        return reaction_entry.products

    def get_reaction_template(self, reaction: str):
        """
        Conveniently get the ReactionTemplate of a reaction

        Parameters
        ----------
        reaction: str
            Key in self.reactions

        Returns
        -------
        None
        """

        reaction_entry = self.get_reaction(reaction)
        if reaction_entry is None:
            print("Reaction not found in Network")
            return None

        return reaction_entry.reaction_template

    def get_reaction_SMARTS(self, reaction: str):
        """
        Conveniently get the Reaction SMARTS of a reaction

        Parameters
        ----------
        reaction: str
            Key in self.reactions

        Returns
        -------
        None
        """

        reaction_entry = self.get_reaction(reaction)

        if reaction_entry is None:
            print("Reaction not found in Network")
            return None

        if reaction_entry.reaction_template is not None:
            return reaction_entry.reaction_template.reaction_SMARTS

        return None

    def get_reaction_name(self, reaction: str):
        """
        Conveniently get the Name of a reaction

        Parameters
        ----------
        reaction: str
            Key in self.reactions

        Returns
        -------
        None
        """

        reaction_entry = self.get_reaction(reaction)

        if reaction_entry is None:
            print("Reaction not found in Network")
            return None

        if reaction_entry.reaction_template is not None:
            return reaction_entry.reaction_template.name

        return None
