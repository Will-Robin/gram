from .reaction import Reaction
from .compound import Compound
from .substructure import Substructure
from .reaction_rule import ReactionRule
from gram.chemoinformatics import substructure_match as substr


class RuleNetwork:
    """
    A network designed to show the relationship between functional group
    transformations.

    The network has three kinds of node: reaction rules, compounds, and
    substructures. Substructures connect reactions to compounds. Compounds and
    reactions can only link to substructures.
    """

    def __init__(self, rules, name, description):
        """
        Parameters
        ----------
        rules: list[gram.Classes.ReactionRule]
            List of reaction rules to create the network.
        name: str
            A name for the network.
        description: str
            A description of the network.

        Attributes
        ----------
        name: str
            A name for the network.
        description: str
            A description of the network.
        rules: dict
            A dictionary containing the ReactionRule objects keyed
            by their reaction SMILES: {reaction SMARTS: ReactionRule}
        substructures: dict
            A dictionary containing the Substructure objects keyed
            by their SMARTS: {SMARTS: Substructure}
        reactions: dict
            A dictionary containing the Reaction objects keyed
            by their reaction SMILES: {reactionSMILES: Reaction}
        compounds: dict
            A dictionary containing the Compound objects keyed
            by their SMILES: {SMILES: Compound}
        """

        self.name: str = name

        self.description: str = description

        self.reaction_rules: dict[str, ReactionRule] = dict()

        self.substructures: dict[str, Substructure] = dict()

        self.compounds: dict[str, Compound] = dict()

        self.reactions: dict[str, Reaction] = dict()

        if len(rules) == 0:
            pass
        else:
            self.add_rules(rules)

    def add_compound(self, compound: Compound):
        """
        Add a compound to the network.

        The compound will only be added to the SubstructureNetwork if it has at
        least one substructure match to one of the substructures already in the
        network.

        Parameters
        ----------
        compound: gram.Classes.Compound

        Returns
        -------
        None
        """

        # Only perform insertion if the compound is not already in the
        # SubstructureNetwork
        if compound.smiles not in self.compounds:
            # find the subtructures in the SubstructureNetwork which match the
            # compound

            add_to_net = False  # changes to True if a match is found
            for sub in self.substructures:
                substructure = self.substructures[sub]
                if substr.has_substructure_match(compound, substructure):
                    compound.reactive_substructures.append(sub)
                    substructure.matching_compounds.append(compound.smiles)
                    if not add_to_net:
                        add_to_net = True

            if add_to_net:
                self.compounds[compound.smiles] = compound

    def add_compounds(self, compounds: list[Compound]):
        """

        Parameters
        ----------
        compounds: list[gram.Classes.Compound]
            compounds to be added

        Returns
        -------
        None
        """
        for compound in compounds:
            self.add_compound(compound)

    def remove_compound(self, compound: Compound):
        """
        Remove a compound from the substructure network.


        Parameters
        ----------
        compound: gram.Classes.Compound

        Returns
        -------
        None
        """

        compound_token = compound.smiles
        # remove the compound connection to substructures
        for substruct in compound.reactive_substructures:
            substructure = self.substructures[substruct]
            if compound_token in substructure.matching_compounds:
                substructure.matching_compounds.remove(compound_token)

        del self.compounds[compound_token]

    def remove_compounds(self, compounds: list[Compound]):
        """
        Remove list of compounds from the SubstructureNetwork.

        Parameters
        ----------
        compounds: list[gram.Classes.Compound]
            compounds to be removed

        Returns
        -------
        None
        """

        for compound in compounds:
            self.remove_compound(compound)

    def add_reaction_rule(self, rule: ReactionRule):
        """
        Add a reaction rule to the network.

        Parameters
        ----------
        rule: gram.Classes.ReactionRule

        Returns
        -------
        None
        """

        rxn_key = rule.reaction_smarts

        if rxn_key not in self.reaction_rules:
            self.reaction_rules[rxn_key] = rule

        for rxn_substr in rule.reactant_substructures:
            # connect to reaction
            if rxn_substr not in self.substructures:
                self.substructures[rxn_substr.smarts] = rxn_substr

            self.substructures[rxn_substr.smarts].reaction_participations.append(
                rxn_key
            )

            working_substruct = self.substructures[rxn_substr.smarts]

            # connect to compounds
            for comp in self.compounds:
                compound = self.compounds[comp].mol
                if substr.has_substructure_match(compound, working_substruct):
                    self.compounds[comp].reactive_substructures.append(
                        rxn_substr.smarts
                    )
                    self.substructures[rxn_substr.smarts].matching_compounds.append(
                        comp
                    )

        for p_substruct in rule.product_substructures:
            # connect to reaction
            if p_substruct not in self.substructures:
                self.substructures[p_substruct.smarts] = p_substruct

            self.substructures[p_substruct.smarts].producing_reactions.append(rxn_key)

            working_substruct = self.substructures[p_substruct.smarts]

            # connect to compounds
            for comp in self.compounds:
                compound = self.compounds[comp]
                if substr.has_substructure_match(compound, working_substruct):
                    self.compounds[comp].reactive_substructures.append(
                        p_substruct.smarts
                    )
                    self.substructures[p_substruct.smarts].matching_compounds.append(
                        comp
                    )

    def add_reaction_rules(self, rules: list[ReactionRule]):
        """
        Add reaction rules to the network.

        Parameters
        ----------
        rules: list[gram.Classes.ReactionRule]

        Returns
        -------
        None
        """

        for t in rules:
            if t.reaction_smarts in self.reaction_rules:
                pass
            else:
                self.add_reaction_rule(t)

    def remove_reaction_rule(self, rule: ReactionRule):
        """
        Remove a reaction rule from a SubstructureNetwork.

        Parameters
        ----------
        rule: gram.Classes.ReactionRule

        Returns
        -------
        None
        """

        # Disconnect the reaction rule from substructures
        for r_substr in rule.reactant_substructures:
            substructure = self.substructures[r_substr.smarts]
            rxn_partic = substructure.reaction_participations
            if rule.reaction_smarts in rxn_partic:
                rxn_partic.remove(rule.reaction_smarts)

        for r_substr in rule.reactant_substructures:
            substructure = self.substructures[r_substr.smarts]
            prod_rxns = substructure.producing_reactions
            if rule.reaction_smarts in prod_rxns:
                prod_rxns.remove(rule.reaction_smarts)

        # Check for disconnected substructures
        remove_substructs: list[str] = []
        for sub in self.substructures:
            substructure = self.substructures[sub]
            rxn_partic = substructure.reaction_participations
            prod_rxns = substructure.producing_reactions
            if len(prod_rxns) == 0 and len(rxn_partic) == 0:
                remove_substructs.append(sub)

        # disconnect the substructure from any matching compounds.
        for sub in remove_substructs:
            for comp in self.substructures[sub].matching_compounds:
                self.compounds[comp].reactive_substructures.remove(sub)
            del self.substructures[sub]

        # Check for disconnected compounds
        remove_compounds = []
        for comp in self.compounds:
            compound = self.compounds[comp]
            if len(compound.reactive_substructures) == 0:
                remove_compounds.append(comp)

        for comp in remove_compounds:
            del self.compounds[comp]

        # Finally, remove the reaction rule
        del self.reaction_rules[rule.reaction_smarts]

    def remove_reaction_rules(self, rules: list[ReactionRule]):
        """
        Remove a list of ReactionTemplates from the SubstructureNetwork.

        Parameters
        ----------
        rules: list[gram.Classes.ReactionRule]

        Returns
        -------
        None
        """

        for temp in rules:
            self.remove_reaction_rule(temp)

    def add_reaction(self, reaction: Reaction):
        """
        Add a reaction to the network.

        Parameters
        ----------
        reaction: gram.Classes.Reaction

        Returns
        -------
        None
        """

        r_key = reaction.reaction_rule.reaction_smarts

        self.reaction_rules[r_key] = reaction.reaction_rule

        for r_subst in reaction.reaction_rule.reactant_substructures:
            # connect to reaction
            if r_subst.smarts in self.substructures:
                pass
            else:
                self.substructures[r_subst.smarts] = r_subst

            self.substructures[r_subst.smarts].reaction_participations.append(r_key)

            working_substruct = self.substructures[r_subst.smarts]

            # connect to compound
            for reac in reaction.reactants:
                if reac in self.compounds:
                    pass
                else:
                    self.compounds[reac] = Compound(reac)

                compound = self.compounds[reac]

                if substr.has_substructure_match(compound, working_substruct):
                    compound.reactive_substructures.append(r_subst.smarts)
                    self.substructures[r_subst.smarts].matching_compounds.append(reac)
                else:
                    pass

        for p_substruct in reaction.reaction_rule.product_substructures:
            # connect to reaction
            if p_substruct.smarts in self.substructures:
                pass
            else:
                self.substructures[p_substruct.smarts] = p_substruct

            self.substructures[p_substruct.smarts].producing_reactions.append(r_key)

            working_substruct = self.substructures[p_substruct.smarts]

            # connect to compound
            for prod in reaction.products:
                if prod in self.compounds:
                    pass
                else:
                    self.compounds[prod] = Compound(prod)

                compound = self.compounds[prod]

                if substr.has_substructure_match(compound, working_substruct):
                    self.compounds[prod].reactive_substructures.append(
                        p_substruct.smarts
                    )
                    self.substructures[p_substruct.smarts].matching_compounds.append(
                        prod
                    )
                else:
                    pass

    def add_reactions(self, reactions: list[Reaction]):
        """
        Add reactions to the network.

        Parameters
        ----------
        reactions: list[gram.Classes.Reaction]

        Returns
        -------
        None
        """

        for r in reactions:
            if r.reaction_rule.name == "none":
                pass
            elif r.reaction_rule.reaction_smarts in self.reaction_rules:
                pass
            else:
                self.add_reaction(r)
