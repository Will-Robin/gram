class OutputProcess:
    """
    Class to store outputs as reactions into reaction network
    """

    def __init__(self, reaction_output_string: str):
        """

        Parameters
        ----------
        reaction_input_string: str
            A token for reaction output should follow the convention
            `SMILES>>#0`

        Attributes
        ----------
        token: string
            Token for the output.
        output_id: str
            ID for the output.
        output_compound: str
            The compound associated with the output.
        product_of: list[str]
            List of processes connected to the output.
        """

        self.token: str = reaction_output_string

        token_sides = reaction_output_string.split(">>")
        self.output_id: str = token_sides[1]
        self.output_compound: str = token_sides[0]

        self.product_of: list[str] = []
