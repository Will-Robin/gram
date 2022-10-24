class InputProcess:
    """
    Class to store inputs as reactions into reaction network
    """

    def __init__(self, reaction_input_string: str):
        """

        Parameters
        ----------
        reaction_input_string: str
            A token for reaction input should follow the convention
            `SMILES_#0>>SMILES`

        Attributes
        ----------
        input_id: str
            ID for the input.
        input_compound: str
            Compound parameters associated with the input.
        out: list[str]
            List of compounds arising from the input.
        """

        self.token: str = reaction_input_string

        token_sides = reaction_input_string.split(">>")
        self.input_id: str = token_sides[0]
        self.input_compound: list[str] = token_sides[1].split(".")
        self.out: list[str] = []
