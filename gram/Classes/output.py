class Output:
    """
    Class to store inputs into reaction network
    """

    def __init__(self, output: str):
        """

        Parameters
        ----------
        output: str
            A token for reaction input should follow the convention `SMILES_#0`

        Attributes
        ----------
        token: str

        product_of: list[str]
            List of Reaction Output processes connected to the output.
        """

        self.token: str = output
        self.product_of: list[str] = []
