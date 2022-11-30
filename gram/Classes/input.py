class Input:
    """
    Class to store inputs into reaction network
    """

    def __init__(self, net_input: str):
        """

        Parameters
        ----------
        net_input: str
            A token for reaction input. It should follow the convention
            `SMILES_#0`.

        Attributes
        ----------
        token: str
        out: list[str]
            List of reaction input tokens which this NetworkInput feeds into.
        """

        self.token = net_input
        self.out: list[str] = []
