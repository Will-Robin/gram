from gram.Classes import Input


def input_from_string(input_str: str) -> Input:
    """
    Create an Input object from a string.

    Parameters
    ----------
    input_str: str

    Returns
    -------
    Input
    """

    return Input(input_str)
