from gram.Classes import InputProcess


def input_process_from_string(input_process: str) -> InputProcess:
    """
    Create an InputProcess object from a string.

    Parameters
    ----------
    input_process: str

    Returns
    -------
    InputProcess
    """
    return InputProcess(input_process)
