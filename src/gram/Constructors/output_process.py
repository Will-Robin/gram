from gram.Classes import OutputProcess


def output_process_from_string(output_process: str) -> OutputProcess:
    """
    Create an OutputProcess object from a string.

    Parameters
    ----------
    output_process: str

    Returns
    -------
    OutputProcess
    """
    return OutputProcess(output_process)
