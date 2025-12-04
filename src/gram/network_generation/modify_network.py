"""
These functions modify reaction networks.
"""
from gram.Classes import Network
from gram.Classes import InputProcess
from gram.Classes import OutputProcess


def add_flow_terms(network: Network, inputs: list[str]) -> None:
    """
    Add flow inputs and outputs into network.

    Parameters
    ----------
    network: Network
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
