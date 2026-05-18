from .core import Transcript, DEFAULT_BEHAVIOUR
from .graph import RiboSkeleton, FluxGraph
from .render import RiboGraphVis

def quick_graph(sequence: str):
    """
    Quickly generate a flux graph from a nucleotide sequence using default settings
    """
    tcript = Transcript(sequence)
    tcript.auto_stop_starts(reinitiation_prob=1)
    skeleton = RiboSkeleton(tcript.pipes, behaviours=DEFAULT_BEHAVIOUR)
    return FluxGraph(skeleton)

def quick_plot(sequence: str, log_scale = 3):
    """
    Quickly generate a plot from a nucleotide sequence
    Run result.show() to display.

    Example:

    >>> from dynordg import quick_plot
    >>> plot = quick_plot("AUCGAUCUGCACUGACGCGUACGUGAUCGUAGC")
    >>> plot.show()
    """
    graph = quick_graph(sequence)
    return RiboGraphVis(graph, log_scale=log_scale)