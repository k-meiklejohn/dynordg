from ..core import  RiboNode
from networkx import DiGraph

class RiboGraph(DiGraph):
    """
    A NetworkX Digraph that only accepts RiboNode instances as nodes
    """

    def __init__(self, incoming_graph_data = None, **attr):
        super().__init__(incoming_graph_data, **attr)
        self.bulk_node = RiboNode((-1,-1))
        self.add_node(self.bulk_node)

    def add_node(self, node_for_adding: RiboNode, **attr):
        """
        Add a new RiboNode to the graph
        """
        if isinstance(node_for_adding, RiboNode):
            pass
        elif isinstance(node_for_adding, tuple):
            node_for_adding = RiboNode(node_for_adding)
        else:
            raise TypeError(f"Ribograph only accepts RiboNodes or RiboNode-like tuples, got {type(node_for_adding).__name__!r}")
        super().add_node(node_for_adding, **attr)
    
    
    def add_nodes_from(self, nodes_for_adding, **attr):
        """
        Add multiple RiboNodes to a graph
        """
        for node in nodes_for_adding:
            n = node[0] if isinstance(node, tuple) else node
            if not isinstance(n, RiboNode):
                if isinstance(n, tuple):
                    pass  # will be coerced by add_node
                else:
                    raise TypeError(f'RiboGraph only accepts RiboNodes or RiboNode-like tuples, got {type(n).__name__!r}')
        super().add_nodes_from(nodes_for_adding, **attr)
    
    
    def add_edge(self, u, v, **attr):
        """
        Add an edge consisting of 2 RiboNodes to a graph
        """
        if not isinstance(u, RiboNode):
            u = RiboNode(u)
        if not isinstance(v, RiboNode):
            v = RiboNode(v)

        if self.has_edge(u, v):
            existing = self.edges[u, v]
            merged = {}
            for key in set(existing) | set(attr):
                existing_val = existing.get(key)
                new_val = attr.get(key)
                if existing_val is None:
                    merged[key] = new_val
                elif new_val is None:
                    merged[key] = existing_val
                elif key.startswith('flux'): #this may change
                    merged[key] = existing_val + new_val
                else:
                    merged[key] = new_val  # overwrite non-flux attributes
            super().add_edge(u, v, **merged)
        else:
            super().add_edge(u, v, **attr)

    def horizontal_in_edge(self, node: RiboNode, data:bool|str|list=False):
        """
        Returns the in edge in the same phase as the node
        """
        if data:
            for u, v, data in self.in_edges(node, data=data):
                if u.phase == node.phase:
                    return (u,v,data)
        else:    
            for u, v in self.in_edges(node):
                if u.phase == node.phase:
                    return (u,v)
    def horizontal_out_edge(self, node: RiboNode, data:bool|str|list=False):
        """
        Returns the out edge in the same phase as the node
        """

        if data:
            for u, v, data in self.out_edges(node, data=data):
                if v.phase == node.phase:
                    return (u,v, data)
            
        else:    
            for u, v in self.in_edges(node):
                if u.phase == node.phase:
                    return (u,v)
