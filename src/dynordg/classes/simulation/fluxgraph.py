from ..graph import RiboGraph
from .skeleton import RiboSkeleton
from ..core import RiboNode, State
from networkx import DiGraph, all_simple_paths, topological_sort
import warnings
from copy import deepcopy
from collections import defaultdict
from functools import cached_property


class RiboGraphFlux(RiboGraph):
    """
    A directed graph representing ribosomal flux through phase space, built from a TransitionMap.

    RiboGraphFlux models the movement of ribosomes along an mRNA transcript by propagating
    flux from a bulk cytoplasmic pool through a network of ribosomal phase-space nodes.
    Each edge carries flux values that account for initiation, elongation, termination,
    reinitiation, and decay during scanning and translation.

    The graph is constructed automatically on instantiation via construct(), which:
      - Removes edges below the weight cutoff from the transition map
      - Normalises entry flux across all initiation sites
      - Propagates flux iteratively through the graph, applying decay and branching
      - Collapses redundant intermediate nodes

    Parameters
    ----------
    transition_map : TransitionMap
        The transition map from which the flux graph is built. Encodes the probability
        of moving between ribosomal phase-space nodes.
    half_life_scanning : float or None
        The half-life (in nucleotides) of a scanning 40S subunit. Controls the rate of
        decay applied to flux along scanning edges (phase == 0). If None, no scanning
        decay is applied.
    half_life_translation : float or None
        The half-life (in nucleotides) of an elongating ribosome. Controls
        decay along translation edges (phase > 0). If None, no translational decay is applied.
    weight_cutoff : float
        Edges in the transition map with weight below this value are removed before
        flux propagation. Default is 0.0.
    reinitiation_half_life : float or None
        Controls the decay of reinitiation potential as a function of distance from the
        initiation node. Larger values allow reinitiation over greater distances.
        If None, reinitiation potential does not decay.
    ternary_complex_half_life : float or None
        Half-life governing the replenishment of ternary complex during scanning after
        a termination event. Used to scale initiation probability following reinitiation.
        If None, ternary complex is assumed to be fully available at all times.
    flux_cutoff : float
        Flux values below this threshold are not propagated further, pruning negligible
        branches from the graph. Default is 0.001.
    retention_limit : int or None
        Maximum number of consecutive retention events (ribosome returning to a scanning state on the same
        mRNA after termination) allowed per path. If None, retention is unlimited.
        Default is 1.

    Attributes
    ----------
    transitions : TransitionMap
        The underlying transition map used to build the graph.
    ribopaths : list of list of tuple
        All paths representing continuous 40S association with the mRNA, from the
        first loading node to the bulk node. Each path is a list of edge tuples (u, v).
    translons : list of list of tuple
        All sub-paths within ribopaths representing continuous 60S association
        (i.e., active translation). Each translon is a list of RiboNodes
        where both nodes have phase > 0.

    Edge Data
    ---------
    Each edge (u, v) in the graph carries the following data:
      flux_start : float
          Flux entering the edge at node u.
      flux_end : float
          Flux remaining at node v after decay along the edge.
      decay : float
          Flux lost along the edge due to ribosome drop-off.


    Methods
    -------
    flux_proportion(u, v) -> float
        Returns the total proportion of bulk flux that travels between nodes u and v,
        summed across all simple paths connecting them.
    flux_proportion_path(path) -> float
        Returns the proportion of bulk flux that travels along a specific path,
        given as an ordered list of RiboNodes.
    edge_weight(u, v) -> float
        Returns the fraction of flux leaving node u that is carried by the edge to v.
    edge_decay(u, v) -> float
        Returns the fraction of flux lost to ribosome drop-off along the edge from u to v.
    rein_decay(u, v) -> float
        Returns the decay factor applied to reinitiation potential between nodes u and v.
    ternary_complex_proportion(u, v) -> float
        Returns the fraction of ternary complex available at v, given prior scanning from u.
    add_transition(source, target, probability)
        Adds a single transition to the underlying map and reconstructs the flux graph.
    add_transitions_from(tbunch)
        Adds multiple transitions from an iterable of (source, target, weight) tuples
        and reconstructs the flux graph.

    Notes
    -----
    - Node phase conventions: phase == -1 represents the bulk cytoplasmic pool (off-mRNA),
      phase == 0 represents a scanning 40S subunit, and phase > 0 represents an elongating
      80S ribosome.
    - Flux is normalised so that the maximum edge flux equals 1.0.
    - Floating-point accumulation errors in complex graphs may cause small discrepancies
      between total inbound and outbound flux; a warning is raised if this exceeds the
      internal tolerance (flux_error = 1e-15).
    """

    def __init__(self, skeleton: RiboSkeleton, 
                 flux_cutoff = 0.0,
                 incoming_graph_data=None, 
                 **attr):
        

        super().__init__(incoming_graph_data, **attr)

        self.skeleton = skeleton
        if incoming_graph_data is not None:
            raise ValueError('Incoming graph data must be left empty, graph is calculated from skeleton')
        self.flux_cutoff = flux_cutoff
        self.flux_error = 0.000000000000001

        if self.skeleton:
            self.construct()   
    
    def construct(self):
        self._iterate_graph_topo()
        self._normalize_flux()
        self._is_valid()

    def _iterate_graph_topo(self):
        accumulated = defaultdict(float)
        accumulated[self.bulk_node] = 1
        

        # Copy to plain DiGraph to avoid subgraph_view instantiation issues
        dag = DiGraph(
            (u, v, d) for u, v, d in self.skeleton.edges(data=True)
            if v != self.bulk_node
        )

        for node in topological_sort(dag):
            flux = accumulated[node]
            if not flux:
                continue
            out_edge_weights = {(u,v):w for u,v,w in self.skeleton.out_edges(node, data='weight')}
            dispatchable = self._calc_flux_dispatch(out_edge_weights, flux)
            for edge, eflux in dispatchable.items():
                self.add_edge(edge[0], edge[1], flux=eflux)
                accumulated[edge[1]] += eflux


    def _calc_flux_dispatch(self, weight_dict: dict[tuple[RiboNode, RiboNode], float], flux: float) -> dict:
        if not weight_dict:
            return {}
        while weight_dict:
            # Find the single minimum edge
            min_edge = min(weight_dict, key=lambda e: weight_dict[e])

            if flux * weight_dict[min_edge] >= self.flux_cutoff:
                # All edges above cutoff
                return {edge: flux * w for edge, w in weight_dict.items()}

            # Remove minimum and redistribute its weight proportionally
            dead_weight = weight_dict.pop(min_edge)
            total_surviving = sum(weight_dict.values())

            for edge in weight_dict:
                weight_dict[edge] += dead_weight * (weight_dict[edge]/ total_surviving)
    
        return {}

                
    def _normalize_flux(self):
        flux_keys = ('flux', 'flux_start', 'flux_end', 'decay')
        max_flux = max(d[k] for k in flux_keys for _,_,d in self.edges(data=True) if k in d)
        factor = 1 / max_flux if max_flux > 1 else 1

        if factor == 1:
            return
        for _, _, data in self.edges(data=True):
            for k in flux_keys:
                if k in data:
                    data[k] *= factor
            


    def _is_valid(self):
        self._valid_in_out()

    @cached_property
    def ribopaths(self) -> list[list[RiboNode]]:
        """
        Returns a list of the paths with continued 40S association, each as a list of edge tuples.
        """
        paths = []
        for loading in self.successors(self.bulk_node):
            for path in all_simple_paths(self, loading, self.bulk_node):
                paths.append(path)
        return paths
    
    @cached_property
    def translons(self) -> list[list[tuple[RiboNode, RiboNode]]]:
        """
        Returns a list of all translons in the graph (continued 60S association) 
        """
        translon_list = []

        for path in self.ribopaths:
            translon = False
            current_translon = []

            for node in path:
                if translon:
                    if node.phase < 1:
                        translon = False
                        current_translon.append((source, prev_node.simple))
                        translon_list.append(current_translon)
                    
                    else:
                        if phase != node.phase:
                            phase = node.phase
                            current_translon.append((source, prev_node.simple))
                            source = node.simple

                elif node.phase > 0:
                    phase = node.phase
                    translon = True
                    source = node.simple
                prev_node = node

        return translon_list

    def flux_proportion(self, u:RiboNode, v:RiboNode) -> float:
        total_proportion = 0
        for path in all_simple_paths(self, u , v):
            total_proportion += self.flux_proportion_path(path)
        return total_proportion
    

    def flux_proportion_path(self, path: list[RiboNode]) -> float:
        if not path:
            raise ValueError('Cannot compute flux proportion of an empty path')

        # Total flux proportion arriving at path[0] from bulk_node
        entry_proportion = self.node_flux(path[0])

        # Proportion of that flux which travels along the given path
        path_proportion = 1.0
        for i in range(len(path) - 1):
            path_proportion *= self.edge_weight(path[i], path[i + 1])

        return entry_proportion * path_proportion

    def node_flux(self, nbunch: RiboNode) -> float:
        total_flux = 0.0
        for _, _, flux in self.out_edges(nbunch=nbunch, data='flux'):
            total_flux += flux
        return total_flux 

    def _valid_in_out(self):
        out_flux = 0
    
        for _,_, flux in self.out_edges(self.bulk_node, data='flux'):
            out_flux += flux

        in_flux = 0
        for node in self.successors(self.bulk_node):
            for _, _, flux in self.out_edges(node, data='flux'):
                in_flux += flux

            
        if abs(out_flux - in_flux) > self.flux_error:

            warnings.warn(f'Flux in: {in_flux} does not equal Flux out: {out_flux}, '
                          'this may occur due to accumulated errors in floating point numbers,' 
                          'especially in very complex graphs, ' 
                          'and can be ignored to your deisred level of accuracy.')
            
    def edge_weight(self, u: RiboNode,v: RiboNode):
        return self[u][v]['flux'] / self.node_flux(u)
    
    @property
    def simple(self):
        out = deepcopy(self)
        out.clear()
        out.bulk_node = self.bulk_node.simple
        for u, v, flux in self.edges(data='flux'):
            u:RiboNode
            v:RiboNode

            if u.simple == v.simple:
                continue
            elif v.phase == -1 and u.position != v.position and u.phase != -1:
                out.add_edge(u.simple, RiboNode(v.position, State(u.phase)), flux_start=flux, flux_end=0, decay=flux)
                out.add_edge(RiboNode(v.position, State(u.phase)), v.simple, flux_start=flux, flux_end=flux)
    
            else:
                out.add_edge(u.simple, v.simple,
                            flux_start=flux,
                            flux_end=flux)
                
        changed = True
        topo_nodes:list[RiboNode] = list(topological_sort(out.dag))

        while changed:
            changed = False

            for node in topo_nodes:
                

                if node.phase == -1 or node not in out:
                    continue

                in_edges = list(out.in_edges(node, data=True))
                if len(in_edges) != 1:
                    continue


                in_u, _, in_data = in_edges[0]
                if in_u.phase != node.phase:
                    continue

                out_edges = list(out.out_edges(node, data=True))
                h_out = [edge for edge in out_edges if edge[1].phase == node.phase]
                if not h_out:
                    continue
                else:
                    _, out_v, out_data = h_out[0]


                # Flux must be continuous across the pass-through node
                if abs(in_data['flux_end'] - out_data['flux_start']) > out.flux_error:
                    continue

                # Gather everything before touching the out
                in_flux_start = in_data['flux_start']
                in_decay      = in_data.get('decay', 0.0)
                out_decay     = out_data.get('decay', 0.0)
                total_decay   = in_decay + out_decay
                flux_end      = in_flux_start - total_decay
                drop_node = RiboNode(out_v.position, State(-1))

                # Read existing drop edge data before any mutation
                existing_drop = (out[out_v][drop_node].copy()
                                if out.has_edge(out_v, drop_node) else {})

                # Now mutate: remove the intermediate node (drops both its edges)
                
                out.remove_node(node)

                # Reconnect: add_edge accumulates, so safe to call even if edge exists

                out.add_edge(in_u, out_v,
                            flux_start = in_flux_start,
                            flux_end   = flux_end,
                            **({'decay': total_decay} if total_decay > 0 else {}))

                # Update drop edge if there was decay
                if total_decay > 0:
                    if out.has_edge(out_v, drop_node):
                        out.remove_edge(out_v, drop_node)
                    if out.has_edge(drop_node, out.bulk_node):
                        out.remove_edge(drop_node, out.bulk_node)

                    # Preserve any drop flux that existed independently of this merge
                    prior_drop = existing_drop.get('flux_start', 0.0)
                    # Subtract the out_decay portion we already accounted for —
                    # the remainder was from other paths, keep it
                    independent_drop = max(prior_drop - out_decay, 0.0)
                    new_drop = total_decay + independent_drop

                    out.add_edge(out_v, drop_node,
                                flux_start = new_drop,
                                flux_end   = new_drop)
                    out.add_edge(drop_node, out.bulk_node,
                                flux_start = new_drop,
                                flux_end   = new_drop)

                changed = True
                break  # restart topo walk with updated graph
        return(out)

    
    def prune_recycle_edges(self) -> None:
        """
        Remove bulk-to-bulk edges and any nodes they leave isolated.

        Edges where both u and v have phase == -1 are recycling arcs internal
        to the bulk pool.  They carry no mRNA-positional information and
        produce degenerate geometry, so they are stripped before layout runs.
        """
        dead = []

        dead = [(u, v) for u, v in self.edges
                if u.phase == -1 and v.phase == -1]
        self.remove_edges_from(dead)
        isolated = [n for n, d in self.degree() if d < 1]
        self.remove_nodes_from(isolated)
