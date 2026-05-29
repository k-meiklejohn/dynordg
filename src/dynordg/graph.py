from .core import  RiboNode, State, Pipe, FactorBehaviour, Transition
from collections import defaultdict
import networkx as nx
from typing import Literal
from dataclasses import dataclass
import networkx as nx
import warnings
from functools import cached_property
from copy import deepcopy


class RiboGraph(nx.DiGraph):
    """
    A NetworkX Digraph that only accepts RiboNode instances as nodes

    This is the base class for all graph objects in dynordg,
    it ensurse that any node added is properly marked within ribosomal phase space.

    It inherits all original methods from networkx DiGraph except:
    add_node       -- now only accepts RiboNodes
    add_nodes_from -- uses the overwritten add_node
    add_edge       -- tries to coerce nodes to RiboNode, any attribute
                      on an edge beginning with 'flux' or 'decay'
                      is additive, all others overwrite 

    Attributes:
    
    dag -- returns a new graph removing any edge leading to the bulk node, usually 
           returning a directed acyclic graph
    bulk_node -- alias for (Pos:-1, Phase:-1) in ribosomal phase space
    
    Methods:

    remove_weakly_connected -- iteratively removes every node with no in edges
                               essentially removes any node not part of a cycle
    """

    def __init__(self, incoming_graph_data = None, **attr):
        super().__init__(incoming_graph_data, **attr)
        self.bulk_node = RiboNode(-1,State(-1))
        self.add_node(self.bulk_node)



    def add_node(self, node_for_adding: RiboNode, **attr):
        """
        Add a new RiboNode to the graph
        """
        if not isinstance(node_for_adding, RiboNode):
            raise TypeError(f"Ribograph only accepts RiboNodes or RiboNode-like tuples, got {type(node_for_adding).__name__!r}")
        super().add_node(node_for_adding, **attr)
    
    
    def add_nodes_from(self, nodes_for_adding, **attr):
        """
        Add multiple RiboNodes to a graph
        """
        for node in nodes_for_adding:
            if not isinstance(node, RiboNode):
                if isinstance(node, tuple):
                    pass  # will be coerced by add_node
                else:
                    raise TypeError(f'RiboGraph only accepts RiboNodes or RiboNode-like tuples, got {type(node).__name__!r}')
        super().add_nodes_from(nodes_for_adding, **attr)
    
    
    def add_edge(self, u:RiboNode, v:RiboNode, **attr):
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
                elif key.startswith(('flux', 'decay')): #this may change
                    merged[key] = existing_val + new_val
                else:
                    merged[key] = new_val  # overwrite non-flux attributes
            super().add_edge(u, v, **merged)
        else:
            super().add_edge(u, v, **attr)
                
    def remove_weakly_connected(self):
        """Remove all nodes not part of a cycle"""

        while True:
            zero_indegree = [n for n, d in self.in_degree() if d == 0]

            if not zero_indegree:
                break

            self.remove_nodes_from(zero_indegree)

    @property
    def dag(self):
        return nx.DiGraph((u, v) for u, v in self.edges() if v != self.bulk_node)


@dataclass(frozen=True)
class PipeIndexEntry:
    pipe: Pipe
    kind: Literal["entry", "exit"]

    @property
    def position(self) -> int|None:
        if self.kind == "entry":
            return self.pipe.input_position
        return self.pipe.output_position
    


class RiboSkeleton(RiboGraph):
    """
    RiboGraph subclass that holds information about all possible choices of a ribosome

    All edges consist of two RiboNodes and a weight
    
    The graph is cyclical through the 'bulk node' (Pos:-1,Phase:-1) by convention

    The parent classes for adding edges are not implemented and instead add_transition is used
    to ensure that edges of the correct form are added to the graph. 

    Parameters:
    pipelist: list[Pipe]
        The pipes from which the choices are to be constructed
    behaviours: list[FactorBehaviour]
        The behaviours which apply to ribosomes as they move along the transcript

    """

    def __init__(
        self,
        pipelist: list[Pipe] | None = None,
        behaviours: list[FactorBehaviour] | None = None,
        **attr,
    ):
        incoming_graph_data=None
        super().__init__(incoming_graph_data, **attr)
        self.pipelist = pipelist
        
        self.stack: list[tuple[RiboNode, float]] = []
        self.behaviours: list[FactorBehaviour] | None = sorted(behaviours) if behaviours else None
        self.visited: set[RiboNode] = set()
        self._pipe_index = defaultdict(list)
        self.expanded: set[RiboNode] = set()
        self.queued: set[RiboNode] = set()
        if self.pipelist:
            self.pipelist.sort(reverse=True)
            for pipe in self.pipelist:

                phases = pipe.phase_condition.phases

                entry = PipeIndexEntry(pipe, "entry")
                exit_ = PipeIndexEntry(pipe, "exit")

                if isinstance(phases, set):
                    for p in phases:
                        self._pipe_index[p].append(entry)
                else:
                    self._pipe_index[phases].append(entry)

                self._pipe_index[pipe.output_phase].append(exit_)
            self._construct()

    def _construct(self):
        if not self.pipelist:
            raise RuntimeError('No Pipelist to construct from')
        for pipe in self.pipelist:
            if pipe.phase_condition.phases == -1:
                factors = pipe.subphase_condition.required or None
                f_dict = {factor:True for factor in factors} if factors is not None else None
                source = RiboNode(pipe.input_position, State(-1, **f_dict)) if f_dict \
                    else RiboNode(pipe.input_position, State(-1))
                target = pipe.target(source.state)
                self.add_transition(Transition(self.bulk_node, source, 1))
                self.add_transition(Transition(source, target, pipe.probability))
                self._build_graph(target)
        self._is_valid()

    def _build_graph(self, starting_node: RiboNode):
        self.stack.append((starting_node, 1))
        while self.stack:
            node, weight = self.stack.pop()
            self.queued.discard(node)
            self.visited.add(node)

            pipe_indexes = self.next_pipes(node)
            
            if not pipe_indexes:
                raise RuntimeError(
                    f"Continue flux at node: {node} runs to infinity (no downstream pipes)"
                )
            pipes = [idx.pipe for idx in pipe_indexes]

            if pipe_indexes[0].position:
                node_at_pipe = RiboNode(pipe_indexes[0].position, node.state)
            else:
                raise RuntimeError("Pipe found without position")
            remaining_weight = weight

            if self.behaviours:
                for behaviour in self.behaviours:
                    if behaviour.applies(node, node_at_pipe):
                        f_trans = behaviour.transitions(node, node_at_pipe, remaining_weight)
                        if f_trans.weight > 0:
                            remaining_weight -= f_trans.weight
                            self.add_transition(f_trans)
                            self._expand_node(f_trans.target, pipes)

            if remaining_weight > 0:
                self.add_transition(Transition(node, node_at_pipe, remaining_weight))
                left_over = self._pipe_out(node_at_pipe, pipes, weight_scale=1)
                if left_over > 0:
                    self._add_to_stack(node_at_pipe, left_over)

    def _expand_node(self, node: RiboNode, pipes: list[Pipe]):
        if node in self.expanded:
            return

        self.expanded.add(node)

        if node.phase == -1:
            self.add_transition(Transition(node, self.bulk_node, 1))
            return

        unpiped = self._pipe_out(node, pipes, weight_scale=1)

        if unpiped > 0:
            self._add_to_stack(node, unpiped)

    def _pipe_out(self, node: RiboNode, pipes: list[Pipe], weight_scale: float) -> float:
        """
        Add transitions from node through any pipes whose enter() condition is met.
        Returns the remaining (unpiped) weight fraction.
        """
        unpiped = weight_scale
        for pipe in pipes:
            if not pipe.enter(node.state):
                continue
            pipe_target = pipe.target(node.state)
            self.add_transition(
                Transition(node, pipe_target, pipe.probability * unpiped)
            )
            flow = pipe.probability * unpiped
            unpiped -= flow
            if pipe_target.phase == -1:
                self.add_transition(Transition(pipe_target, self.bulk_node, 1))
            else:
                self._add_to_stack(pipe_target, 1)
        return unpiped

    def _add_to_stack(self, node, weight):
        if node not in self.visited and node not in self.queued:
            self.stack.append((node, weight))
            self.queued.add(node)

    def next_pipes(self, node: RiboNode) -> list[PipeIndexEntry] | None:
        candidates = self._pipe_index.get(node.phase)

        if not candidates:
            return None

        min_pos = None
        result = []

        for entry in candidates:
            pos = entry.position

            if pos <= node.position:
                continue

            if min_pos is None or pos < min_pos:
                min_pos = pos
                result = [entry]

            elif pos == min_pos:
                result.append(entry)

        return result or None
    

    def add_transition(self, transition: Transition):
        if isinstance(transition, Transition):
            if transition.weight == 0:
                return
            super().add_edge(transition.source, transition.target, weight=transition.weight)
        else:
            raise TypeError(f"add_transition requires type 'Transition', got {type(transition).__name__}")

    def add_transitions_from(self, transitions):
        for t in transitions:
            self.add_transition(t)



    def _is_valid_weight(self):

        for node in self.nodes:

            if node.phase == -1:
                continue

            if not any(True for _ in self.successors(node)):
                continue

            total_weight = 0

            for _, _, w in self.out_edges(node, data='weight'):
                total_weight += w

            if total_weight > 1:
                raise RuntimeError(f"Weight from node: {node} exceeds 1")

    def _is_valid(self):
        self._is_valid_weight()

    def smooth_all_weights(self, by_factor:float|None=None):
        """Redistribute the weight from all nodes more evenly"""
        for node in self.nodes:
            node:RiboNode
            self.smooth_node_weight(node=node, by_factor=by_factor)
        self._is_valid()
            

    def smooth_node_weight(self, node:RiboNode, by_factor:float|None=None):
        """Redistribute weight more evenly from a node"""
        if by_factor != None:
            if by_factor < 0:
                raise ValueError(f"by_factor must be greater than or equal to 0, got {by_factor}")
            out_num = self.out_degree(node)
            target = 1/out_num
        
            for u,v,weight in self.out_edges(node, data='weight'):
                distance = abs(target - weight)
                distance = - distance if weight > target else distance
                self[u][v]['weight'] = weight + distance * (1 - 1/by_factor)
        else:
            for u,v,weight in self.out_edges(node, data='weight'):
                self[u][v]['weight'] = target



class FluxGraph(RiboGraph):
    """
    A directed graph representing ribosomal flux through phase space, built from a TransitionMap.

    RiboGraphFlux models the movement of ribosomes along an mRNA transcript by propagating
    flux from a bulk cytoplasmic pool through a network of ribosomal phase-space nodes base on
    a RiboSkeleton

    The graph is constructed automatically on instantiation via construct(), which:
        -Propagates flux along the skeleton nodes, ignoring edges that dont meet the flux cutoff


    Parameters
    ----------
    skeleton : RiboSkeleton
        The map which provides the possible paths of the ribosomes
    flux_cutoff : float
        Flux values below this threshold are not propagated further, pruning negligible
        branches from the graph. Default is 0.0.


    Attributes
    ----------
    skeleton : RiboSkeleton
        The underlying map used to build the graph.
    ribopaths : list of list of RiboNode
        All paths representing continuous 40S association with the mRNA, from the
        first loading node to the bulk node. Each path is a list of RiboNodes.
    translons : list of list of tuple
        All sub-paths within ribopaths representing continuous 60S association
        (i.e., active translation). Each translon is a list of RiboNode tuples
        where both nodes have phase > 0.


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

    def __init__(self, skeleton: RiboSkeleton| None = None, 
                 flux_cutoff = 0.0,
                 **attr):
        

        super().__init__(None, **attr)
        self.skeleton = skeleton
        self.flux_cutoff = flux_cutoff
        self.flux_error = 0.000000000000001

        if self.skeleton:
            self.construct()   
    
    def construct(self):
        self._iterate_graph_topo()
        self._normalize_flux()
        self._is_valid()

    def _iterate_graph_topo(self):
        if not self.skeleton:
            raise RuntimeError('Skeleton does not exist to build from')
        accumulated = defaultdict(float)
        accumulated[self.bulk_node] = 1
        

        # Copy to plain DiGraph to avoid subgraph_view instantiation issues
        dag = nx.DiGraph(
            (u, v, d) for u, v, d in self.skeleton.edges(data=True)
            if v != self.bulk_node
        )

        for node in nx.topological_sort(dag):
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
            for path in nx.all_simple_paths(self, loading, self.bulk_node):
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
        for path in nx.all_simple_paths(self, u , v):
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

    def node_flux(self, node: RiboNode) -> float:
        total_flux = 0.0
        for _, _, flux in self.out_edges(nbunch=node, data='flux'):
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
    def simple(self) -> "SimpleFluxGraph":
        out = SimpleFluxGraph()
        out.bulk_node = self.bulk_node.simple
        out.flux_error = self.flux_error
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
        topo_nodes:list[RiboNode] = list(nx.topological_sort(out.dag))

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
        return out

    
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
    
    @property
    def weight_skeleton(self) -> RiboSkeleton:
        """Return a Riboskeleton based on the fluxes of this graph"""
        out = RiboSkeleton()
        for u,v in self.edges:
            out.add_transition(Transition(u,v,self.edge_weight(u,v)))
        return out
    

class SimpleFluxGraph(FluxGraph):
    def __init__(self, skeleton = None, flux_cutoff=0, **attr):
        super().__init__(skeleton, flux_cutoff, **attr)

    def edge_weight(self, u, v):
        return self[u][v]['flux_start'] / self.node_flux(u)
    
    def node_flux(self, node):
        total_flux = 0.0
        for _, _, flux in self.out_edges(nbunch=node, data='flux_start'):
            total_flux += flux
        return total_flux 