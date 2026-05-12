from ..graph import RiboGraph
from .transitionmap import TransitionMap
from ..core import RiboNode, Transition, State
import networkx as nx
import matplotlib.pyplot as plt
from abc import ABC, abstractmethod
from dataclasses import dataclass

# ------------------------------------------------------------------ #
#  Base class for factor behaviours                                    #
# ------------------------------------------------------------------ #
from functools import total_ordering

@total_ordering
class FactorBehaviour(ABC):
    """
    Defines what happens to a ribosome (or its associated factors)
    as it moves between two nodes. Subclass this to define custom behaviour.
    """
    def __init__(self, half_life:int|None=None) -> None:
        super().__init__()
        self.half_life = half_life

    @abstractmethod
    def applies(self, u: RiboNode, v: RiboNode) -> bool:
        """Return True if this behaviour should fire for the u -> v transition."""
        ...

    @abstractmethod
    def fraction(self) -> float:
        """
        Returns the flux fraction that undergoes this state change
        """
        ...

    @abstractmethod
    def transitions(self, u: RiboNode, v: RiboNode, weight: float) -> Transition:
        """
        Return the transitions produced by this behaviour.
        """
        ...

    @property
    @abstractmethod
    def priority(self) -> int:
        """
        Lower value = applied first.
        Define a consistent priority scheme, e.g.:
          0  : structural / decay (always first)
          10 : factor acquisition
          20 : factor loss
        """
        ...

    def __eq__(self, other: object) -> bool:
        return isinstance(other, FactorBehaviour) and self.priority == other.priority

    def __lt__(self, other: "FactorBehaviour") -> bool:
        return self.priority < other.priority


# ------------------------------------------------------------------ #
#  Built-in behaviours                                                 #
# ------------------------------------------------------------------ #
class PhaseDecay(FactorBehaviour):
    
    @property
    def priority(self) -> int:
        return 0
    
    def fraction(self, u: RiboNode, v: RiboNode) -> float:
        if self.half_life:
            return 1 - 0.5 ** (abs(u.position - v.position) / self.half_life)
        return 0
    
    def transitions(self, u: RiboNode, v: RiboNode, weight: float) -> Transition:
        decay = weight * self.fraction(u, v)
        return Transition(u, RiboNode(v.position, State(-1)), decay)

class ScanningDecay(PhaseDecay):
    """Ribosome decay during scanning"""
    def applies(self, u: RiboNode, v: RiboNode) -> bool:
        return u.phase == v.phase and u.phase == 0

class TranslationDecay(PhaseDecay):
    """Ribosome decay during translation"""
    def applies(self, u: RiboNode, v: RiboNode) -> bool:
        return u.phase == v.phase and u.phase > 0
    
class TernaryComplexAssociation(FactorBehaviour):
    """Loss of ternary complex during scanning (phase 0)."""
    @property
    def priority(self) -> int:
        return 10

    def applies(self, u: RiboNode, v: RiboNode) -> bool:
        return u.phase == 0 and v.phase == 0 and "ternary_complex" not in u.factors

    def fraction(self, u: RiboNode, v: RiboNode) -> float:
        if self.half_life is None:
            return 1.0
        return 1 - 0.5 ** (abs(u.position - v.position) / self.half_life)

    def transitions(self, u: RiboNode, v: RiboNode, weight: float) -> Transition:
        prop = weight * self.fraction(u, v)
        return Transition(u, RiboNode(v.position, u.state + "ternary_complex"), prop)


class ScanningFactorDissociation(FactorBehaviour):
    """Loss of scanning factors during translation (phase > 0)."""

    @property
    def priority(self) -> int:
        return 20

    def applies(self, u: RiboNode, v: RiboNode) -> bool:
        return u.phase == v.phase and u.phase > 0 and "scanning_factors" in u.factors

    def fraction(self, u: RiboNode, v: RiboNode) -> float:
        if self.half_life is None:
            return 0.0
        return 1 - 0.5 ** (abs(u.position - v.position) / self.half_life)

    def transitions(self, u: RiboNode, v: RiboNode, weight: float) -> Transition:
        lost = weight * self.fraction(u, v)
        return Transition(u, RiboNode(v.position, u.state - "scanning_factors"), lost)


class RiboSkeleton(TransitionMap):

    def __init__(
        self,
        transition_map: TransitionMap,
        incoming_graph_data=None,
        behaviours: list[FactorBehaviour] | None = None,
        **attr,
    ):
        super().__init__(incoming_graph_data, **attr)
        self.transitions = transition_map
        self.cont_weight = self._continuation_weights()
        self.next_node = {
            n: self.transitions.downstream_node(n)
            for n in self.transitions.nodes
            if n.phase != -1
        }

        self.behaviours: list[FactorBehaviour]|None = sorted(behaviours) if behaviours else None

        self._construct()

    def _construct(self):

        for u, v, w in self.transitions.edges(data="weight"):
            self.add_transition(Transition(u, v, w))
            self._attach_node(u, is_source=True)
            self._attach_node(v, is_source=False)

        self.remove_all_indegree_zero_nodes()

        


    def _continuation_weights(self) -> dict[RiboNode,float]:
        """Weight remaining after all explicit out-edges (i.e. the 'continue' share)."""
        weights = {}
        for node in self.transitions.nodes:
            out_total = sum(w for _, _, w in self.transitions.out_edges(node, data="weight"))
            weights[node] = 1 - out_total
        return weights
    
    def _attach_node(self, node: RiboNode, is_source: bool):
        """Wire up same_phase transitions for a single node."""
        if node.phase == -1:
            if is_source:
                self.add_transition(Transition(self.bulk_node, node, 1))
            else:
                self.add_transition(Transition(node, self.bulk_node, 1))
            return

        tnode = self.next_node.get(node)

        if tnode is None:
            drop = RiboNode(node.position, State(-1))
            self.add_transition(Transition(node, drop, self.cont_weight[node]))
            self.add_transition(Transition(drop, self.bulk_node, 1))
            return

        remaining_weight = self.cont_weight[node]

        if self.behaviours:
            for behaviour in self.behaviours:
                if behaviour.applies(node, tnode):
                    t = behaviour.transitions(node, tnode, remaining_weight)
                    if t.weight > 0:
                        remaining_weight -= t.weight
                        self.add_transition(t)
                        if t.target.phase == -1:
                            self.add_transition(Transition(t.target, self.bulk_node, 1))
                        elif t.target not in self.transitions.nodes:
                            # New node introduced by this behaviour — extend the graph from it
                            self.next_node[t.target] = self.transitions.downstream_node(t.source) 
                            self.cont_weight[t.target] = 1
                            self._attach_node(t.target, is_source=True)

        if remaining_weight:
            self.add_transition(Transition(node, tnode, remaining_weight))