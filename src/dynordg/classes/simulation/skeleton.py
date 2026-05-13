from ..graph import RiboGraph
from .transitionmap import TransitionMap
from ..core import RiboNode, Transition, State
from ..core.events import Pipe
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
        pipelist: list[Pipe],
        incoming_graph_data=None,
        behaviours: list[FactorBehaviour] | None = None,
        **attr,
    ):
        super().__init__(incoming_graph_data, **attr)
        self.pipelist = pipelist
        print('pipelist:', self.pipelist)
        self.stack: list[tuple[RiboNode, float]] = []
        self.behaviours: list[FactorBehaviour] | None = sorted(behaviours) if behaviours else None
        self.visited: list[RiboNode] = []
        self._construct()

    def _construct(self):
        for pipe in self.pipelist:
            if isinstance(pipe.input_phase, int) and pipe.input_phase == -1:
                factors = pipe.required_factors or []
                source = RiboNode(pipe.position, State(pipe.input_phase, *factors))
                target = pipe.target(source.state)
                self.add_transition(Transition(self.bulk_node, source, 1))
                self.add_transition(Transition(source, target, pipe.probability))
                self.build_graph(target)
        self._is_valid()

    def build_graph(self, starting_node: RiboNode):
        self.stack.append((starting_node, 1))
        while self.stack:
            node, weight = self.stack.pop()
            self.visited.append(node)

            pipes = self.next_pipes(node)
            if not pipes:
                raise RuntimeError(
                    f"Continue flux at node: {node} runs to infinity (no downstream pipes)"
                )

            node_at_pipe = RiboNode(pipes[0].position, node.state)
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
        """Route a node either back to bulk (phase == -1) or through downstream pipes."""
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
        viable = [p for p in pipes if p.enter(node.state)]
        for pipe in reversed(sorted(viable)):
            if unpiped < 0:
                break
            pipe_target = pipe.target(node.state)
            self.add_transition(
                Transition(node, pipe_target, pipe.probability * unpiped)
            )
            unpiped -= pipe.probability * unpiped
            if pipe_target.phase == -1:
                self.add_transition(Transition(pipe_target, self.bulk_node, 1))
            else:
                self._add_to_stack(pipe_target, 1)
        return unpiped

    def _add_to_stack(self, node: RiboNode, weight: float):
        if node not in self.visited:
            self.stack.append((node, weight))

    def next_pipes(self, node: RiboNode) -> list[Pipe] | None:
        def pipe_matches(pipe: Pipe) -> bool:
            return (
                pipe.input_phase == node.phase
                or pipe.input_phase is None
                or pipe.output_phase == node.phase
            )

        downstream_positions = [
            pipe.position for pipe in self.pipelist
            if pipe_matches(pipe) and pipe.position > node.position
        ]
        if not downstream_positions:
            return None

        min_pos = min(downstream_positions)
        return [
            pipe for pipe in self.pipelist
            if pipe_matches(pipe) and pipe.position == min_pos
        ]
    

# class RiboSkeleton(TransitionMap):

#     def __init__(
#         self,
#         pipelist: list[Pipe],
#         incoming_graph_data=None,
#         behaviours: list[FactorBehaviour] | None = None,
#         **attr,
#     ):
#         super().__init__(incoming_graph_data, **attr)
#         self.pipelist = pipelist


#         self.stack: list[tuple[RiboNode, float]] = []
#         self.behaviours: list[FactorBehaviour]|None = sorted(behaviours) if behaviours else None
#         self.visited: list[RiboNode] = []
#         self._construct()

#     def _construct(self):
#         for pipe in self.pipelist:
#             if isinstance(pipe.input_phase, int) and pipe.input_phase == -1:
#               factors = pipe.required_factors if pipe.required_factors else []
#               source = RiboNode(pipe.position, State(pipe.input_phase, *factors))
#               target = pipe.target(source.state)
#               from_bulk = Transition(self.bulk_node, source, 1)
#               self.add_transition(from_bulk)
#               self.add_transition(Transition(source, target, pipe.probability))
#               self.build_graph(target)
#         self._is_valid()

#     def build_graph(self, starting_node: RiboNode):
#         self.stack.append((starting_node, 1))
#         while self.stack:
#             node, weight = self.stack.pop()
#             self.visited.append(node)
#             pipes = self.next_pipes(node)
#             if not pipes:
#                 raise RuntimeError(f"Continue flux at node: {node} runs to infinity (no downstream pipes)")
#             remaining_weight = weight

#             node_at_pipe = RiboNode(pipes[0].position, node.state)
#             if self.behaviours:
#                 for behaviour in self.behaviours:
#                     if behaviour.applies(node, node_at_pipe):
#                         f_trans = behaviour.transitions(node, node_at_pipe, remaining_weight)
#                         if f_trans.weight > 0:
#                             remaining_weight -= f_trans.weight
#                             self.add_transition(f_trans)
#                             if f_trans.target.phase == -1:
#                                 self.add_transition(Transition(f_trans.target, self.bulk_node, 1))
#                             else:
#                                 unpiped_weight = 1
#                                 for pipe in pipes:
#                                     if pipe.enter(f_trans.target.state):
#                                         pipe_target = pipe.target(f_trans.target.state)
#                                         self.add_transition(Transition(f_trans.target,
#                                                                     pipe_target,
#                                                                     pipe.probability))
#                                         unpiped_weight -= pipe.probability
#                                         if pipe_target.phase == -1:
#                                             self.add_transition(Transition(pipe_target, self.bulk_node, 1))
#                                         else:  
#                                             self._add_to_stack(pipe_target, 1)
#                             if unpiped_weight > 0:
#                                 self._add_to_stack(f_trans.target, unpiped_weight)
#             if remaining_weight > 0:
#                 self.add_transition(Transition(node, node_at_pipe, remaining_weight))
#                 unpiped_weight = 1
#                 viable_pipes =  [pipe for pipe in pipes if pipe.enter(node_at_pipe.state)]
#                 if len(viable_pipes):
#                     for pipe in reversed((viable_pipes)):
#                         pipe_target = pipe.target(node_at_pipe.state)
#                         pipe_transition = Transition(node_at_pipe,
#                                                     pipe_target,
#                                                     pipe.probability*unpiped_weight)
#                         self.add_transition(pipe_transition)
#                         unpiped_weight -= pipe.probability*unpiped_weight
#                         if pipe_target.phase == -1:
#                             self.add_transition(Transition(pipe_target, self.bulk_node, 1))
#                         else:  
#                             self._add_to_stack(pipe_target, 1)
#                 if unpiped_weight > 0:
#                     self._add_to_stack(node_at_pipe, unpiped_weight)



#     def _add_to_stack(self, node:RiboNode, weight: float):
#         if node not in self.visited:
#             self.stack.append((node, weight))

#     def next_pipes(self, node: RiboNode) -> list[Pipe]|None:
#         same_phase = [pipe.position for pipe in self.pipelist
#                       if (pipe.input_phase == node.phase or pipe.input_phase == None or pipe.output_phase == node.phase)
#                       and pipe.position > node.position]
#         if not len(same_phase):
#             return None
#         min_pos = min(same_phase)

#         return [pipe for pipe in self.pipelist 
#                 if (pipe.input_phase == node.phase 
#                     or pipe.input_phase==None 
#                     or pipe.output_phase==node.phase) 
#                     and pipe.position == min_pos]
