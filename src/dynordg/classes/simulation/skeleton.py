from ..core import RiboNode, Transition, State
from ..core.events import Pipe
from abc import ABC, abstractmethod
from dataclasses import dataclass
from collections import defaultdict
from ..graph import RiboGraph
from typing import Literal

@dataclass(frozen=True)
class PipeIndexEntry:
    pipe: Pipe
    kind: Literal["entry", "exit"]

    @property
    def position(self) -> int:
        if self.kind == "entry":
            return self.pipe.input_position
        return self.pipe.output_position
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
    def fraction(self, u:RiboNode, v:RiboNode) -> float:
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
    """Gain of ternary complex during scanning (phase 0). half_life=None returns a fraction of 1, instant reassocaition"""
    @property
    def priority(self) -> int:
        return 10

    def applies(self, u: RiboNode, v: RiboNode) -> bool:
        return u.phase == 0 and v.phase == 0 and "ternary_complex" not in u.state

    def fraction(self, u: RiboNode, v: RiboNode) -> float:
        if not self.half_life:
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
        return u.phase == v.phase and u.phase > 0 and "scanning_factors" in u.state

    def fraction(self, u: RiboNode, v: RiboNode) -> float:
        if self.half_life is None:
            return 0.0
        return 1 - 0.5 ** (abs(u.position - v.position) / self.half_life)

    def transitions(self, u: RiboNode, v: RiboNode, weight: float) -> Transition:
        lost = weight * self.fraction(u, v)
        return Transition(u, RiboNode(v.position, u.state - "scanning_factors"), lost)

class RiboSkeleton(RiboGraph):

    def __init__(
        self,
        pipelist: list[Pipe],
        incoming_graph_data=None,
        behaviours: list[FactorBehaviour] | None = None,
        **attr,
    ):
        super().__init__(incoming_graph_data, **attr)
        self.pipelist = pipelist
        self.pipelist.sort(reverse=True)
        
        self.stack: list[tuple[RiboNode, float]] = []
        self.behaviours: list[FactorBehaviour] | None = sorted(behaviours) if behaviours else None
        self.visited: set[RiboNode] = set()
        self._pipe_index = defaultdict(list)
        self.expanded: set[RiboNode] = set()
        self.queued: set[RiboNode] = set()
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
        for pipe in self.pipelist:
            if pipe.phase_condition.phases == -1:
                factors = pipe.subphase_conditions.required or None
                f_dict = {factor:True for factor in factors} if factors is not None else None
                source = RiboNode(pipe.input_position, State(-1, **f_dict)) if f_dict \
                    else RiboNode(pipe.input_position, State(-1))
                target = pipe.target(source.state)
                self.add_transition(Transition(self.bulk_node, source, 1))
                self.add_transition(Transition(source, target, pipe.probability))
                self.build_graph(target)
        self._is_valid()

    def build_graph(self, starting_node: RiboNode):
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

            node_at_pipe = RiboNode(pipe_indexes[0].position, node.state)
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