from .nodes import RiboNode, State
from .transitions import RiboTransition
import re
from abc import ABC, abstractmethod
from dataclasses import dataclass


        
@dataclass
class Transition:
    source: RiboNode
    target: RiboNode
    weight: float

class Event:
    FRAME_DICT = {1:1,
                  2:2,
                  0:3}
    
    def __init__(self, position: int, probability: float):
        if not isinstance(position, int):
            raise ValueError(f"Event position must be type 'int', got '{type(position).__name__}'")
        if not isinstance(probability, (float, int)):
            raise ValueError(f"Event probability must be float or int, got '{type(probability).__name__}'")

        self.position = position
        self.probability = probability

    @abstractmethod
    def transitions(self) -> list[Transition]:
        """Return transitions for this event"""
        ...

    @property
    def frame(self):
        return self.FRAME_DICT[self.position%3]

    def __repr__(self):
        return f"{self.__class__.__name__}(pos={self.position}, prob={self.probability})"

class Reading(Event):
    def __init__(self, position: int, probability: float):
        if probability > 1 or probability <= 0:
            raise ValueError(f"Event probability must be greater than 0 and not more than 1, got {probability}")
        super().__init__(position, probability)


class Initiation(Reading):
    def transitions(self) -> list[Transition]:
        source = RiboNode(self.position, State(0, 'ternary_complex'))
        target = RiboNode(self.position, State(self.frame, 'scanning_factors'))
        return [Transition(source, target, self.probability)]
    
class Loading(Event):
    def __init__(self, position: int, probability: float):
        super().__init__(position, probability)
        self.source = RiboNode(self.position, State(-1))

class IRES(Loading):
    def transitions(self) -> list[Transition]:
        target = RiboNode(self.position, State(self.frame))
        return [Transition(self.source, target, self.probability)]
    
class DropOff(Reading):
    def __init__(self, position: int, probability: float):
        super().__init__(position, probability)
        self.target = RiboNode(self.position, State(-1))

class Termination(DropOff):
    def transitions(self) -> list[Transition]:
        source = 
        
    
class RiboEvent(tuple):
    """
    Tuple in the form position:int, type: str, probability: float, drop_probability: float
    Probabilities must be 0 <= p <= 1 and must not add to be greater than 1 (except termination which \
    may have a probability (of reinitiation) and a drop_probability).
    """


    def _initiation_transition(self):
        transitions = []
        if self.probability > 0:
            transitions.append(RiboTransition(RiboNode(self.position, 0, True), 
                                              RiboNode(self.position, self.frame, True), 
                                              self.probability))
        return transitions
    
    def _termination_transition(self):
        transitions = []
        if self.probability > 0:
            transitions.append(RiboTransition(RiboNode(self.position, self.frame, False),
                                            RiboNode(self.position, -1, False),
                                            self.probability))

        return transitions
    
    def _frameshift_transition(self):
        transitions = []
        shift_to_pos = self.position + self.shift
        shift_to_frame = shift_to_pos % 3 + 1
        if self.probability > 0:
            for sub in (False, True):
                transitions.append(RiboTransition(RiboNode(self.position, self.frame, sub),
                                                RiboNode(shift_to_pos, shift_to_frame, sub ), 
                                                self.probability))
        return transitions
    
    def _ires_transition(self):
        transitions = []
        if self.probability > 0:
            transitions.append(RiboTransition(RiboNode(self.position, -1, True), 
                                            RiboNode(self.position, self.frame, True), 
                                            self.probability))
        return transitions

    def _load_transition(self):
        transitions = []
        if self.probability > 0:
            transitions.append(RiboTransition(RiboNode(self.position, -1, True), 
                                            RiboNode(self.position, 0, True), 
                                            self.probability))
        return transitions
    

    def _end_transition(self):
        transitions = []
        if self.probability > 0:
            for phase in range(4):
                for subphase in (False, True):
                    transitions.append(RiboTransition(RiboNode(self.position, phase, subphase),
                                                      RiboNode(self.position, -1, False),
                                                      1))
        return transitions
    
    def _retention_transition(self):
        transitions = []
        if self.probability > 0:
            transitions.append(RiboTransition(RiboNode(self.position, self.frame, True), 
                                              RiboNode(self.position, 0, False), 
                                              self.probability))
        return transitions

    def _scanning_drop(self):
        transitions = []
        if self.probability > 0:
            for factor in (True, False):
                transitions.append(RiboTransition(RiboNode(self.position, 0, factor),
                                                RiboNode(self.position, -1, factor),
                                                self.probability))
        return transitions

    TRANSITION_MAP = {
        'initiation':  _initiation_transition,
        'termination':    _termination_transition,
        'frameshift+1': _frameshift_transition,
        'frameshift-1': _frameshift_transition,
        'ires':  _ires_transition,
        'cap': _load_transition,
        'load_scanning': _load_transition,
        'all_drop': _end_transition,
        'end': _end_transition,
        'retention': _retention_transition
    }




    @property
    def frame(self):
        """
        Provides the frame (1,2,3) of the event
        """
        return (self.position % 3) + 1
    
    

    @property
    def shift(self):
        """
        If event is a frameshift type, this gives the direction and magnitude of the frameshift
        """
        match = re.search(r'shift([+-]?\d+)', self.type)
        return int(match.group(1)) if match else 0

    def __new__(cls, *args):
            
            if len(args) == 1 and isinstance(args[0], (tuple, RiboEvent)):
                data = args[0]
            elif len(args) == 3:
                data = args
            else:
                raise ValueError(f'RiboEvent requires either or a length-3 tuple, got: {args}')

            if len(data) != 3:
                raise ValueError(f'RiboEvent tuple must be of length 3, got length: {len(data)}')
            

            if not isinstance(data[0], int):
                raise ValueError(f"Position of RiboEvent tuple must be 'int', got {type(data[0]).__name__!r}")
            elif not isinstance(data[1], str):
                raise ValueError(f"Type of RiboEvent tuple must be 'str', got {type(data[1]).__name__!r}")
            elif not isinstance(data[2], (float, int)):
                raise ValueError("Probabilities must be float or int")
            elif data[2] <= 0:
                    raise ValueError(f'Probability must be greater than 0, got {data[2]}')

            elif data[2] >1 and data[1] not in ('load_scanning', 'cap', 'ires'):
                raise ValueError (f'Non-loading probabilities cannot exceed 1, got {data[2]}')
            return super().__new__(cls, data)
    
    def __init__(self, *args):
        super().__init__()
        self.position = self[0]
        self.type = self[1]
        self.probability = self[2]

    def __repr__(self):
        return f"(Pos:{self.position}, type:{self.type}, prob:{self.probability})"
    
    def to_transition(self) -> list[RiboTransition]:
        """
        Returns a list of RiboTransitions coresponding to the event
        """
        handler = self.TRANSITION_MAP.get(self.type)
        if handler is None:
            raise ValueError(f"No transition defined for event type {self.type!r}")
        return handler(self)