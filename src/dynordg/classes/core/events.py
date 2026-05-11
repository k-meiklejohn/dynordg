from .nodes import RiboNode, State
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
            raise TypeError(f"Event position must be type 'int', got '{type(position).__name__}'")
        if not isinstance(probability, (float, int)):
            raise TypeError(f"Event probability must be float or int, got '{type(probability).__name__}'")

        self.position = position
        self.probability = probability
        self.type = self.__class__.__name__

    @abstractmethod
    def transitions(self) -> list[Transition]:
        """Return transitions for this event"""
        ...
    @property
    def __members(self):
        return (self.position, self.probability, self.type)
    
    def __eq__(self, value):
        if  isinstance(value, Event):
            return self.__members == value.__members

    def __hash__(self):
        return hash(self.__members)

    def __lt__(self, other):
        if not isinstance(other, Event):
            raise TypeError(f"< not supported between instances of '{self.type}' and '{type(other).__name__}")
        if self.position > other.position:
            return False
        elif self.position < other.position:
            return True
        elif self.type < other.type:
            return True
        elif self.probability < other.probability:
            return True
        else:
            return False

    @property
    def frame(self):
        return self.FRAME_DICT[self.position%3]

    def __repr__(self):
        return f"{self.type}(pos={self.position}, prob={self.probability})"

class Reading(Event):
    def __init__(self, position: int, probability: float):
        if probability > 1 or probability <= 0:
            raise ValueError(f"Event probability must be greater than 0 and not more than 1, got {probability}")
        super().__init__(position, probability)

class Initiation(Reading):
    def transitions(self) -> list[Transition]:
        source = RiboNode(self.position, State(0, 'ternary_complex', 'scanning_factors'))
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
        source = RiboNode(self.position, State(self.frame))
        return [Transition(source, self.target, self.probability)]
    
class Retention(Reading):
    def transitions(self):
        source = RiboNode(self.position, State(self.frame, 'scanning_factors'))
        target = RiboNode(self.position, State(0, 'scanning_factors'))
        return [Transition(source, target, self.probability )]
    
class Frameshift(Reading):
    def __init__(self, position, probability, amount: int):
        super().__init__(position, probability)
        self.amount = amount
        symbol = '+' if amount >= 0 else ''
        self.type = self.type + symbol + str(self.amount)

    def transitions(self):
        frame_to = self.FRAME_DICT[(self.position + self.amount)%3]
        source = RiboNode(self.position, State(self.frame))
        target = RiboNode(self.position + self.amount, State(frame_to))
        return
    
class LoadScanning(Loading):
    def transitions(self):
        target = RiboNode(self.position, State(0, 'scanning_factors', 'ternary_complex'))
        return [Transition(self.source, target, self.probability)]
    
class AllDrop(DropOff):
    def transitions(self):
        out:list[Transition] = []
        for phase in range(4):
            source = RiboNode(self.position, State(phase))
            out.append(Transition(source, self.target, self.probability))
        return out
    
