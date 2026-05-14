from .nodes import RiboNode, State
from abc import ABC, abstractmethod
from dataclasses import dataclass


        
@dataclass
class Transition:
    source: RiboNode
    target: RiboNode
    weight: float

class Pipe:
    output_phase: int
    probability: float
    position: int
    input_phase: int|tuple[int,...]|None = None
    remove_factors: str|tuple[str,...]|None = None
    add_factors: str|tuple[str,...]|None = None
    required_factors: str|tuple[str,...]|None = None
    exclude_factors: str|tuple[str,...]|None = None
    position_change: int|None = None

    def enter(self, state:State) -> bool:
        """Return true if state fulfils input conditions of this pipe"""
        if self.input_phase is not None:
            if isinstance(self.input_phase, int):
                phase = state.phase == self.input_phase
            elif isinstance(self.input_phase, tuple):
                phase = state.phase in self.input_phase
        else:
            phase = True

        if self.required_factors:
            factors = self.required_factors in state.factors
        else:
            factors = True
        if self.exclude_factors:
            if isinstance(self.exclude_factors, tuple):
                exclude = any(e in state.factors for e in self.exclude_factors)
            elif isinstance(self.exclude_factors, str):
                exclude = self.exclude_factors in state.factors
        else:
            exclude = False
        
        return phase and factors and not exclude
    
    def out_state(self, state:State) -> State:
        """Gives the finishing state of a ribosome with the inputted state passing through this pipe"""
        out = State(state.phase, *state.factors)
        if self.remove_factors:
            if self.remove_factors == 'all':
                out = State(out.phase)
            else:
                out -= self.remove_factors
        if self.add_factors:
            out += self.add_factors
        out.phase = self.output_phase
        return out 
    
    def target(self, state:State) -> RiboNode:
        """Returns the transition from a given node"""
        shift = self.position_change if self.position_change else 0
        return RiboNode(self.position + shift, self.out_state(state))
    
    def __repr__(self) -> str:
        if self.required_factors:
            return f"Pos:{self.position}, Phase:{self.input_phase}, Factors:{self.required_factors}"
        else:
            return f"Pos:{self.position}, Phase:{self.input_phase}"
        
    def __lt__(self, other):
        if not isinstance(other, Pipe):
            raise TypeError(f"< not supported between 'Pipe' and {type(other).__name__}'" )
        own_factors = self.required_factors if self.required_factors else []
        if isinstance(own_factors, str):
            own_factors = [1]
        other_factors = other.required_factors if other.required_factors else []
        if isinstance(other_factors, str):
            other_factors = [1]
        return len(own_factors) < len(other_factors)
    
class Event:
    FRAME_DICT = {1:1,
                  2:2,
                  0:3}
    
    def __init__(self, position: int, probability: float = 1):
        if not isinstance(position, int):
            raise TypeError(f"Event position must be type 'int', got '{type(position).__name__}'")
        if not isinstance(probability, (float, int)):
            raise TypeError(f"Event probability must be float or int, got '{type(probability).__name__}'")

        self.position = position
        self.probability = probability
        self.type = self.__class__.__name__

    @property
    @abstractmethod
    def pipe(self) -> Pipe:
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
    @property
    def pipe(self) -> Pipe:
        return Pipe(output_phase=self.frame,
                    probability=self.probability,
                    position=self.position,
                    input_phase=0,
                    remove_factors="ternary_complex",
                    required_factors="ternary_complex",
                    )
    
class Loading(Event):
    def __init__(self, position: int, probability: float):
        super().__init__(position, probability)
        self.source = RiboNode(self.position, State(-1))

class IRES(Loading):
    @property
    def pipe(self) -> Pipe:
        return Pipe(output_phase=self.frame,
                    probability=self.probability,
                    position=self.position,
                    input_phase=-1,
                    required_factors='ternary_complex',
                    remove_factors='ternary_complex')
    
class DropOff(Reading):
    def __init__(self, position: int, probability: float):
        super().__init__(position, probability)
        self.target = RiboNode(self.position, State(-1))

class Termination(DropOff):
    @property
    def pipe(self) -> Pipe:
        return Pipe(output_phase=-1,
                    probability=self.probability, 
                    position=self.position,
                    input_phase=self.frame,
                    )
    
class Retention(Reading):
    @property
    def pipe(self) -> Pipe:
        return Pipe(output_phase=0,
                    probability=self.probability,
                    position=self.position,
                    input_phase=self.frame,
                    required_factors='scanning_factors')
    
class Frameshift(Reading):
    def __init__(self, position, probability, amount: int):
        super().__init__(position, probability)
        self.amount = amount
        symbol = '+' if amount >= 0 else ''
        self.type = self.type + symbol + str(self.amount)

    @property
    def pipe(self) -> Pipe:
        frame_to = self.FRAME_DICT[(self.position + self.amount)%3]
        return Pipe(output_phase=frame_to,
                    probability=self.probability,
                    position=self.position,
                    input_phase=self.frame,
                    position_change=self.amount)
    
class LoadScanning(Loading):
    @property
    def pipe(self) -> Pipe:
        return Pipe(output_phase=0,
                    probability=self.probability,
                    position=self.position,
                    input_phase=-1,
                    required_factors=('ternary_complex', 'scanning_factors'),
                    add_factors=('ternary_complex', 'scanning_factors')
                    )
    
class AllDrop(DropOff):
    @property
    def pipe(self) -> Pipe:
        return Pipe(output_phase=-1,
                    probability=self.probability,
                    position=self.position)    
