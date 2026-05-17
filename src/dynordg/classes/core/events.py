from .nodes import RiboNode, State
from dataclasses import dataclass, field
from typing import Callable

        
@dataclass
class Transition:
    source: RiboNode
    target: RiboNode
    weight: float



@dataclass
class PhaseCondition:
    """Matches a set of phases. Pass a single int, a set, or a callable."""
    phases: int | set[int] | Callable[[int], bool]
    
    def match(self, phase: int) -> bool:
        if callable(self.phases):
            return self.phases(phase)
        if isinstance(self.phases, set):
            return phase in self.phases
        return phase == self.phases

# Pre-built helpers
ANY_TRANSLATING = PhaseCondition({1, 2, 3})
SCANNING       = PhaseCondition(0)
IN_SOLUTION    = PhaseCondition(-1)
SAME_FRAME     = None  # resolved at pipe-creation time using event.frame


@dataclass
class SubphaseCondition:
    """
    Class that holds conditions of subphase that must to return self.match=True
    minimum is inclusive
    maximum is exclusive
    """
    required:  tuple[str, ...]|str = field(default_factory=tuple)   # must ALL be present
    excluded:  tuple[str, ...]|str = field(default_factory=tuple)   # must NONE be present
    minimum:   dict[str, int|None]  = field(default_factory=dict)    # key >= value
    maximum:   dict[str, int|None]  = field(default_factory=dict)    # key <= value
    
    def __post_init__(self):
        if isinstance(self.required, str):
            self.required = (self.required,)
        if isinstance(self.excluded, str):
            self.excluded = (self.excluded,)

    def match(self, state: State) -> bool:
        factor_dict = dict(state.subphases)
        
        for f in self.required:
            if f not in factor_dict:
                return False
        for f in self.excluded:
            if f in factor_dict:
                return False
        for f, min_val in self.minimum.items():
            if min_val is not None:
                val = factor_dict.get(f)
                if val is None or (isinstance(val, int) and val < min_val):
                    return False
        for f, max_val in self.maximum.items():

            if max_val is not None:
                val = factor_dict.get(f)
                if max_val == 0:
                    return False
                if val is None:
                    continue

                if (isinstance(val, int) and val >= max_val):
                    return False
        return True
    
    def __lt__(self, other):
        if not isinstance(other, SubphaseCondition):
            raise TypeError(f"'<' not supported between instances of 'SubphaseCondition' and '{type(other).__name__}'")
        return len(self.required) < len(other.required)

    
@dataclass
class Pipe:
    output_phase:     int
    probability:      float
    input_position:   int
    output_position:  int|None = None
    phase_condition:  PhaseCondition   = field(default_factory=lambda: PhaseCondition({-1,0,1,2,3}))
    subphase_conditions: SubphaseCondition  = field(default_factory=SubphaseCondition)
    remove_factors:   str | tuple[str, ...] | None = None
    add_factors:      str | tuple[str, ...] | None = None

    def __post_init__(self):
        if self.output_position is None:
            self.output_position = self.input_position
    

    def enter(self, state: State) -> bool:
        return (
            self.phase_condition.match(state.phase)
            and self.subphase_conditions.match(state)
        )

    def out_state(self, state: State) -> State:
        out = State(self.output_phase, **dict(state.subphases))
        
        if self.remove_factors == 'all':
            out = State(self.output_phase)
        elif self.remove_factors:
            factors = (self.remove_factors,) if isinstance(self.remove_factors, str) else self.remove_factors
            for f in factors:
                out = out - f
                
        if self.add_factors:
            factors = (self.add_factors,) if isinstance(self.add_factors, str) else self.add_factors
            for f in factors:
                out = out + f
        return out

    def target(self, state: State) -> RiboNode:
        if not isinstance(self.output_position, int):
            raise RuntimeError(f'Output position was not set for this pipe: {self}')
        return RiboNode(self.output_position, self.out_state(state))
    
    def __repr__(self) -> str:
        inpos = f"Pos:{self.input_position}"
        outpos = f"Pos:{self.output_position}"
        inphase = f", Phase:{self.phase_condition}" if self.phase_condition else ""
        outphase = f", Phase:{self.output_phase}"
        infactors = f", Factors:{self.subphase_conditions}" if self.subphase_conditions else ""
        add = f", Add:{self.add_factors}" if self.add_factors else ""
        remove = f", Remove{self.remove_factors}" if self.remove_factors else ""
        return f"\n>PIPE>\n" \
               f"IN | {inpos}{inphase}{infactors}\n" \
               f"OUT| {outpos}{outphase}{add}{remove}"

    def __lt__(self, other):
        if not isinstance(other, Pipe):
            raise TypeError(f"'<' not supported between instances of 'Pipe' and '{type(other).__name__}'")
        return self.subphase_conditions < other.subphase_conditions
    
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
    def pipe(self) -> Pipe:
        """Return transitions for this event"""
        return Pipe(output_phase=-1,
                    probability=self.probability,
                    input_position=self.position)


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
        pc = PhaseCondition(0)
        spc = SubphaseCondition(required="ternary_complex")
        return Pipe(self.frame,
                    self.probability,
                    self.position,
                    phase_condition=pc,
                    subphase_conditions=spc,
                    remove_factors='ternary_complex')
    

class IRES(Event):
    @property
    def pipe(self) -> Pipe:
        return Pipe(output_phase=self.frame,
                    probability=self.probability,
                    input_position=self.position,
                    phase_condition=PhaseCondition(-1),
                    subphase_conditions=SubphaseCondition('ternary_complex'),
                    remove_factors='ternary_complex')
    

class Termination(Reading):
    @property
    def pipe(self) -> Pipe:
        return Pipe(output_phase=-1,
                    probability=self.probability,
                    input_position=self.position,
                    phase_condition=PhaseCondition(self.frame),
                    remove_factors='all'
                    )
    
class Retention(Reading):
    def __init__(self, position: int, probability: float, limit: int|None = None):
        """Limit is the max number of times a ribosome may have previously initiated"""
        super().__init__(position, probability)
        self.limit = limit
        

    @property
    def pipe(self) -> Pipe:
        return Pipe(output_phase=0,
                    probability=self.probability,
                    input_position=self.position,
                    phase_condition=PhaseCondition(self.frame),
                    subphase_conditions=SubphaseCondition('scanning_factors', maximum={'retained':self.limit}),
                    add_factors='retained'
                    )
    
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
                    input_position=self.position,
                    output_position=self.position+self.amount,
                    phase_condition=PhaseCondition(self.frame),
                    )
    
class LoadScanning(Event):
    @property
    def pipe(self) -> Pipe:
        return Pipe(output_phase=0,
                    probability=self.probability,
                    input_position=self.position,
                    phase_condition=PhaseCondition(-1),
                    subphase_conditions=SubphaseCondition(tuple(('ternary_complex', 'scanning_factors')))
                    )
    
class AllDrop(Reading):
    @property
    def pipe(self) -> Pipe:
        return Pipe(output_phase=-1,
                    probability=self.probability,
                    input_position=self.position)    
