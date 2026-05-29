"""
Core objects for describing Ribosomal phase space and the movemnt of ribsomes through it

Classes:

State

"""

from dataclasses import dataclass, field
from typing import Callable
from abc import ABC, abstractmethod
from functools import total_ordering
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from .scoring import noderer_start_score
import Levenshtein

class State:
    """
    Class object representing a State in ribosomal phase space

    
    Args:
    phase -- an integer between -1 and 3 inclusive 
        * -1 refers to a ribosome being in solution 
        * 0 refers to a scanning ribosome
        * 1,2,3 refer to translation in each frame when the first subcodon position that the ribosome reads
        is in frame with the first, second and thrid nucleotide of the transcript respecitvely

    subphases -- keywords can be supplied to State to specify the subphase of the State. 
        Used by default by dynordg are ternary_complex and scanning_factors when using the defaul Event and FactorBehaviour types

    Attributes:
    phase -- returns the phase of the state
    subphase -- returns the subphase dictionary of the state

    Methods:
    '+'  operator -- takes a string and adds a new subphase specification
                     or increments the value of the specified key by 1 if
                     the value is an integer.
                     returns a new State object
    '-'  operator -- takes a string and removes the subphase from the State 
                     or decremnts the corresponding value by 1 if it is an integer
                     returns a new State object

    '<' operator -- sorts by phase and then subphase
    
    'in' operator -- returns true if the specified string is in the subphase keys
    """

    ___slots__ = ("phase", "_subphases", "_keys", "_hash")  
    def __init__(self, phase: int, **subphases):

        if not isinstance(phase, int):
            raise ValueError(f'RiboNode phase must be int, got {phase}')
        if phase < -1 or phase > 3:
            raise ValueError(f'Phase must be int, between -1 and 3 inclusive, got {phase}')
        self.phase = phase
        canonical = tuple(sorted(subphases.items()))
        self._subphases = canonical
        self._hash = hash((phase, canonical))
        self._keys = frozenset(subphases)
        self._map = dict(subphases)

    @property
    def subphases(self):
        return self._subphases

    def __repr__(self):
        if len(self.subphases):
            return f"Phase:{self.phase} Factors:{self.subphases}"
        else:
            return f"Phase:{self.phase}"
        
    def __eq__(self, value):
        if isinstance(value, State):
            return self._hash == value._hash and self._subphases == value._subphases and self.phase == value.phase
        return False

    def __hash__(self) -> int:
        return self._hash
        

    def __add__(self, other):
        other = str(other)
        out_dict = dict(self.subphases)
        if other in self._map:
            if isinstance(self._map[other], int):
                out_dict[other] += 1
                return State(self.phase, **out_dict)
            else:
                return
    
        out_dict[other] = 1
        return State(self.phase, **out_dict)

    def __sub__(self, other):
        other = str(other)
        out_dict = dict(self.subphases)
        if other in self._map:
            if isinstance(self._map[other], int):
                out_dict[other] -= 1
                if out_dict[other] == 0:
                    out_dict.pop(other)
                return State(self.phase, **out_dict)
        
            out_dict.pop(other)
            return State(self.phase, **out_dict)
        return self
    
    def __contains__(self, item):
        return item in self._keys
    
    def __lt__(self, other):
        if isinstance(other, State):
            if self == other:
                return False
            elif self.phase < other.phase:
                return True
            elif self.subphases < other.subphases:
                return True 
            else:
                return False
        else:
            raise TypeError(f"Operand < not supported between 'State' and '{type(other).__name__}'")
            
    def __gt__(self, other):
        if isinstance(other, State):
            if self == other:
                return False
            elif self.phase > other.phase:
                return True
            elif self.subphases > other.subphases :
                return True 
            else:
                return False
        else:
            raise TypeError(f"Operand > not supported between 'State' and '{type(other).__name__}'")

class RiboNode:
    """
    Class object representing a position in Ribosomal Phase Space.

    Args:

    position -- 1-indexed nucleotide position in the transcript
    state    -- State object representing the state of the ribosome

    Attributes:

    position  -- returns the transcript position of the node
    phase     -- returns the phase of the state of the node
    subphases -- returns the subphases of the state of the node
    simple    -- returns a position, phase only version of the node

    For sorting; position takes precedent over state

    """
    def __init__(self, position: int, state: State):
        if not isinstance(position, int):
            raise TypeError(f'RiboNode position must be int, got {position}')
        self.position: int = position
        self.state: State = state
        self._members = (position, state)      # computed once
        self._hash = hash(self._members)       # computed once
    
        
    @property
    def phase(self) -> int:
        """Return the phase of the node"""

        return self.state.phase
    
    @property
    def subphases(self) -> tuple[tuple,...]:
        """Return the subphase of the node"""
        return self.state.subphases

    @property
    def simple(self):
        """Return a new node object with only position and phase"""
        return RiboNode(self.position, State(self.phase))

    def __repr__(self):
        if len(self.state.subphases):
            return f"(Pos:{self.position}, Phase:{self.phase}, F:{self.subphases})"
        else:
            return f"(Pos:{self.position}, Phase:{self.phase})"
        
    def __eq__(self, value):
        if isinstance(value, RiboNode):
            return (self._hash == value._hash and 
                    self.position == value.position and 
                    self.state._hash == value.state._hash and
                    self.state._subphases == value.state._subphases)
        return False

    
    def __hash__(self) -> int:
        return self._hash
    
    def __add__(self, other):
        return RiboNode(self.position, self.state + other)
    
    def __lt__(self, other):
        if isinstance(other, RiboNode):
            if self == other:
                return False
            elif self.position < other.position:
                return True
            elif self.state < other.state:
                return True
            else:
                return False
            
        else:
            raise TypeError(f"Operand < not supported between 'RiboNode' and '{type(other).__name__}'")
            
    def __gt__(self, other):
        if isinstance(other, RiboNode):
            if self == other:
                return False
            elif self.position > other.position:
                return True
            elif self.state > other.state:
                return True
            else:
                return False
            
        else:
            raise TypeError(f"Operand > not supported between 'RiboNode' and '{type(other).__name__}'")



@dataclass
class Transition:
    """Defines a weighted edge in Ribosomal Phase Space"""
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
    Class that holds conditions of subphase that must be true to return self.match=True

    Parameters:
    required -- a tuple of keys that must be present
    excluded -- a tuple of keys that cannot be present
    minimum  -- a dictionary of keys and integers, where the value
                of the test must be at least the value of the condition
    ceiling  -- a dictionary of keys and integers, where the value
                of the test cannot be greater than or equal to
                the value of the condition

    Methods
    match  -- returns true if the provided state matches the the conditions
    """
    required:  tuple[str, ...]|str = field(default_factory=tuple)   # must ALL be present
    excluded:  tuple[str, ...]|str = field(default_factory=tuple)   # must NONE be present
    minimum:   dict[str, int|None]  = field(default_factory=dict)    # key >= value
    ceiling:   dict[str, int|None]  = field(default_factory=dict)    # key <= value
    
    def __post_init__(self):
        if isinstance(self.required, str):
            self.required = (self.required,)
        if isinstance(self.excluded, str):
            self.excluded = (self.excluded,)

    def match(self, state: State) -> bool:
        """Returns True if the provided state matches the condtions"""
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
        for f, max_val in self.ceiling.items():

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
    """
    Class representing the effect of an Event on a ribosome

    A pipe takes in a ribosome at a certain position or a range of positions in Ribosomal Phase Space
    and outputs it to a different position in RPS

    Parameters:
    output_phase        -- an integer representing the phase of the pipes exit
    probability         -- the proportion of ribosomes that reach the pipe that enter it (0<x<=1)
    input_position      -- the nucleotide position of the pipes entry point
    output_position     -- the nucleotide position of the pipes exit point
    phase_condition     -- A PhaseCondition object that representing the phase(s) of the  entry point
    subphase_conditions -- A SubphaseCondition object representing the subphase(s) of the entry point
    remove_factors      -- The factors that are consumed when moving through the pipe
    add_factors         -- The factors that are added when moving through the pipe

    Methods:

    enter - checks if a state matches the pipe entry point
        params:
            state -- The State object to be checked
        returns:
            bool

    out_state - returns the state after piping
        params:
            state -- state before entering the pipe
        returns:
            State

    target - returns the target of the pipe given a state
        params:
            state -- state entering the pipe
        returns:
            RiboNode
    """
    output_phase:     int
    probability:      float
    input_position:   int
    output_position:  int|None = None
    phase_condition:  PhaseCondition   = field(default_factory=lambda: PhaseCondition({-1,0,1,2,3}))
    subphase_condition: SubphaseCondition  = field(default_factory=SubphaseCondition)
    remove_factors:   str | tuple[str, ...] | None = None
    add_factors:      str | tuple[str, ...] | None = None

    def __post_init__(self):
        if self.output_position is None:
            self.output_position = self.input_position
    

    def enter(self, state: State) -> bool:
        """Return True if the provided state matches the pipe entry"""
        return (
            self.phase_condition.match(state.phase)
            and self.subphase_condition.match(state)
        )

    def out_state(self, state: State) -> State:
        """Returns the state of the ribosome after moving through the pipe"""
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
        """Returns the final RiboNode of a ribosome entering the pipe """
        if not isinstance(self.output_position, int):
            raise RuntimeError(f'Output position was not set for this pipe: {self}')
        return RiboNode(self.output_position, self.out_state(state))
    
    def __repr__(self) -> str:
        inpos = f"Pos:{self.input_position}"
        outpos = f"Pos:{self.output_position}"
        inphase = f", Phase:{self.phase_condition}" if self.phase_condition else ""
        outphase = f", Phase:{self.output_phase}"
        infactors = f", Factors:{self.subphase_condition}" if self.subphase_condition else ""
        add = f", Add:{self.add_factors}" if self.add_factors else ""
        remove = f", Remove{self.remove_factors}" if self.remove_factors else ""
        return f"\n>PIPE>\n" \
               f"IN | {inpos}{inphase}{infactors}\n" \
               f"OUT| {outpos}{outphase}{add}{remove}"

    def __lt__(self, other):
        if not isinstance(other, Pipe):
            raise TypeError(f"'<' not supported between instances of 'Pipe' and '{type(other).__name__}'")
        return self.subphase_condition < other.subphase_condition
    
class Event:
    """
    Base class for all event types, subclass to create new event types
    
    Parameters:
    position    -- nucleotide position of the event
    probability -- probabilty that a ribosome of suitable state
                   reaching the event will interact with the event

    Attributes:
    pipe        -- the Pipe object associated with the event
    frame       -- the frame of the event, based on position
    """
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
    """Base class for events on the transcript, enforecs probability between 0 and 1"""
    def __init__(self, position: int, probability: float):
        if probability > 1 or probability <= 0:
            raise ValueError(f"Event probability must be greater than 0 and not more than 1, got {probability}")
        super().__init__(position, probability)



# ------------------------------------------------------------------ #
#  Base class for factor behaviours                                    #
# ------------------------------------------------------------------ #

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

class PhaseDecay(FactorBehaviour):
    """Base class for decay in a particular phase"""
    
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

 
class Initiation(Reading):

    @property
    def pipe(self) -> Pipe:
        pc = PhaseCondition(0)
        spc = SubphaseCondition(required="ternary_complex")
        return Pipe(self.frame,
                    self.probability,
                    self.position,
                    phase_condition=pc,
                    subphase_condition=spc,
                    remove_factors='ternary_complex')
    

class IRES(Event):
    @property
    def pipe(self) -> Pipe:
        return Pipe(output_phase=self.frame,
                    probability=self.probability,
                    input_position=self.position,
                    phase_condition=PhaseCondition(-1),
                    subphase_condition=SubphaseCondition('ternary_complex'),
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
        if self.limit is not None:
            return Pipe(output_phase=0,
                        probability=self.probability,
                        input_position=self.position,
                        phase_condition=PhaseCondition(self.frame),
                        subphase_condition=SubphaseCondition('scanning_factors', ceiling={'retained':self.limit}),
                        add_factors='retained'
                        )
        else:
            return Pipe(output_phase=0,
                        probability=self.probability,
                        input_position=self.position,
                        phase_condition=PhaseCondition(self.frame),
                        subphase_condition=SubphaseCondition('scanning_factors'),
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
                    subphase_condition=SubphaseCondition(tuple(('ternary_complex', 'scanning_factors')))
                    )
    
class AllDrop(Reading):
    @property
    def pipe(self) -> Pipe:
        return Pipe(output_phase=-1,
                    probability=self.probability,
                    input_position=self.position)    

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
    
DEFAULT_BEHAVIOUR = [ScanningFactorDissociation(25), TernaryComplexAssociation(20)]


class Transcript(SeqRecord):
    """
    An RNA transcript record that calculates and holds information about Events.

    Inherits all BioPython SeqRecord functionality (reading/writing FASTA,
    GenBank, feature annotation, slicing, …) and adds dynamic RDG specific helpers.

    Parameters
    ----------
    sequence : str | Bio.Seq.Seq
        The RNA sequence. DNA (T) is silently converted to RNA (U).
    id : str
        Transcript identifier (e.g. "ENST00000123456").
    name : str
        Short name (e.g. gene symbol).
    description : str
        Free-text description.
    **kwargs
        Any additional keyword arguments accepted by SeqRecord.

    """
    def __init__(
        self,
        sequence: str | Seq,
        id: str = "<unknown id>",
        name: str = "<unknown name>",
        description: str = "<unknown description>",
        auto = False,
        blank = False,
        *args,
        **kwargs,
    ) -> None:

        raw = str(sequence).upper().replace("T", "U")
        bad = set(raw) - set("AUCGN")
        if bad:
            raise ValueError(f"Invalid RNA characters detected: {bad}")

        # Only set defaults if not already provided
        kwargs.setdefault("id", id)
        kwargs.setdefault("name", name)
        kwargs.setdefault("description", description)

        super().__init__(Seq(raw), *args, **kwargs)

        self._events: list[Event] = []
        if not blank:
            self.add_event(LoadScanning(0, 1))
            self.add_event(AllDrop(len(self), 1))
 
        if auto:
            self.auto_stop_starts()
    
    @property
    def pipes(self, weight_cutoff = 0.0) -> list[Pipe]:
        
        """
        Returns a TransitionMap instance based on the events stored on the transcript
        """
        pipe_list: list[Pipe] = []
        for event in self._events:
            if event.probability > weight_cutoff:
                pipe_list.append(event.pipe)
        return pipe_list


    # ------------------------------------------------------------------
    # Dunder helpers
    # ------------------------------------------------------------------

    def __str__(self) -> str:
        preview = str(self.seq)[:40]
        ellipsis = "..." if len(self) > 40 else ""
        return (
            f"RNASequence(id={self.id!r}, name={self.name!r}, "
            f"len={len(self)}, seq={preview}{ellipsis})"
        )

    def __getitem__(self, index):
        parent = super().__getitem__(index)
        
        new = Transcript(
            str(parent.seq),
            id=parent.id,
            name=parent.name,
            description=parent.description,
        )
        
        return new


    # Public access to event list fo the transcript
    def add_event(self, event:Event):
        """
        Add an event to the transcript, if an identical event exists already at that position,
        it will be silently removed.
        """
        self.remove_event(event)
        self._events.append(event)

    def remove_event(self, event: Event):
        """
        Remove an Event from the transcript with the same type and position as the input.
        """
        for e in self._events:
            if e.type == event.type and e.position == event.position:
                self._events.remove(e)


    def auto_stop_starts(self, cutoff=0.0, reinitiation_prob = 0.0, reinitiation_limit = None):
        """
        Automatically adds initiation and termination events to the transcript based on sequence.
        Uses scores based on Noderer 2014 to calculate probabilities of initation.
        Cutoff determines the minimum probability of initiation of non-AUG codons in order to be considered.
        """
        for i in range(len(self)):
            if self.seq:
                codon = str(self.seq[i:i+3])
                if codon == 'AUG':
                    if 6 < i < len(self) - 5:
                        prob=noderer_start_score(sequence=str(self.seq[i-6:i+5]), aug=True)
                        if prob >= cutoff:
                            self.add_event(Initiation(i+1, prob))
                        
                elif Levenshtein.distance(codon, 'AUG') == 1:
                    if 4 < i < len(self) - 4:
                        prob = noderer_start_score(sequence=str(self.seq[i-4:i+4]), aug=False)
                        if prob >= cutoff:
                            self.add_event(Initiation(i+1, prob))
                            
                elif codon in ['UAG', 'UGA', 'UAA']:
                    self.add_event(Termination(i+1, 1))
                    if reinitiation_prob > 0:
                        self.add_event(Retention(i+1, reinitiation_prob, reinitiation_limit))

    
    def translon_product(self, translon:list[tuple[RiboNode, RiboNode]]):
        if not self.seq:
            raise ValueError(f"Sequence for this transcript does not exist")
        product = ''
        exclude = 0
        for segment in translon:
            start = segment[0].position - 1 + exclude
            exclude = 3
            end   = segment[1].position + 2
            seqment = self.seq[start: end]
            product += Seq.translate(seqment)
        return product
