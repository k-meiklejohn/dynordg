class State:
    def __init__(self, phase: int, **subphases):

        if not isinstance(phase, int):
            raise ValueError(f'RiboNode phase must be int, got {phase}')
        if phase < -1 or phase > 3:
            raise ValueError(f'Phase must be int, between -1 and 3 inclusive, got {phase}')
        self.phase = phase
        self._subphases = {}
        self._subphases.update(subphases)

    @property
    def subphases(self):
        return tuple(sorted(self._subphases.items()))
    
    @property
    def __members(self):
        return (self.phase, self.subphases)

    def __repr__(self):
        if len(self.subphases):
            return f"Phase:{self.phase} Factors:{self.subphases}"
        else:
            return f"Phase:{self.phase}"
        
    def __eq__(self, value):
        if isinstance(value, State):
            return self.__members == value.__members
        else:
            return False
        
    def add_subphase(self, key:str, value:... = True) -> 'State':
        """Returns a new state with the subphase added, will overwrite existing values silently"""
        out_dict = dict(self.subphases)
        out_dict[key] = value
        return State(self.phase, **out_dict)

    def __add__(self, other):
        other = str(other)
        out_dict = dict(self.subphases)
        if other in self._subphases:
            if isinstance(self._subphases[other], int):
                out_dict[other] += 1
                return State(self.phase, **out_dict)
    
        out_dict[other] = True
        return State(self.phase, **out_dict)

    def __sub__(self, other):
        other = str(other)
        out_dict = dict(self.subphases)
        print(out_dict)
        if other in self._subphases:
            if isinstance(self._subphases[other], int):
                out_dict[other] -= 1
                if out_dict[other] == 0:
                    out_dict.pop(other)
                return State(self.phase, **out_dict)
        
            out_dict.pop(other)
            return State(self.phase, **out_dict)
        return self
    
    def __hash__(self) -> int:
        return hash(self.__members)
    
    def __contains__(self, item):
        if not isinstance(item, str):
            return TypeError(f"'in' not supported between {type(item).__name__} and 'State")
        return item in self.factors
    
    def __lt__(self, other):
        if isinstance(other, State):
            if self == other:
                return False
            elif self.phase < other.phase:
                return True
            elif self.factors < other.factors:
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
            elif self.factors > other.factors :
                return True 
            else:
                return False
        else:
            raise TypeError(f"Operand > not supported between 'State' and '{type(other).__name__}'")

class RiboNode:
    def __init__(self, position: int, state: State):

        if not isinstance(position, int):
            raise TypeError(f'RiboNode position must be int, got {position}')

        self.position:int = position
        self.state:State = state
    
    @property
    def __members(self):
        return (self.position, self.state)
        


    @property
    def phase(self) -> int:
        """
        Simplified phase of the node -1 is in solution, 0 is scanning, 1, 2, 3 are translating in one of those frames
        """

        return self.state.phase

    
    @property
    def factors(self) -> tuple[str, ...]:
        return self.state.factors

    @property
    def simple(self):
        return RiboNode(self.position, State(self.phase))

    def __repr__(self):
        if len(self.state.factors):
            return f"(Pos:{self.position}, Phase:{self.phase}, F:{self.factors})"
        else:
            return f"(Pos:{self.position}, Phase:{self.phase})"
        
    def __eq__(self, value):
        if isinstance(value, RiboNode):
            return self.__members == value.__members
        else:
            return False

    
    def __hash__(self) -> int:
        return hash(self.__members)
    
    def __add__(self, *other):
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

class Ribosome:
    UPSTREAM_LENGTH = 15
    DOWNSTREAM_LENGTH = 15

    def __init__(self, coords:RiboNode) -> None:
        self.coords = coords

    def advance(self, distance):
        self.coords = RiboNode(self.position + distance, self.coords.state)

    def initiate(self):
        pass
    
    @property
    def position(self) -> int:
        return self.coords.position
    
    @property
    def heel(self) -> int:
        return self.position - self.DOWNSTREAM_LENGTH
    
    @property
    def toe(self) -> int:
        return self.position + self.UPSTREAM_LENGTH


    def clashes(self, other) -> bool:
        if not isinstance(other, Ribosome):
            raise ValueError(f"Only Ribosomes may clash, got type '{type(other).__name__}'")
        return self.toe >= other.heel and self.heel <= other.heel \
                or self.heel <= other.toe and self.toe >= other.toe
