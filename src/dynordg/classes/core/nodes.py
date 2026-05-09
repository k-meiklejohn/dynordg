class State:
    def __init__(self, phase: int, *args:str):

        if not isinstance(phase, int):
            raise ValueError(f'RiboNode phase must be int, got {phase}')
        if phase < -1 or phase > 3:
            raise ValueError(f'Phase must be int, between -1 and 3 inclusive, got {phase}')
        for arg in args:
            if not isinstance(arg, str):
                raise ValueError
        self.phase = phase
        self.factors: tuple[str,...] = tuple({*args})

    def __repr__(self):
        if len(self.factors):
            return f"Phase:{self.phase} Factors:{self.factors}"
        else:
            return f"Phase:{self.phase}"
        
    def __eq__(self, value):
        if not isinstance(value, State):
            return False
        phases = self.phase == value.phase
        factors = self.factors == value.factors
        return factors and phases
    
    def __add__(self, *other):
        for arg in other:
            if not isinstance(other, str):
                raise ValueError(f"+ not implemented between 'State' and '{type(arg).__name__}'")
        new_factors = self.factors + other
        print(new_factors)
        print(self.phase)
        return State(self.phase, *new_factors)
    
    def __hash__(self) -> int:
        return hash((self.phase, self.factors))
    
    def __contains__(self, item):
        return item in self.factors


class RiboNode:
    def __init__(self, position: int, state: State):
        if not isinstance(position, int):
            raise ValueError(f'RiboNode position must be int, got {position}')


        self.position = position
        self.state = state
        
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
        if not isinstance(value, RiboNode):
            return False
        pos = self.position == value.position
        phases = self.phase == value.phase
        factors = self.factors == value.factors
        return factors and phases and pos
    
    def __hash__(self) -> int:
        return hash((self.position, self.phase, self.factors))
    
    def __add__(self, *other):
        return RiboNode(self.position, self.state + other)

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
