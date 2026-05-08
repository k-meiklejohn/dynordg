class State:
    def __init__(self, phase, *args:str):

        if not isinstance(phase, int):
            raise ValueError(f'RiboNode phase must be int, got {phase}')
        if phase < -1 or phase > 3:
            raise ValueError(f'Phase must be int, between -1 and 3 inclusive, got {phase}')
        for arg in args:
            if not isinstance(arg, str):
                raise ValueError
        self.phase = phase
        self.factors = args

    def __repr__(self):
        if len(self.factors):
            return f"Phase:{self.phase} Factors:{self.factors}"
        else:
            return f"Phase:{self.phase}"
        
    def __eq__(self, value: State):
        if not isinstance(value, State):
            return False
        phases = self.phase == value.phase
        factors = self.factors == value.factors
        return factors and phases



        

class RiboNode:
    """
    Special 2-tuple that refers to a position in simplified Ribosomal phase space, the first integer \
    refers to the nucleotide position on a transcript, while the second refers to the phase \
    of the ribosome: -1 for not associated, 0 for scanning and 1,2,3 for translating in the frame \
    where frame = position % 3 + 1
    
    """
    def __init__(self, position: int, state: State):
        if not isinstance(position, int):
            raise ValueError(f'RiboNode position must be int, got {position}')


        self.coords: RiboNode = (position, state.phase)
        self.state = state
        
    @property
    def position(self) -> int:
        """
        Nucleotide position of the node
        """
        return self.coords[0]

    @property
    def phase(self) -> int:
        """
        Simplified phase of the node -1 is in solution, 0 is scanning, 1, 2, 3 are translating in one of those frames
        """
        return self.coords[1]
    
    @property
    def factors(self) -> list[str]:
        return self.state.factors
    
    @property
    def simple(self):
        return RiboNode(self.position, self.phase)
    
    def __repr__(self):
        if len(self.state.factors):
            return f"(Pos:{self.position}, Phase:{self.phase}, F:{self.factors})"
        else:
            return f"(Pos:{self.position}, Phase:{self.phase})"
        
    def __eq__(self, value: RiboNode):
        if not isinstance(value, RiboNode):
            return False
        pos = self.position == value.position
        phases = self.phase == value.phase
        factors = self.factors == value.factors
        return factors and phases and pos