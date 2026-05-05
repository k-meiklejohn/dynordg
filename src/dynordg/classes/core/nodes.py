class RiboNode(tuple):
    """
    Special 2-tuple that refers to a position in simplified Ribosomal phase space, the first integer \
    refers to the nucleotide position on a transcript, while the second refers to the phase \
    of the ribosome: -1 for not associated, 0 for scanning and 1,2,3 for translating in the frame \
    where frame = position % 3 + 1
    
    """
    def __new__(cls, *args):
        if len(args) == 1 and isinstance(args[0], (tuple, RiboNode)):
            coords = args[0]
        elif len(args) >  1:
            coords = args
        else:
            raise ValueError(f'RiboNode requires 2 ints or a length-2 tuple, got: {args}')

        if not len(coords) >= 2:
            raise ValueError(f'RiboNode tuple must be at least length 2, got length: {len(coords)}')



        if not isinstance(coords[0], int) or not isinstance(coords[1], int):
            raise ValueError("RiboNode coordinates must be ints")


        return super().__new__(cls, coords)

    @property
    def position(self) -> int:
        """
        Nucleotide position of the node
        """
        return self[0]

    @property
    def phase(self) -> int:
        """
        Simplified phase of the node -1 is not reading, 0 is scanning, 1, 2, 3 are translating in one of those frames
        """
        return self[1]
    
    @property
    def factors(self) -> bool:
        """
        Subphase of node:
        No Extra association = False
        Ternary Complex Associated = True
        Scanning factors associated = True
        """
        if len(self) > 2:
            return self[2]
        else:
            return False

    def __repr__(self):
        return f"(Pos:{self.position}, Phase:{self.phase}, F:{self.factors})"
    
    @property
    def simple(self):
        return RiboNode(self.position, self.phase)