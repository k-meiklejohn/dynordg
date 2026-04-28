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
        elif len(args) == 2:
            coords = args
        else:
            raise ValueError(f'RiboNode requires 2 ints or a length-2 tuple, got: {args}')

        if len(coords) != 2:
            raise ValueError(f'RiboNode tuple must be of length 2, got length: {len(coords)}')

        x, y = coords

        if not isinstance(x, int) or not isinstance(y, int):
            raise ValueError("RiboNode coordinates must be ints")

        if not (-1 <= y < 4):
            raise ValueError("position must be between -1 and 3")

        return super().__new__(cls, (x, y))

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

    def __repr__(self):
        return f"(Pos:{self.position}, Phase:{self.phase})"