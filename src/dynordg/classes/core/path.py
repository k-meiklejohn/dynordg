from typing import Iterable
from .nodes import RiboNode


class RiboPath(list):
    """
"""

    def __init__(self, iterable: Iterable[tuple[RiboNode, RiboNode]]) -> None:
        super().__init__(iterable)
        self.valid()

    def valid(self):
        for i in self:
            if not isinstance(i, tuple) :
                raise ValueError(f'RiboPath must be an interable of RiboNode 2-tuples, got {i}')
            elif len(i) != 2:
                raise ValueError(f'Ribopath tuple must be length 2, got {len(i)}')
            for j in i:
                if not isinstance(j, RiboNode):
                    raise ValueError(f'RiboPath tuples must contain only 2 RiboNodes, got {j}')

        for i in range(len(self)):
            if i == len(self)- 1:
                continue
            if self[i][1] != self[i+1][0]:
                raise ValueError('RiboPath must be continuous, each target must be the the source of the next edge.' \
                f'target: {self[i][1]}, next source: {self[i+1][0]}')