from ..data import SEQUENCE_TO_EFF_AUG, SEQUENCE_TO_EFF_NON_AUG

def start_score(sequence: str, aug: bool):
    """
    Calculate probability of initiation of a start site.
    AUG: length 11
    non-AUG: length 8
    """

    if aug and len(sequence) != 11:
        raise ValueError(f"Expected length 11, got {len(sequence)}: {sequence}")
    if not aug and len(sequence) != 8:
        raise ValueError(f"Expected length 8, got {len(sequence)}: {sequence}")

    table = SEQUENCE_TO_EFF_AUG if aug else SEQUENCE_TO_EFF_NON_AUG

    efficiency = table.get(sequence)
    if efficiency is None:
        return 0

    return min(efficiency * 0.84 / 100, 1)