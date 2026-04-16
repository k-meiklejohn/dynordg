from ..data import AUG_SCORE, NON_AUG_SCORE

def start_score(sequence: str, aug: bool):
    """
    Calculate the probabiltily of initation of a start site. If aug is True, 
    """

    if len(sequence) != 11 and aug:
        raise ValueError(f'Sequence must be length 11 when aug is True, got {len(sequence)}: {sequence}')
    elif len(sequence) != 8 and not aug:
        raise ValueError(f'Sequence must be length 8 when aug is false, got {len(sequence)}: {sequence}')
    
    df = AUG_SCORE if aug else NON_AUG_SCORE

    
    match = df[df['sequence'] == sequence]
    
    if match.empty:
        return 0  # or raise an error / default value
    
    efficiency = match['efficiency'].values[0]
    return min(efficiency * 0.84/100, 1)
