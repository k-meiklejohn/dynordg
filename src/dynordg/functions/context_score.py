from ..data import AUG_SCORE, NON_AUG_SCORE

def start_score(sequence: str, aug: bool):
    df = AUG_SCORE if aug else NON_AUG_SCORE
    
    match = df[df['sequence'] == sequence]
    
    if match.empty:
        return 0  # or raise an error / default value
    
    efficiency = match['efficiency'].values[0]
    return min(efficiency * 0.84/100, 1)
