from importlib.resources import files, as_file
import csv

def load_score_dict(filename: str) -> dict:
    with as_file(files("dynordg.data") / filename) as path:
        with open(path, newline="") as f:
            reader = csv.DictReader(f)
            return {row["sequence"]: float(row["efficiency"]) for row in reader}
        

SEQUENCE_TO_EFF_AUG = load_score_dict("aug.csv")
SEQUENCE_TO_EFF_NON_AUG = load_score_dict("non_aug.csv")

