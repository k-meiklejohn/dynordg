from importlib.resources import files, as_file
import csv

def load_score_dict(filename: str) -> dict:
    with as_file(files("dynordg.data") / filename) as path:
        with open(path, newline="") as f:
            reader = csv.DictReader(f)
            return {row["sequence"]: float(row["efficiency"]) for row in reader}
        
# The AUG and NON-AUG data are derived from:
# Noderer, W. L. et al. Quantitative analysis of mammalian translation
#  initiation sites by FACS‐seq. Mol. Syst. Biol. 10, 748 (2014) and
# Diaz de Arce, A. J., Noderer, W. L. & Wang, C. L. Complete motif analysis of 
# sequence requirements for translation initiation at non-AUG start codons. 
# Nucleic Acids Res. 46, 985–994 (2018) respectively.

SEQUENCE_TO_EFF_AUG = load_score_dict("aug.csv")
SEQUENCE_TO_EFF_NON_AUG = load_score_dict("non_aug.csv")

