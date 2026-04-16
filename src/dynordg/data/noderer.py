from pandas import read_csv
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent

AUG_SCORE = read_csv(BASE_DIR / 'aug.csv')
NON_AUG_SCORE = read_csv(BASE_DIR / 'non_aug.csv')

