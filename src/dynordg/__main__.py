import csv
from Bio import SeqIO
import argparse as ap
from .core import Event, EVENT_REGISTRY, Transcript, TernaryComplexAssociation, ScanningFactorDissociation, ScanningDecay, TranslationDecay
from .graph import RiboSkeleton, FluxGraph
from .render import RiboGraphVis
parser = ap.ArgumentParser(
    description='Generate a Dynamic RDG plot from a csv of events',
)

parser.add_argument('filename', nargs='?', default=None, help="A CSV file of the format position,event,probability")
parser.add_argument('-l', '--transcript_length', help="Integer | The end of the transcript, default is the largest event position +10")
parser.add_argument('-o', '--log_scale', default=0, help='Integer | Factor to reduce the disctance between nodes logarithmically  (larger gaps become shorter by a greater amount)')
parser.add_argument('-e', '--loading_efficiency', help="Float | The loading efficiency of the cap, i.e. the arbitrarty proportion of ribosomes that load at position 0")
parser.add_argument('-a', '--tc_association_halflife', help='Integer | The distance in nucleotides for half of a given proportion of scanning ribosomes to reassociate with the ternary complex')
parser.add_argument('-d', '--sf_dissociation_halflife', help='Integer | The distance, measured in nucleotides, for half a given proportion of ribosomes to lose scanning factors after initiation')
parser.add_argument('-t', '--translation_decay', help='Integer | The halflife of dissociation of translating ribosomes measured in nucleotides')
parser.add_argument('-s', '--scanning_decay', help='Integer | The halflife of dissociation of scanning ribosomes measured in nucleotides')
parser.add_argument('-c', '--flux_cutoff', help='Float | The level below which the flux of a given edge will be distributed to its sisters')
parser.add_argument('--fasta', help="File |  A fasta file to associate to the transcript, overrides --length option")
parser.add_argument('-g', '--guess', action='store_true', help='Boolean |EXPERIMENTAL, requires --fasta: if used, the probabilities of initiations will be guessed based on data from Noderer et al. 2014 and Diaz de Arce et al. 2018')
parser.add_argument('-r', '--reinit_prob', default=0, help='Float | EXPERIMENTAL:¸requires --fasta Probability of retention of 40S subunit after elongation termination, to be used when assigning probabilities at stop codons')
args = parser.parse_args()

if args.filename is None and args.fasta is None:
    parser.error("filename is required unless --fasta is specified")

def parse_event(pos: int, event_str: str, prob: float) -> Event:
    for cls in EVENT_REGISTRY:
        cls : Event
        m = cls.pattern.match(event_str)
        if m:
            return cls.from_match(m, pos, prob)
    raise ValueError(f"Unrecognized event string: {event_str!r}")


events = []
max_pos = 0

if args.filename:
    with open(args.filename) as paramfile:
        reader = csv.reader(paramfile)
        first = True
        for row in reader:
            if first:
                first = False
                if row == ["position", "event", "probability"]:
                    continue
            pos = int(row[0])
            max_pos = max(max_pos,pos)
            event_str = str(row[1])
            prob = float(row[2])
            events.append(parse_event(pos, event_str, prob))

if args.fasta:
    record = SeqIO.read(args.fasta, format='fasta')
    transcript = Transcript(record.seq)
    if args.guess:
        transcript.auto_stop_starts(reinitiation_prob=args.reinit_prob)
else:
    length = args.transcript_length if args.transcript_length else max_pos + 10
    transcript = Transcript("N"*length)

for event in events:
    transcript.add_event(event)

behaviours = []
scan_decay = ScanningDecay(int(args.scanning_decay)) if args.scanning_decay else None
scan_fact  = ScanningFactorDissociation(int(args.sf_dissociation_halflife)) if args.sf_dissociation_halflife else None
tran_decay = TranslationDecay(int(args.translation_decay)) if args.translation_decay else None
tc_assoc   = TernaryComplexAssociation(int(args.tc_associaiton_halflife)) if args.tc_association_halflife else None

for i in [scan_decay,scan_fact,tran_decay,tc_assoc]:
    if i:
        behaviours.append(i)


skel = RiboSkeleton(transcript.pipes, behaviours=behaviours)

flux_cutoff = args.flux_cutoff if args.flux_cutoff else 0
graph = FluxGraph(skel, flux_cutoff=float(flux_cutoff))

plot = RiboGraphVis(graph, log_scale=float(args.log_scale))

plot.show()