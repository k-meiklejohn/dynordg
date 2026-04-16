"""
RNASequence - A BioPython-backed class for RNA transcript analysis.

Inherits from Bio.SeqRecord.SeqRecord, so all standard BioPython methods
(IO, alignment, feature annotation, slicing, …) work out of the box.

Custom additions:
    - find_orfs()       : locate all open reading frames
    - translate_orfs()  : translate all ORFs to protein sequences
    - find_motifs()     : search for sequence motifs (string or regex)
    - gc_content()      : GC percentage
    - nucleotide_freq() : per-base frequency counts
"""

from __future__ import annotations
from collections import defaultdict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ..core.events import RiboEvent
from . import TransitionMap
import Levenshtein as lv
from ..core.transitions import RiboTransition
from ...functions import start_score



class Transcript(SeqRecord):
    """
    An RNA transcript record that holds information about RiboEvents.

    Inherits all BioPython SeqRecord functionality (reading/writing FASTA,
    GenBank, feature annotation, slicing, …) and adds RNA-specific helpers.

    Parameters
    ----------
    sequence : str | Bio.Seq.Seq
        The RNA sequence. DNA (T) is silently converted to RNA (U).
    id : str
        Transcript identifier (e.g. "ENST00000123456").
    name : str
        Short name (e.g. gene symbol).
    description : str
        Free-text description.
    **kwargs
        Any additional keyword arguments accepted by SeqRecord.

    """
    def __init__(
        self,
        sequence: str | Seq,
        id: str = "<unknown id>",
        name: str = "<unknown name>",
        description: str = "<unknown description>",
        auto = False,
        *args,
        **kwargs,
    ) -> None:

        raw = str(sequence).upper().replace("T", "U")
        bad = set(raw) - set("AUCGN")
        if bad:
            raise ValueError(f"Invalid RNA characters detected: {bad}")

        # Only set defaults if not already provided
        kwargs.setdefault("id", id)
        kwargs.setdefault("name", name)
        kwargs.setdefault("description", description)

        super().__init__(Seq(raw), *args, **kwargs)

        self.events= defaultdict(lambda: defaultdict(dict))
        if auto:
            self.auto_stop_starts()
    
    def add_event(self, pos: int, type: str, prob: float = 1, drop_prob: float = 0):

        """
        Add event to the transcript of with types: initiation, termination, ires, shift+/-n
        (May also add cap, end but these are added anyway when transition_map is called)
        For termination, prob refers to the probability of 40S retention, i.e. reinitaiton. For all others, \
        it the probability is as normally expected.
        """

        if pos > len(self):
            raise ValueError('')
        self.events[pos][type] = {'probability': prob, 'drop_probability': drop_prob}
    
    @property
    def transition_map(self) -> TransitionMap:
        
        """
        Returns a TransitionMap instance based on the events stored on the transcript
        """

        list_of_events: list[RiboEvent] = []
        list_of_transitions: list[RiboTransition] = []
        for pos in self.events:
            for event in self.events[pos]:
                prob = self.events[pos][event]['probability'] if 'probability' in self.events[pos][event] else 0
                drop_prob = self.events[pos][event]['drop_probability'] if 'drop_probability' in self.events[pos][event] else 0
                list_of_events.append(RiboEvent(pos,
                                                event,
                                                prob,
                                                drop_prob))
        for event in list_of_events:
            list_of_transitions.extend(event._to_transition())
        tmap = TransitionMap()
        tmap.add_weighted_edges_from(list_of_transitions)
        return tmap


    # ------------------------------------------------------------------
    # Dunder helpers
    # ------------------------------------------------------------------

    def __len__(self) -> int:
        return len(self.seq)

    def __str__(self) -> str:
        preview = str(self.seq)[:40]
        ellipsis = "..." if len(self) > 40 else ""
        return (
            f"RNASequence(id={self.id!r}, name={self.name!r}, "
            f"len={len(self)}, seq={preview}{ellipsis})"
        )

    def __getitem__(self, index):
        parent = super().__getitem__(index)
        
        new = Transcript(
            str(parent.seq),
            id=parent.id,
            name=parent.name,
            description=parent.description,
        )
        
        return new


    def auto_stop_starts(self, cutoff=0.001):
        """
        Automatically adds initiation and termination events to the transcript based on sequence.
        Uses scores based on Noderer 2014 to calculate probabilities of initation.
        Cutoff determines the minimum probability of initiation of non-AUG codons in order to be considered.
        """
        for i in range(len(self)):
            codon = str(self.seq[i:i+3])
            if codon == 'AUG':
                if 6 < i < len(self) - 5:
                    prob=start_score(sequence=str(self.seq[i-6:i+5]), aug=True)
                    if prob >= cutoff:
                        self.add_event(i+1, 'initiation', prob)
                    
            elif lv.distance(codon, 'AUG') == 1:
                if 4 < i < len(self) - 4:
                    prob = start_score(sequence=str(self.seq[i-4:i+4]), aug=False)
                    if prob >= cutoff:
                        self.add_event(i+1, 'initiation', prob)
                        
            elif codon in ['UAG', 'UGA', 'UAA']:
                self.add_event(i+1, 'termination', 0, 1)

        self.events[1]['cap'] = {'probability': 1,
                                 'drop_probability': 0}
        self.events[len(self.seq)]['end'] = {'probability': 0,
                                             'drop_probability': 1}        



            

        





