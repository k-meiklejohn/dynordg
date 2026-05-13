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
from typing import Any, NoReturn
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ..core.events import Event, Initiation, Termination, Retention, LoadScanning, AllDrop, Transition, Pipe
from . import TransitionMap
import Levenshtein as lv
from ...functions import start_score

class Transcript(SeqRecord):
    """
    An RNA transcript record that holds information about RiboEvents.

    Inherits all BioPython SeqRecord functionality (reading/writing FASTA,
    GenBank, feature annotation, slicing, …) and adds dynamic RDG specific helpers.

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
        blank = False,
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

        self._events: list[Event] = []
        if not blank:
            self.add_event(LoadScanning(0, 1))
            self.add_event(AllDrop(len(self), 1))
 
        if auto:
            self.auto_stop_starts()
    
    @property
    def pipes(self, weight_cutoff = 0.0) -> list[Pipe]:
        
        """
        Returns a TransitionMap instance based on the events stored on the transcript
        """
        pipe_list: list[Pipe] = []
        for event in self._events:
            if event.probability > weight_cutoff:
                pipe_list.append(event.pipe)
        return pipe_list


    # ------------------------------------------------------------------
    # Dunder helpers
    # ------------------------------------------------------------------

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


    # Public access to event list fo the transcript
    def add_event(self, event:Event):
        """
        Add an event to the transcript, if an identical event exists already at that position,
        it will be silently removed.
        """
        self.remove_event(event)
        self._events.append(event)

    def remove_event(self, event: Event):
        """
        Remove an Event from the transcript with the same type and position as the input.
        """
        for e in self._events:
            if e.type == event.type and e.position == event.position:
                self._events.remove(e)


    def auto_stop_starts(self, cutoff=0.0, reinitiation_prob = 0.0):
        """
        Automatically adds initiation and termination events to the transcript based on sequence.
        Uses scores based on Noderer 2014 to calculate probabilities of initation.
        Cutoff determines the minimum probability of initiation of non-AUG codons in order to be considered.
        """
        for i in range(len(self)):
            if self.seq:
                codon = str(self.seq[i:i+3])
                if codon == 'AUG':
                    if 6 < i < len(self) - 5:
                        prob=start_score(sequence=str(self.seq[i-6:i+5]), aug=True)
                        if prob >= cutoff:
                            self.add_event(Initiation(i+1, prob))
                        
                elif lv.distance(codon, 'AUG') == 1:
                    if 4 < i < len(self) - 4:
                        prob = start_score(sequence=str(self.seq[i-4:i+4]), aug=False)
                        if prob >= cutoff:
                            self.add_event(Initiation(i+1, prob))
                            
                elif codon in ['UAG', 'UGA', 'UAA']:
                    self.add_event(Termination(i+1, 1))
                    if reinitiation_prob > 0:
                        self.add_event(Retention(i+1, reinitiation_prob))

  



            

        





