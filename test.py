from Bio import SeqIO
import src.dynordg as dr
import matplotlib.pyplot as plt
seq = SeqIO.read('gene.fasta', 'fasta')

tcript = dr.Transcript(seq.seq)

tcript.auto_stop_starts(reinitiation_prob=1)


flux = dr.RiboGraphFlux(transition_map=tcript.transition_map,
                            weight_cutoff=0.01,
                            reinitiation_half_life=25,
                            half_life_scanning=1000,
                            ternary_complex_half_life=100,
                            flux_cutoff=0.1,
                            retention_limit=None
                            )



dr.RiboGraphVis(flux, log_scale=10).show()