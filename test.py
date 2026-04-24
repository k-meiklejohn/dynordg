from Bio import SeqIO
import src.dynordg as dr
import matplotlib.pyplot as plt
import networkx as nx
seq = SeqIO.read('gene.fasta', 'fasta')

tcript = dr.Transcript(seq.seq)

tcript.auto_stop_starts(reinitiation_prob=1)

    
flux = dr.RiboGraphFlux(transition_map=tcript.transition_map,
                            weight_cutoff=0.01,
                            reinitiation_half_life=25,
                            half_life_scanning=None,
                            half_life_translation=None,
                            ternary_complex_half_life=100,
                            flux_cutoff=0.01,
                            retention_limit=1
                            )

# nx.draw(flux, with_labels=True)
# plt.show()

dr.RiboGraphVis(flux, log_scale=10).show()