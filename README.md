# DYNORDG

## Description

A [ribosome decision graph (RDG)](https://www.genome.org/cgi/doi/10.1101/gr.278810.123) represents the paths of the ribosome as it moves along a transcript. At points where the ribosome can alter its state, i.e. when initiating, the path branches and, representing the two different states that the ribosome can now occupy. 

Dynamic RDGs build on this idea by representing the proportional flux of the ribosome through a particular path by the thickness of the path. However, there is another important difference with traditional RDGs. While an the traditional RDG maintains an explicit distinction between different overlapping translons, a dynamic RDG implies their existence by the varaible thickness of the paths. To clarify, the phase of two ribosomes is identical if they have the same downstream potential, i.e. two scanning ribosomes can both recognise a start codon and initiate, changing from the scanning phase to the translating phase. While a translating ribosome cannot initiate at this start codon, and therefore lacks the potential and is subsequently in a different phase.

This package provides a way to build, simulate and render dynamic RDGs

## Usage

There are 4 main classes that dynordg uses in order to generate and render graphs. 

### class Transcript(SeqRecord)

This class inherits properties from biopython's SeqRecord and acts as a way to store information about a transcript that will be used to generate a graph.

#### Input

The main required input for the Transcript class is a nucleotide sequence. T's will be silently converted to U's'. 

#### Properties
events - a dict in a dict in a dict structure contain all events on the transcript index by position and the type, and finally with probability/drop_probability.

#### Methods
auto_stop_starts() - finds all AUGs and near cognate start codons and calculates a probability score for them based on kozak scores from [Noderer 2014](https://link.springer.com/article/10.15252/msb.20145136), as well as stop codons whose drop_probablity is set to 1.

#### 