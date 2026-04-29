# DYNORDG

Dynamic Ribosome Decision Graphs (RDGs) for simulating and visualizing ribosome flux along transcripts.

## Overview

A Ribosome Decision Graph (RDG) models the possible paths a ribosome can take along an mRNA transcript.

**Dynamic RDGs** extend this by:
- Representing ribosome flux using edge thickness
- Implicitly encoding overlapping translons via flow rather than explicit separation
- Modeling ribosomal phase states based on downstream potential

This package provides tools to:
- Build RDGs from transcript sequences, or from user defined phase transistions
- Simulate ribosome movement
- Render dynamic flux graphs

## Installation

```bash
pip install dynordg
```

## Quick Start

```python
from dynordg import Transcript, RiboGraphFlux, RiboGraphVis

# Create transcript
t = Transcript("AUGGCCAUGGCGCCCAGAACUGGGUAA")

# Automatically detect start/stop events
t.auto_stop_starts()

# Build flux graph
graph = RiboGraphFlux(t.transition_map())

# Create render object
plot = RiboGraphVis(graph)

# Render
plot.show()
```

## Example Output

Below is example dynamic RDG (if not a realistic one):

![Example RDG](docs/example.png)




