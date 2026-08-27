# DynoRDG API Reference

DynoRDG builds and renders **Dynamic Ribosome Decision Graphs (RDGs)** — models
of the paths a ribosome can take along an mRNA transcript, with flux
(proportion of ribosomes) encoded as edge thickness.

This reference covers the public API (everything importable from `dynordg`
directly) plus the extensibility hooks intended for users who want to define
custom initiation/termination behaviour or customize rendering.

## Contents

- [Installation & Quick Start](#installation--quick-start)
- [Conceptual model](#conceptual-model)
- [Public API](#public-api)
  - [`Transcript`](#transcript)
  - [`State`](#state)
  - [`RiboNode`](#ribonode)
  - [`Transition`](#transition)
  - [Event types](#event-types)
  - [Factor behaviours](#factor-behaviours)
  - [`RiboSkeleton`](#riboskeleton)
  - [`FluxGraph` / `SimpleFluxGraph`](#fluxgraph--simplefluxgraph)
  - [`RiboGraphVis`](#ribographvis)
  - [Convenience functions](#convenience-functions)
- [Extensibility](#extensibility)
  - [Custom event types](#custom-event-types)
  - [Custom factor behaviours](#custom-factor-behaviours)
  - [Customizing rendering](#customizing-rendering)
- [Appendix: phase & frame reference](#appendix-phase--frame-reference)

---

## Installation & Quick Start

```bash
pip install dynordg
```

```python
from dynordg import quick_plot

plot = quick_plot("AUGGCCAUGGCGCCCAGAACUGGGUAA")
plot.show()
```

`quick_plot` is a shortcut over the full pipeline: it builds a `Transcript`,
auto-detects start/stop codons, propagates flux through a `RiboSkeleton` into
a `FluxGraph`, and renders it with `RiboGraphVis`. The sections below cover
each stage individually so you can customize any part of it.

---

## Conceptual model

Every position in riosomal phase space is represented as a **`RiboNode`**: a
 `(position, State)` pair.

- `position` — 1-indexed nucleotide position on the transcript. Position `-1`
  is reserved for the **bulk node**, the cytoplasmic pool ribosomes occupy
  when off the mRNA.
- `State` — carries a `phase` and optional keyword **subphases** (factors).
  - `phase = -1` → in solution (bulk pool)
  - `phase = 0`  → scanning (40S subunit)
  - `phase = 1, 2, 3` → elongating (80S) in the corresponding reading frame
  - Subphases such as `ternary_complex` and `scanning_factors` track which
    initiation/scanning factors a ribosome currently carries.

A **`Pipe`** connects a matching `(phase, subphase)` entry condition at one
position to an output `(phase, subphase)` at another position — it's how an
**`Event`** (start codon, stop codon, IRES, frameshift, etc.) is translated
into graph transitions. A **`FactorBehaviour`** additionally models gradual
factor gain/loss (e.g. re-acquiring ternary complex while scanning) that
happens continuously rather than at a single event.

The pipeline is:

```
Transcript (+ Events)
      │  .pipes
      ▼
RiboSkeleton (+ FactorBehaviours)   — all theoretically possible transitions
      │
      ▼
FluxGraph                            — actual flux propagated from the bulk pool
      │  .simple
      ▼
RiboGraphVis                         — rendered matplotlib Figure
```

---

## Public API

### `Transcript`

`dynordg.Transcript(sequence, id="<unknown id>", name="<unknown name>", description="<unknown description>", auto=False, blank=False, *args, **kwargs)`

An RNA transcript record. Subclasses BioPython's `SeqRecord`, so all standard
`SeqRecord` functionality (FASTA/GenBank I/O, slicing, annotation) is
available in addition to the members below.

- **Sequence handling**: input is uppercased and DNA `T` is silently
  converted to RNA `U`. Any character outside `AUCGN` raises `ValueError`.
- **`blank`**: if `False` (default), the transcript automatically gets a
  `LoadScanning(0, 1)` event (ribosomes enter scanning at position 0) and an
  `AllDrop(len(self), 1)` event (anything reaching the end drops to bulk).
  Pass `blank=True` to build the event list entirely yourself.
- **`auto`**: if `True`, calls `auto_stop_starts()` immediately after
  construction.

**Attributes / properties**

| Member | Description |
|---|---|
| `pipes` | Property. Returns the `list[Pipe]` derived from every added `Event`. *Note: because this is a property, not a method, the `weight_cutoff` parameter defined on it cannot actually be passed — it always uses the default `0.0`.* |

**Methods**

- `add_event(event: Event)` — Adds an event. If an event of the same type
  already exists at that position, it is silently replaced.
- `remove_event(event: Event)` — Removes any event matching the given
  event's type and position.
- `auto_stop_starts(cutoff=0.0, reinitiation_prob=0.0, reinitiation_limit=None)`
  — Scans the sequence and automatically adds `Initiation` events at AUG and
  near-cognate (single substitution) start codons, scored using the Noderer
  et al. (2014) / Diaz de Arce et al. (2018) initiation-efficiency tables,
  and `Termination` events at in-frame stop codons (`UAG`, `UGA`, `UAA`).
  - `cutoff` — minimum initiation probability required to add an
    `Initiation` event.
  - `reinitiation_prob` — if `> 0`, also adds a `Retention` event at every
    stop codon with this probability.
  - `reinitiation_limit` — passed through as the `limit` on those
    `Retention` events (max prior initiations before reinitiation stops).
- `translon_product(translon: list[tuple[RiboNode, RiboNode]]) -> Seq` —
  Translates the amino acid sequence for a translon (a list of
  `(start_node, end_node)` segment pairs, as produced by
  `FluxGraph.translons`).

### `State`

`dynordg.State(phase: int, **subphases)`

An immutable, hashable representation of a ribosome's phase and factor
composition.

- `phase` — must be an `int` in `[-1, 3]`.
- `subphases` — arbitrary keyword flags/counters, e.g.
  `State(0, ternary_complex=True, retained=2)`.

**Attributes**: `phase`, `subphases` (canonical sorted tuple of the keyword
items).

**Operators**

| Operator | Behaviour |
|---|---|
| `state + "factor"` | Returns a new `State` with `"factor"` added. If the existing value is an `int`, increments it instead. |
| `state - "factor"` | Returns a new `State` with `"factor"` removed (or decremented, dropped once it hits 0). No-op if absent. |
| `"factor" in state` | `True` if `"factor"` is a key in the subphases. |
| `<`, `>` | Orders first by `phase`, then by `subphases`. |
| `==`, `hash()` | Structural equality/hash on `(phase, subphases)`. |

### `RiboNode`

`dynordg.RiboNode(position: int, state: State)`

A node in ribosomal phase space — the fundamental unit of every graph in
this package.

**Properties**

- `phase -> int` — shortcut for `state.phase`.
- `subphases -> tuple` — shortcut for `state.subphases`.
- `simple -> RiboNode` — a copy with the same `position` and `phase` but an
  empty (subphase-free) `State`. Used to collapse factor bookkeeping detail
  out of a graph (see [`FluxGraph.simple`](#fluxgraph--simplefluxgraph)).

Supports `==`, `hash()`, `<`, `>` (ordered by `position` then `state`), and
`node + "factor"` (shortcut for `RiboNode(position, state + "factor")`).

### `Transition`

`dynordg.Transition(source: RiboNode, target: RiboNode, weight: float)`

A plain dataclass describing a single weighted edge: ribosomes move from
`source` to `target` with a probability defined by `weight`.

### Event types

All event types share the base contract from `Event` (see
[Custom event types](#custom-event-types) for the full base class and how to
add your own). Each event exposes a `.pipe` property that resolves it into
one or more `Pipe` transitions when the transcript's `.pipes` is evaluated.

| Class | Base | Constructor | Behaviour |
|---|---|---|---|
| `Initiation(position, probability)` | `Reading` (`0 < probability <= 1`) | start codon | 40S → 80S: requires phase `0` **and** `ternary_complex`; consumes `ternary_complex`; outputs into the reading frame implied by `position`. |
| `Termination(position, probability)` | `Reading` | stop codon | 80S → bulk: requires the ribosome to be translating in the matching frame; drops **all** factors. |
| `Retention(position, probability, limit=None)` | `Reading` | reinitiation after stop | 80S → 40S (no positional movement): requires `scanning_factors`; adds a `retained` counter. `limit`, if set, caps how many times a ribosome may have already retained (via a `ceiling` on `retained`). |
| `AllDrop(position, probability)` | `Reading` | end-of-transcript catch-all | Any phase → bulk, unconditionally. Added automatically at `len(transcript)` unless `Transcript(..., blank=True)`. |
| `LoadScanning(position, probability)` | `Event` (any `probability >= 0`) | scanning entry point | Bulk → scanning (phase `0`): requires **both** `ternary_complex` and `scanning_factors` present on the bulk-pool node. Added automatically at position `0` unless `blank=True`. |
| `IRES(position, probability)` | `Event` | internal ribosome entry | Bulk → translating (cap-independent initiation): requires `ternary_complex`; consumes it; outputs into the frame implied by `position`. |
| `Frameshift(position, probability, amount)` | `Reading` | frameshift site | position, current frame x →  position + amount, frame y where y = (position + amount) % 3 | 

Every event's `.frame` property derives the reading frame from its
`position % 3` (see the [appendix](#appendix-phase--frame-reference)).

### Factor behaviours

Factor behaviours model gradual, distance-dependent factor changes that
happen *between* events rather than at a single position. Pass a list of
them to `RiboSkeleton(..., behaviours=[...])`.

| Class | `priority` | Fires when | Effect |
|---|---|---|---|
| `ScanningDecay(half_life)` | 0 | phase `0` → phase `0` | Diverts a fraction of flux to the bulk pool, modelling ribosome loss during scanning. Fraction = `1 - 0.5 ** (distance / half_life)`. |
| `TranslationDecay(half_life)` | 0 | phase `> 0` → same phase | Same decay formula, applied during elongation. |
| `TernaryComplexAssociation(half_life)` | 10 | phase `0` → phase `0`, `ternary_complex` not yet present | Adds `ternary_complex` to a fraction of flux. If `half_life` is `None`, the *entire* fraction re-associates instantly; otherwise it follows the same exponential formula. |
| `ScanningFactorDissociation(half_life)` | 20 | phase `> 0` → same phase, `scanning_factors` present | Removes `scanning_factors` from a fraction of flux (`0` if `half_life` is `None`, i.e. never). |

`half_life` is passed positionally as the first constructor argument, e.g.
`ScanningDecay(50)`.

**`DEFAULT_BEHAVIOUR`** — the list used internally by `quick_graph`:

```python
DEFAULT_BEHAVIOUR = [ScanningFactorDissociation(25), TernaryComplexAssociation(20)]
```

Lower `priority` behaviours are applied first at each node when a graph is
constructed; `PhaseDecay` subclasses (structural decay) always run before
factor gain/loss.

### `RiboSkeleton`

`dynordg.RiboSkeleton(pipelist: list[Pipe] | None = None, behaviours: list[FactorBehaviour] | None = None, **attr)`

The **theoretical map** of every reachable `RiboNode` and weighted
`Transition`, built by walking outward from each pipe's entry point and
applying the supplied factor behaviours at every step. This is an
intermediate structure — most users go straight from `Transcript.pipes` to
`RiboSkeleton` to `FluxGraph`, but `RiboSkeleton` is exposed so you can:

- Inspect/modify the full transition set before propagating flux.
- Call `smooth_all_weights(by_factor=None)` or
  `smooth_node_weight(node, by_factor=None)` to redistribute a node's
  outgoing weight more evenly across its edges (`by_factor=None` sets all
  equal; otherwise it partially interpolates toward equal weighting).
- Manually add transitions with `add_transition(transition)` /
  `add_transitions_from(transitions)`.

Raises `RuntimeError` if any node's flux is asked to continue indefinitely
(no matching downstream pipe) or if total outgoing weight from a node
exceeds `1`.

### `FluxGraph` / `SimpleFluxGraph`

`dynordg.FluxGraph(skeleton: RiboSkeleton | None = None, flux_cutoff: float = 0.0, **attr)`

Propagates actual flux (starting with `1.0` at the bulk node) through a
`RiboSkeleton` in topological order, pruning any edge whose flux falls below
`flux_cutoff` (redistributing it proportionally among surviving edges at
that node) and normalizing so the maximum edge flux is `1.0`.

**Properties**

- `ribopaths -> list[list[RiboNode]]` — every simple path representing
  continuous 40S association with the mRNA, from a loading node to the bulk
  node.
- `translons -> list[list[tuple[RiboNode, RiboNode]]]` — sub-segments of
  `ribopaths` representing continuous 80S/translating association; each
  translon is a list of `(start, end)` node pairs (per contiguous frame),
  suitable for `Transcript.translon_product`.
- `simple -> SimpleFluxGraph` — collapses subphase detail (via
  `RiboNode.simple`), merges pass-through nodes, and produces a
  decay-annotated graph ready for rendering. This is what `RiboGraphVis`
  consumes internally.
- `weight_skeleton -> RiboSkeleton` — rebuilds a `RiboSkeleton` whose edge
  weights are this graph's *conditional* `edge_weight` values rather than
  absolute flux.

**Methods**

- `flux_proportion(u, v) -> float` — total proportion of bulk flux
  travelling from `u` to `v`, summed over all simple paths between them.
- `flux_proportion_path(path: list[RiboNode]) -> float` — flux proportion
  along one specific path.
- `node_flux(node) -> float` — total outgoing flux from a node.
- `edge_weight(u, v) -> float` — fraction of `u`'s outgoing flux carried by
  the edge to `v`.
- `prune_recycle_edges()` — removes bulk-to-bulk edges (and any nodes left
  isolated by their removal); called automatically before rendering.

A `UserWarning` is raised if inbound and outbound flux at the bulk node
disagree by more than `flux_error` (`1e-15`), which can happen from
floating-point accumulation in complex graphs and is generally safe to
ignore.

`SimpleFluxGraph` is the same structure with `edge_weight`/`node_flux`
computed from `flux_start` rather than `flux`; you normally only encounter
it via `FluxGraph.simple`. It is necessary for renderings to be performed with
this structure but that transformation is performed automatically when calling RiboGraphVis.

### `RiboGraphVis`

`dynordg.RiboGraphVis(graph_to_render: FluxGraph | SimpleFluxGraph, fig_size=(12, 6), dpi=150, log_scale=1, engine=None, renderer=None, **attr)`

Renders a `FluxGraph` (or `SimpleFluxGraph`) as a matplotlib `Figure`,
computing layout automatically on construction.

- `log_scale` — base for logarithmic compression of the x-axis between
  distant nodes. `<= 1` disables compression and uses raw nucleotide
  positions.
- `engine` / `renderer` — override points; see
  [Customizing rendering](#customizing-rendering).

**Methods / properties**

- `compute_layout()` — re-runs the layout pipeline and re-renders `self.fig`.
  Call again after manually editing graph edges.
- `get_figure() -> Figure`
- `show()` — displays inline in a Jupyter/IPython environment, or via
  `plt.show(block=True)` in a script.
- `save(filename="output.png", dpi=150, format=None, **kwargs)` — saves via
  `Figure.savefig`.
- `positions -> list[tuple[float, float]]` — every world-space coordinate
  point in the current layout.

### Convenience functions

```python
quick_graph(sequence: str) -> FluxGraph
```
Builds a `Transcript`, calls `auto_stop_starts(reinitiation_prob=1)`, and
propagates flux through a `RiboSkeleton` built with `DEFAULT_BEHAVIOUR`.

```python
quick_plot(sequence: str, log_scale: float = 3) -> RiboGraphVis
```
Calls `quick_graph` and wraps the result in `RiboGraphVis`. Call
`.show()` on the result to display it.

---

## Extensibility

### Custom event types

Subclass `dynordg.core.Event` (or `Reading`, which additionally enforces
`0 < probability <= 1`) and override the `pipe` property to return a `Pipe`
describing the transition your event causes.

```python
from dynordg.core import Reading, Pipe, PhaseCondition, SubphaseCondition

class MyStopCodonReadthrough(Reading):
    """Stop-codon readthrough: stays translating instead of terminating."""

    @property
    def pipe(self) -> Pipe:
        return Pipe(
            output_phase=self.frame,
            probability=self.probability,
            input_position=self.position,
            phase_condition=PhaseCondition(self.frame),
            subphase_condition=SubphaseCondition(),
        )
```

Then add it with `transcript.add_event(MyStopCodonReadthrough(position, prob))`.

**Key building blocks**

- `Pipe(output_phase, probability, input_position, output_position=None, phase_condition=..., subphase_condition=..., remove_factors=None, add_factors=None)`
  — describes entry conditions and the resulting output state/position.
  `output_position` defaults to `input_position` (no positional movement).
  `remove_factors`/`add_factors` accept a string or tuple of strings;
  `remove_factors='all'` clears every subphase.
- `PhaseCondition(phases)` — wraps an `int`, a `set[int]`, or a
  `Callable[[int], bool]` and exposes `.match(phase)`. Prebuilt constants:
  `ANY_TRANSLATING` (`{1,2,3}`), `SCANNING` (`0`), `IN_SOLUTION` (`-1`).
- `SubphaseCondition(required=(), excluded=(), minimum={}, ceiling={})` —
  matches a `State`'s subphases: all of `required` must be present, none of
  `excluded` may be present, `minimum`/`ceiling` bound integer-valued
  subphases.

### Custom factor behaviours

Subclass the abstract `dynordg.core.FactorBehaviour` and implement:

- `applies(u: RiboNode, v: RiboNode) -> bool` — whether this behaviour fires
  for the `u → v` step.
- `fraction(u, v) -> float` — the flux fraction that undergoes the change.
- `transitions(u, v, weight) -> Transition` — the resulting transition.
- `priority: int` (property) — lower runs first; the built-in convention is
  `0` for structural decay, `10`–`20` for factor gain/loss.

Pass instances of your subclass in the `behaviours` list to `RiboSkeleton`.

### Customizing rendering

`RiboGraphVis` delegates layout math to a `LayoutEngine` and patch drawing
to a `RiboRenderer` (both in `dynordg.render`), and you can override either
independently:

- **Colours/styles** — subclass `RiboRenderer` and override the
  `COLOR_DICT` (maps edge-type strings — `'0'`–`'3'`, `'initiation'`,
  `'40s_retention'`, `'drop'`, `'load'`, `'frameshift'` — to matplotlib
  colours) and/or `STYLE_OVERRIDES` (per-type `alpha`/`linewidth`/`zorder`
  overrides), or override `edge_style(geom) -> EdgeStyle` entirely for full
  control.
- **Patch shapes** — subclass `EdgePainter` and override its `primitives()`
  method (or the individual `_body_primitives`/`_decay_primitives`/
  `_taper_primitives`/`_bulk_arrow_primitives` builders).
- **Layout math** — subclass `LayoutEngine` and override individual phase
  methods (`classify_edges`, `order_nodes`, `compute_geometries`,
  `align_layout`) as needed; each is documented in-line in `render.py` with
  its role in the four-phase pipeline (classify → order → local geometry →
  global alignment).

Pass your custom instances in:

```python
from dynordg import RiboGraphVis
from dynordg.render import LayoutEngine, RiboRenderer

vis = RiboGraphVis(flux_graph, engine=MyLayoutEngine(), renderer=MyRenderer())
```

---

## Appendix: phase & frame reference

| Phase value | Meaning |
|---|---|
| `-1` | In solution / bulk pool (off the mRNA) |
| `0`  | Scanning (40S subunit) |
| `1`, `2`, `3` | Translating (80S), in the frame where the ribosome's first subcodon position aligns with the 1st, 2nd, or 3rd nucleotide of the transcript, respectively |

`Event.frame` maps a 0-indexed nucleotide `position % 3` to a phase via
`FRAME_DICT = {1: 1, 2: 2, 0: 3}`.

Common subphase keys used by the built-in event/behaviour set:
`ternary_complex` (eIF2·GTP·Met-tRNAi, required to initiate),
`scanning_factors` (required to continue scanning or to retain after
termination), and `retained` (integer count of prior reinitiation events,
used with `Retention(..., limit=...)`).
