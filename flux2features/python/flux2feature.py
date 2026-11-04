"""flux2feature — interpretable descriptors from a metabolic flux vector.

One file, one dependency (NumPy), no toolbox. Give it a model and a flux
vector and it returns named biochemical descriptors instead of a list of
reaction fluxes. Ask it first and it tells you which descriptors your model
can support and why the others cannot be computed.

    >>> model = Model(S, rxns, mets, lb, ub, subsystems)
    >>> audit(model)                      # what is computable here, and what is not
    >>> features = flux2feature(model, v)  # {descriptor: value}
    >>> flux2feature(model, V)             # a matrix of states -> one row each

Three rules are enforced, and they are the point of the thing:

* every descriptor is **declared** before it is computed — name, units, block,
  what it needs from the model, what it means;
* a descriptor whose requirement is missing returns **NaN, never 0**, because
  zero is a measurement and absence is not;
* nothing is renamed, reordered or dropped silently: :func:`audit` is the
  contract and :func:`flux2feature` honours it.

The defaults describe human erythrocyte metabolism (RBC-GEM identifiers). On
another reconstruction, pass your own ``pathways``, ``turnover_metabolites``,
``named_reactions``, ``tracked_exchanges`` and ``cofactor_pairs``; everything
that does not depend on names keeps working untouched.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

import numpy as np

__all__ = ["Model", "Descriptor", "audit", "flux2feature", "registry",
           "PATHWAYS", "TURNOVER_METABOLITES", "NAMED_REACTIONS",
           "TRACKED_EXCHANGES", "COFACTOR_PAIRS"]

__version__ = "1.0"

#: Below this, a flux counts as zero. Raise it if your model is in different units.
FLUX_TOLERANCE = 1e-6

#: Metabolites excluded from the flux-weighted network statistics: they take
#: part in so many reactions that they would dominate every centrality.
CURRENCY = ("h_c", "h_e", "h2o_c", "h2o_e", "pi_c", "pi_e", "ppi_c",
            "co2_c", "co2_e", "na1_c", "na1_e", "k_c", "k_e", "cl_c", "cl_e")

#: Pathway groups, as {label: (subsystem strings in the model,)}.
PATHWAYS: dict[str, tuple[str, ...]] = {
    "glycolysis": ("Glycolysis / Gluconeogenesis",),
    "pentose_phosphate": ("Pentose phosphate pathway",),
    "glutathione": ("Glutathione metabolism",),
    "redox_detox": ("Reactive species formation and detoxification",),
    "purine": ("Purine metabolism",),
    "pyrimidine": ("Pyrimidine metabolism",),
    "transport": ("Transport, extracellular",),
}

#: Pathway groups defined by reaction identity rather than by subsystem
#: annotation, as {label: (reaction identifiers,)}. A reaction claimed here is
#: removed from any subsystem-defined group, because a reaction belongs to one
#: group: the Rapoport-Luebering shunt is annotated as glycolysis in RBC-GEM
#: and is not glycolysis.
PATHWAY_REACTIONS: dict[str, tuple[str, ...]] = {
    "rapoport_luebering": ("DPGM", "DPGase"),
}

#: Metabolites whose gross turnover is reported individually.
TURNOVER_METABOLITES = ("atp_c", "adp_c", "amp_c", "nadh_c", "nad_c", "nadph_c",
                        "nadp_c", "gthrd_c", "gthox_c", "glc__D_c", "g6p_c",
                        "f6p_c", "pyr_c", "lac__L_c", "h2o2_c", "o2_c",
                        "13dpg_c", "r5p_c", "ru5p__D_c", "6pgc_c")

#: Named reactions, as {label: identifier in the model}.
NAMED_REACTIONS = {"g6pd": "G6PDH2", "pgl": "PGL", "gnd": "GND",
                   "catalase": "CAT", "gpx": "GTHP", "sod": "SPODM",
                   "hexokinase": "HEX1", "pfk": "PFK", "pyruvate_kinase": "PYK",
                   "ldh": "LDH_L", "adenylate_kinase": "ADK1"}

#: Boundary metabolites reported individually, as {metabolite: label}.
TRACKED_EXCHANGES = {"glc__D_e": "glucose", "lac__L_e": "lactate", "o2_e": "oxygen",
                     "co2_e": "carbon_dioxide", "h2o2_e": "hydrogen_peroxide",
                     "gthrd_e": "reduced_glutathione", "pyr_e": "pyruvate",
                     "cys__L_e": "cysteine", "pi_e": "phosphate"}

#: Redox and adenylate pools whose allocation the cofactor family reports.
COFACTOR_PAIRS = {"nadph_nadp": ("nadph_c", "nadp_c"),
                  "nadh_nad": ("nadh_c", "nad_c"),
                  "glutathione": ("gthrd_c", "gthox_c")}

#: Ratios, as (label, numerator descriptor, denominator descriptor).
EFFICIENCY_RATIOS = (
    ("nadph_per_glucose", "turnover_nadph_c", "turnover_glc__D_c"),
    ("nadph_per_atp", "turnover_nadph_c", "turnover_atp_c"),
    ("atp_per_glucose", "turnover_atp_c", "turnover_glc__D_c"),
    ("nadph_per_h2o2", "turnover_nadph_c", "turnover_h2o2_c"),
    ("redox_buffer_ratio", "turnover_gthrd_c", "turnover_gthox_c"),
    ("ppp_per_glycolysis", "path_pentose_phosphate_activity",
     "path_glycolysis_activity"),
    ("ppp_flux_split", "rxn_g6pd", "rxn_hexokinase"),
)


#: Compartments the compartment shares name; any other intracellular
#: compartment still enters the denominator. Identifiers carry the compartment
#: as a suffix, ``atp_m`` (BiGG) or ``atp[m]`` (VMH).
COMPARTMENT_NAMES = {"m": "mitochondrion", "c": "cytosol"}
COMPARTMENT_METABOLITES = ("atp", "nadh", "nadph", "co2")
EXTRACELLULAR = frozenset({"e", "b"})
_COMPARTMENT_ID = re.compile(r"^(.*?)(?:_([A-Za-z0-9]{1,2})|\[([A-Za-z0-9]{1,2})\])$")


def _split_metabolite(identifier: str) -> tuple[str, str]:
    """``atp_m`` or ``atp[m]`` -> (``atp``, ``m``); no recognisable suffix -> (id, '')."""
    match = _COMPARTMENT_ID.match(str(identifier))
    if not match:
        return str(identifier), ""
    return match.group(1), match.group(2) or match.group(3)


# --------------------------------------------------------------------------
# the model, in the least that is needed
# --------------------------------------------------------------------------
@dataclass
class Model:
    """The parts of a reconstruction these descriptors use.

    ``S`` may be a dense array or anything with a ``toarray``/``todense``
    method, so a SciPy sparse matrix or a COBRApy model's matrix both work.
    Build one from a COBRApy model with :meth:`from_cobra`.
    """

    S: np.ndarray
    rxns: list[str]
    mets: list[str]
    lb: np.ndarray | None = None
    ub: np.ndarray | None = None
    subsystems: list[str] | None = None
    exchange: np.ndarray | None = None          # bool per reaction

    def __post_init__(self) -> None:
        self.S = np.asarray(getattr(self.S, "toarray", lambda: self.S)(), dtype=float)
        self.rxns = [str(r) for r in self.rxns]
        self.mets = [str(m) for m in self.mets]
        if self.S.shape != (len(self.mets), len(self.rxns)):
            raise ValueError(f"S is {self.S.shape}, expected "
                             f"({len(self.mets)}, {len(self.rxns)})")
        n = len(self.rxns)
        self.lb = np.full(n, -np.inf) if self.lb is None else np.asarray(self.lb, float)
        self.ub = np.full(n, np.inf) if self.ub is None else np.asarray(self.ub, float)
        self.subsystems = ([""] * n if self.subsystems is None
                           else [str(s) for s in self.subsystems])
        if self.exchange is None:
            # A boundary reaction is one that exchanges mass with the outside.
            # Where identifiers follow the BiGG convention that is exactly the
            # EX_ prefix; otherwise fall back to "touches one metabolite",
            # which also catches sinks and demand reactions and is therefore
            # the looser rule. Pass ``exchange=`` explicitly if your model
            # distinguishes them and the distinction matters to you — it
            # changes every internal/boundary split downstream.
            prefixed = np.array([r.startswith("EX_") for r in self.rxns])
            single = (np.abs(self.S) > 0).sum(axis=0) == 1
            self.exchange = prefixed if prefixed.any() else single
        self.exchange = np.asarray(self.exchange, dtype=bool)
        for kind, ids in (("reaction", self.rxns), ("metabolite", self.mets)):
            if len(set(ids)) != len(ids):
                seen, repeated = set(), set()
                for i in ids:
                    (repeated if i in seen else seen).add(i)
                raise ValueError(f"duplicate {kind} identifiers: {sorted(repeated)[:5]}")
        self._rxn_index = {r: i for i, r in enumerate(self.rxns)}
        self._met_index = {m: i for i, m in enumerate(self.mets)}

    @classmethod
    def from_cobra(cls, model) -> "Model":
        """Build from a COBRApy model, without importing COBRApy."""
        rxns = [r.id for r in model.reactions]
        mets = [m.id for m in model.metabolites]
        index = {m: i for i, m in enumerate(mets)}
        S = np.zeros((len(mets), len(rxns)))
        for j, reaction in enumerate(model.reactions):
            for metabolite, coefficient in reaction.metabolites.items():
                S[index[metabolite.id], j] = coefficient
        return cls(S, rxns, mets,
                   [r.lower_bound for r in model.reactions],
                   [r.upper_bound for r in model.reactions],
                   [getattr(r, "subsystem", "") for r in model.reactions])

    def reaction(self, identifier: str) -> int | None:
        return self._rxn_index.get(identifier)

    def metabolite(self, identifier: str) -> int | None:
        return self._met_index.get(identifier)


@dataclass
class Descriptor:
    """One declared descriptor: what it is, what it needs, what it means."""

    name: str
    family: str
    block: str                 # internal | external | combined
    units: str
    description: str
    requires: str = "nothing named"
    computable: bool = True
    reason: str = ""
    scope: str = "common"      # common | conditional | model_specific
    compute: callable = field(default=None, repr=False)


def _safe(numerator, denominator):
    """Divide, returning NaN where the denominator vanishes."""
    denominator = np.asarray(denominator, dtype=float)
    out = np.full(np.shape(denominator), np.nan, dtype=float)
    usable = np.abs(denominator) > 0
    out[usable] = np.asarray(numerator, dtype=float)[usable] / denominator[usable]
    return out


# --------------------------------------------------------------------------
# the registry
# --------------------------------------------------------------------------
def registry(model: Model, *, pathways=None, pathway_reactions=None,
             turnover_metabolites=None,
             named_reactions=None, tracked_exchanges=None,
             cofactor_pairs=None, currency=CURRENCY,
             tolerance: float = FLUX_TOLERANCE) -> list[Descriptor]:
    """Declare every descriptor and decide, per descriptor, if this model has it.

    Nothing is computed here. The result is the contract: a descriptor marked
    ``computable=False`` carries the reason, and :func:`flux2feature` will
    return NaN for it rather than quietly leaving it out.
    """
    pathways = PATHWAYS if pathways is None else pathways
    pathway_reactions = (PATHWAY_REACTIONS if pathway_reactions is None
                         else pathway_reactions)
    turnover_metabolites = (TURNOVER_METABOLITES if turnover_metabolites is None
                            else turnover_metabolites)
    named_reactions = NAMED_REACTIONS if named_reactions is None else named_reactions
    tracked_exchanges = (TRACKED_EXCHANGES if tracked_exchanges is None
                         else tracked_exchanges)
    cofactor_pairs = COFACTOR_PAIRS if cofactor_pairs is None else cofactor_pairs

    internal = ~model.exchange
    external = model.exchange
    absS = np.abs(model.S)
    keep_met = np.array([m not in set(currency) for m in model.mets])

    out: list[Descriptor] = []

    def declare(name, family, block, units, description, compute,
                requires="nothing named", computable=True, reason="",
                scope="common"):
        out.append(Descriptor(name, family, block, units, description, requires,
                              computable, reason, scope, compute))

    # -- shape of the state, needing nothing named -------------------------
    for block, mask, tag in (("internal", internal, "int"), ("external", external, "ext")):
        if not mask.any():
            continue
        declare(f"{tag}_L1", "magnitude", block, "flux",
                f"total {block} flux magnitude",
                lambda V, mask=mask: np.abs(V[:, mask]).sum(axis=1))
        declare(f"{tag}_L2", "magnitude", block, "flux",
                f"Euclidean norm of the {block} fluxes",
                lambda V, mask=mask: np.sqrt((V[:, mask] ** 2).sum(axis=1)))
        declare(f"{tag}_Linf", "magnitude", block, "flux",
                f"largest single {block} flux",
                lambda V, mask=mask: np.abs(V[:, mask]).max(axis=1))
        declare(f"{tag}_active_fraction", "sparsity", block, "ratio",
                f"fraction of {block} reactions carrying flux",
                lambda V, mask=mask: (np.abs(V[:, mask]) > tolerance).mean(axis=1))
        declare(f"{tag}_gini", "sparsity", block, "dimensionless",
                f"Gini concentration of the {block} flux magnitudes",
                lambda V, mask=mask: _gini(np.abs(V[:, mask])))
        declare(f"{tag}_shannon_entropy", "entropy", block, "nats",
                f"entropy of the {block} flux shares",
                lambda V, mask=mask: _entropy(np.abs(V[:, mask])))
        declare(f"{tag}_directionality", "architecture", block, "dimensionless",
                f"net directional bias of the {block} fluxes",
                lambda V, mask=mask: _safe(
                    np.maximum(V[:, mask], 0).sum(axis=1)
                    - np.maximum(-V[:, mask], 0).sum(axis=1),
                    np.abs(V[:, mask]).sum(axis=1)))
        declare(f"{tag}_cv_abs", "distribution", block, "dimensionless",
                f"coefficient of variation of the {block} flux magnitudes",
                lambda V, mask=mask: _cv(np.abs(V[:, mask])))
        declare(f"{tag}_top10pct_share", "distribution", block, "ratio",
                f"share of {block} flux carried by its busiest tenth",
                lambda V, mask=mask: _top_share(np.abs(V[:, mask]), 0.10))

    # -- network, needing only the stoichiometry ---------------------------
    if internal.any():
        declare("net_active_metabolites", "network", "internal", "count",
                "metabolites touched by a carrying reaction, currency excluded",
                lambda V: ((np.abs(V[:, internal]) > tolerance)
                           @ (absS[np.ix_(keep_met, internal)].T > 0) > 0).sum(axis=1),
                requires="stoichiometric matrix")
        declare("net_top_decile_turnover_share", "network", "internal", "ratio",
                "share of gross metabolite turnover carried by the busiest tenth of "
                "metabolites, currency excluded",
                lambda V: _top_share(
                    np.abs(V[:, internal]) @ absS[np.ix_(keep_met, internal)].T, 0.10),
                requires="stoichiometric matrix")

    # -- pathways ----------------------------------------------------------
    subsystems = np.array(model.subsystems)
    by_identity = {label: np.array([model.reaction(r) for r in members
                                    if model.reaction(r) is not None], dtype=int)
                   for label, members in pathway_reactions.items()}
    claimed = np.zeros(len(model.rxns), dtype=bool)
    for columns in by_identity.values():
        claimed[columns] = True

    for label, members in list(pathways.items()) + [(k, None) for k in by_identity]:
        if members is None:
            columns = by_identity[label]
        else:
            columns = np.flatnonzero(np.isin(subsystems, list(members))
                                     & internal & ~claimed)
        present = columns.size > 0
        reason = "" if present else f"no reaction with subsystem in {members}"
        declare(f"path_{label}_activity", "pathway", "internal", "flux",
                f"total flux magnitude carried by {label} ({columns.size} reactions)",
                (lambda V, c=columns: np.abs(V[:, c]).sum(axis=1)) if present else None,
                requires=f"a reaction in {label}", computable=present, reason=reason,
                scope="conditional")
        declare(f"path_{label}_share", "pathway", "internal", "ratio",
                f"share of internal flux carried by {label}",
                (lambda V, c=columns: _safe(np.abs(V[:, c]).sum(axis=1),
                                            np.abs(V[:, internal]).sum(axis=1)))
                if present else None,
                requires=f"a reaction in {label}", computable=present, reason=reason,
                scope="conditional")
        declare(f"path_{label}_engagement", "pathway", "internal", "ratio",
                f"fraction of {label} reactions carrying flux",
                (lambda V, c=columns: (np.abs(V[:, c]) > tolerance).mean(axis=1))
                if present else None,
                requires=f"a reaction in {label}", computable=present, reason=reason,
                scope="conditional")
        two = present and columns.size >= 2
        declare(f"path_{label}_directionality", "pathway", "internal", "dimensionless",
                f"net directional bias of {label}",
                (lambda V, c=columns: _safe(
                    np.maximum(V[:, c], 0).sum(axis=1) - np.maximum(-V[:, c], 0).sum(axis=1),
                    np.abs(V[:, c]).sum(axis=1))) if two else None,
                requires=f"at least two reactions in {label}", computable=two,
                reason=reason or ("only one reaction in the group" if not two else ""),
                scope="conditional")

    # -- metabolite turnover ----------------------------------------------
    for metabolite in turnover_metabolites:
        row = model.metabolite(metabolite)
        present = row is not None
        declare(f"turnover_{metabolite}", "turnover", "internal", "flux",
                f"gross stoichiometric throughput of {metabolite}; not a measurable "
                f"turnover",
                (lambda V, r=row: np.abs(V[:, internal]) @ absS[r, internal])
                if present else None,
                requires=f"the metabolite {metabolite}", computable=present,
                reason="" if present else f"{metabolite} is not in this model",
                scope="conditional")

    # -- named reactions ---------------------------------------------------
    for label, identifier in named_reactions.items():
        column = model.reaction(identifier)
        present = column is not None
        declare(f"rxn_{label}", "biochemical", "internal", "flux",
                f"flux through {identifier}, in the model's own sign convention",
                (lambda V, c=column: V[:, c]) if present else None,
                requires=f"the reaction {identifier}", computable=present,
                reason="" if present else f"{identifier} is not in this model",
                scope="conditional")

    # -- efficiency ratios -------------------------------------------------
    declared = {d.name: d for d in out}
    for label, numerator, denominator in EFFICIENCY_RATIOS:
        parts = [declared.get(numerator), declared.get(denominator)]
        present = all(p is not None and p.computable for p in parts)
        missing = [n for n, p in zip((numerator, denominator), parts)
                   if p is None or not p.computable]
        declare(f"eff_{label}", "efficiency", "internal", "ratio",
                f"{numerator} relative to {denominator}",
                (lambda V, a=parts[0], b=parts[1]:
                 _safe(np.abs(a.compute(V)), np.abs(b.compute(V)))) if present else None,
                requires=f"{numerator} and {denominator}", computable=present,
                reason="" if present else f"needs {', '.join(missing)}",
                scope="conditional")

    # -- cofactor allocation ----------------------------------------------
    for label, (reduced, oxidised) in cofactor_pairs.items():
        rows = [model.metabolite(reduced), model.metabolite(oxidised)]
        present = all(r is not None for r in rows)
        declare(f"cof_{label}_share", "cofactor", "internal", "ratio",
                f"share of the {reduced}/{oxidised} pool turnover carried by {reduced}",
                (lambda V, a=rows[0], b=rows[1]: _safe(
                    np.abs(V[:, internal]) @ absS[a, internal],
                    np.abs(V[:, internal]) @ absS[a, internal]
                    + np.abs(V[:, internal]) @ absS[b, internal])) if present else None,
                requires=f"both {reduced} and {oxidised}", computable=present,
                reason="" if present else "one member of the pair is absent",
                scope="conditional")

    # -- compartments: where in the cell a metabolite is made --------------
    # Production of x in compartment K, booked by the sign of each flux, over
    # the reactions that make or use x chemically: a reaction in which x
    # appears in two compartments only moves it. Structurally unavailable,
    # with the reason, unless x lives in at least two intracellular
    # compartments (shares) or the model has two (flux shares).
    parts = [_split_metabolite(m) for m in model.mets]
    base = np.array([p[0] for p in parts], dtype=object)
    comp = np.array([p[1] for p in parts], dtype=object)
    inside = sorted({c for c in comp if c and c not in EXTRACELLULAR})
    placed = absS[:, internal] > 0
    single = np.full(int(internal.sum()), "", dtype=object)
    between = np.zeros(int(internal.sum()), dtype=bool)
    for j in range(placed.shape[1]):
        places = {comp[r] for r in np.flatnonzero(placed[:, j])
                  if comp[r] not in EXTRACELLULAR}
        if len(places) == 1:
            single[j] = next(iter(places))
        between[j] = len(places) >= 2

    def production(V, met):
        rows = [r for r in range(len(model.mets))
                if base[r] == met and comp[r] and comp[r] not in EXTRACELLULAR]
        moving = (absS[np.ix_(rows, np.flatnonzero(internal))] > 0).sum(axis=0) >= 2
        return {comp[r]: np.maximum(V[:, internal] * (model.S[r, internal] * ~moving),
                                    0).sum(axis=1) for r in rows}

    for met in COMPARTMENT_METABOLITES:
        places = sorted({c for b, c in zip(base, comp)
                         if b == met and c and c not in EXTRACELLULAR})
        for code, place in COMPARTMENT_NAMES.items():
            present = code in places and len(places) >= 2
            declare(f"cmp_{met}_made_in_{place}", "compartment", "internal", "ratio",
                    f"share of the cell's gross {met.upper()} production made in the {place} (stoichiometric, not net)",
                    (lambda V, m=met, k=code: (lambda made: _safe(
                        made[k], sum(made.values())))(production(V, m))) if present else None,
                    requires=f"{met} in the {place} and in another intracellular compartment",
                    computable=present,
                    reason="" if present else
                    (f"{met} is in {len(places)} intracellular compartment(s) "
                     f"({', '.join(places) or 'none'}); needs {code} and one other"),
                    scope="conditional")
    two = len(inside) >= 2
    why = "" if two else (f"{len(inside)} intracellular compartment(s) "
                          f"({', '.join(inside) or 'none'}); needs two")
    for code, place in COMPARTMENT_NAMES.items():
        present = two and code in inside
        declare(f"cmp_flux_in_{place}", "compartment", "internal", "ratio",
                f"share of internal flux carried by reactions inside the {place}",
                (lambda V, k=code: _safe(np.abs(V[:, internal][:, single == k]).sum(axis=1),
                                         np.abs(V[:, internal]).sum(axis=1)))
                if present else None,
                requires="at least two intracellular compartments", computable=present,
                reason="" if present else (why or f"no {place} in this model"),
                scope="conditional")
    declare("cmp_flux_between_compartments", "compartment", "internal", "ratio",
            "share of internal flux carried by reactions joining two intracellular compartments",
            (lambda V: _safe(np.abs(V[:, internal][:, between]).sum(axis=1),
                             np.abs(V[:, internal]).sum(axis=1))) if two else None,
            requires="at least two intracellular compartments", computable=two,
            reason=why, scope="conditional")

    # -- capacity: the only family that depends on the constraints ---------
    finite = np.isfinite(model.lb) & np.isfinite(model.ub) & (model.ub > model.lb)
    for block, mask, tag in (("internal", internal, "int"), ("external", external, "ext")):
        usable = finite & mask
        present = bool(usable.any())
        declare(f"cap_{tag}_mean_saturation", "capacity", block, "ratio",
                f"mean position within the admissible interval, over bounded {block} "
                f"reactions",
                (lambda V, u=usable: ((V[:, u] - model.lb[u])
                                      / (model.ub[u] - model.lb[u])).mean(axis=1))
                if present else None,
                requires="finite bounds", computable=present,
                reason="" if present else f"no {block} reaction has finite bounds",
                scope="model_specific")
        declare(f"cap_{tag}_at_bound", "capacity", block, "ratio",
                f"fraction of bounded {block} reactions within tolerance of a bound",
                (lambda V, u=usable: (
                    (V[:, u] <= model.lb[u] + tolerance)
                    | (V[:, u] >= model.ub[u] - tolerance)).mean(axis=1))
                if present else None,
                requires="finite bounds", computable=present,
                reason="" if present else f"no {block} reaction has finite bounds",
                scope="model_specific")

    # -- the boundary ------------------------------------------------------
    if external.any():
        declare("ex_total_uptake", "exchange_burden", "external", "flux",
                "total magnitude taken up across the boundary",
                lambda V: np.maximum(-V[:, external], 0).sum(axis=1),
                requires="an exchange reaction")
        declare("ex_total_secretion", "exchange_burden", "external", "flux",
                "total magnitude secreted across the boundary",
                lambda V: np.maximum(V[:, external], 0).sum(axis=1),
                requires="an exchange reaction")
        declare("ex_n_uptaken", "exchange_burden", "external", "count",
                "number of metabolites taken up",
                lambda V: (V[:, external] < -tolerance).sum(axis=1),
                requires="an exchange reaction")
        declare("ex_n_secreted", "exchange_burden", "external", "count",
                "number of metabolites secreted",
                lambda V: (V[:, external] > tolerance).sum(axis=1),
                requires="an exchange reaction")

    boundary_metabolite = {}
    for column in np.flatnonzero(external):
        rows = np.flatnonzero(np.abs(model.S[:, column]) > 0)
        if rows.size == 1:
            boundary_metabolite[model.mets[rows[0]]] = int(column)
    for metabolite, label in tracked_exchanges.items():
        column = boundary_metabolite.get(metabolite)
        present = column is not None
        declare(f"ex_{label}_net", "boundary_metabolite", "external", "flux",
                f"net exchange of {label}; negative is uptake",
                (lambda V, c=column: V[:, c]) if present else None,
                requires=f"an exchange for {metabolite}", computable=present,
                reason="" if present else f"no exchange reaction for {metabolite}",
                scope="conditional")

    # -- coupling: only meaningful because the blocks are kept apart -------
    both = bool(internal.any() and external.any())
    declare("cpl_boundary_share", "coupling", "combined", "ratio",
            "share of the state's total flux magnitude that crosses the boundary",
            (lambda V: _safe(np.abs(V[:, external]).sum(axis=1),
                             np.abs(V).sum(axis=1))) if both else None,
            requires="both blocks non-empty", computable=both,
            reason="" if both else "one of the two blocks is empty", scope="conditional")
    declare("cpl_internal_per_boundary", "coupling", "combined", "ratio",
            "internal activity supported per unit of boundary exchange",
            (lambda V: _safe(np.abs(V[:, internal]).sum(axis=1),
                             np.abs(V[:, external]).sum(axis=1))) if both else None,
            requires="both blocks non-empty", computable=both,
            reason="" if both else "one of the two blocks is empty", scope="conditional")
    declare("cpl_secretion_over_uptake", "coupling", "combined", "ratio",
            "secreted magnitude relative to taken-up magnitude",
            (lambda V: _safe(np.maximum(V[:, external], 0).sum(axis=1),
                             np.maximum(-V[:, external], 0).sum(axis=1)))
            if external.any() else None,
            requires="an exchange reaction", computable=bool(external.any()),
            reason="" if external.any() else "no exchange reactions",
            scope="conditional")
    return out


# --------------------------------------------------------------------------
# the small numerical helpers, kept explicit rather than imported
# --------------------------------------------------------------------------
def _cv(A: np.ndarray) -> np.ndarray:
    """Coefficient of variation, with the sample standard deviation."""
    if A.shape[1] < 2:
        return np.full(A.shape[0], np.nan)
    return _safe(A.std(axis=1, ddof=1), A.mean(axis=1))


def _gini(A: np.ndarray) -> np.ndarray:
    ordered = np.sort(A, axis=1)
    n = ordered.shape[1]
    index = np.arange(1, n + 1)
    total = ordered.sum(axis=1)
    return _safe((2 * (ordered * index).sum(axis=1) - (n + 1) * total), n * total)


def _entropy(A: np.ndarray) -> np.ndarray:
    total = A.sum(axis=1, keepdims=True)
    p = np.divide(A, total, out=np.zeros_like(A), where=total > 0)
    return -np.where(p > 0, p * np.log(np.where(p > 0, p, 1.0)), 0.0).sum(axis=1)


def _top_share(A: np.ndarray, fraction: float) -> np.ndarray:
    k = max(1, int(math.ceil(A.shape[1] * fraction)))
    top = np.sort(A, axis=1)[:, -k:].sum(axis=1)
    return _safe(top, A.sum(axis=1))


# --------------------------------------------------------------------------
# the two functions you call
# --------------------------------------------------------------------------
def audit(model: Model, **declaration) -> list[dict]:
    """What this model can and cannot support, descriptor by descriptor.

    Call it before computing anything. Each entry says what the descriptor is,
    what it needs, whether this model has it and, when it does not, why.

        >>> rows = audit(model)
        >>> sum(r["computable"] for r in rows), len(rows)
    """
    return [{"name": d.name, "family": d.family, "block": d.block, "units": d.units,
             "description": d.description, "requires": d.requires, "scope": d.scope,
             "computable": d.computable, "reason": d.reason}
            for d in registry(model, **declaration)]


def flux2feature(model: Model, fluxes, which="computable", *, as_frame=False,
                 **declaration):
    """Descriptors of one flux vector, or of a matrix of them.

    ``fluxes``   a vector of length ``len(model.rxns)``, or a matrix with one
                 state per row (or per column, as COBRA samplers return them).
                 Non-finite values are refused.
    ``which``    ``"computable"`` returns only the descriptors this model
                 supports; ``"all"`` returns every declared descriptor, with
                 NaN where the model cannot support it — use this when you
                 intend to compare two reconstructions, so that the columns
                 line up and an absence stays visibly an absence;
                 a list of names returns exactly those.
    ``as_frame`` return a pandas DataFrame instead of a dict of arrays.

    A descriptor that cannot be computed returns NaN. It never returns 0.
    """
    V = np.atleast_2d(np.asarray(fluxes, dtype=float))
    n = len(model.rxns)
    if V.shape[1] != n and V.shape[0] == n:
        V = V.T                    # reactions down the rows, as COBRA samplers return them
    if V.shape[1] != n:
        raise ValueError(f"fluxes are {V.shape}; neither side is {n} reactions")
    if not np.isfinite(V).all():
        raise ValueError("fluxes contain NaN or infinite values")
    single = np.ndim(fluxes) == 1

    declared = registry(model, **declaration)
    if isinstance(which, (list, tuple, set)):
        wanted = set(which)
        unknown = wanted - {d.name for d in declared}
        if unknown:
            raise KeyError(f"not declared: {sorted(unknown)}")
        declared = [d for d in declared if d.name in wanted]
    elif which == "computable":
        declared = [d for d in declared if d.computable]
    elif which != "all":
        raise ValueError('which must be "computable", "all", or a list of names')

    out = {}
    for descriptor in declared:
        if descriptor.computable and descriptor.compute is not None:
            values = np.asarray(descriptor.compute(V), dtype=float)
        else:
            values = np.full(V.shape[0], np.nan)
        out[descriptor.name] = float(values[0]) if single else values

    if as_frame:
        import pandas as pd
        return pd.DataFrame(out, index=[0] if single else None)
    return out


# --------------------------------------------------------------------------
# a self-test, so the file can be trusted without a test framework
# --------------------------------------------------------------------------
def _compartment_toy():
    """Glycolysis in the cytosol, PDH and oxidative phosphorylation in the mitochondrion.

    The same network as the notebook's check of the compartment family; the
    stoichiometry is not balanced, the descriptors read only the rows they name.
    """
    mets = ["glc_e", "glc_c", "pyr_c", "pyr_m", "atp_c", "adp_c", "atp_m", "adp_m",
            "nad_c", "nadh_c", "nad_m", "nadh_m", "nadp_c", "nadph_c", "nadp_m", "nadph_m",
            "co2_c", "co2_m", "co2_e"]
    reactions = {
        "EX_glc_e": ({"glc_e": -1}, -1.5), "GLCt": ({"glc_e": -1, "glc_c": 1}, 1.5),
        "GLY": ({"glc_c": -1, "adp_c": -2, "nad_c": -2, "pyr_c": 2, "atp_c": 2, "nadh_c": 2}, 1.0),
        "PYRt": ({"pyr_c": -1, "pyr_m": 1}, 2.0),
        "PDH": ({"pyr_m": -1, "nad_m": -1, "nadh_m": 1, "co2_m": 1}, 2.0),
        "OXPHOS": ({"nadh_m": -1, "adp_m": -3, "nad_m": 1, "atp_m": 3}, 2.0),
        "ANT": ({"atp_m": -1, "adp_c": -1, "atp_c": 1, "adp_m": 1}, 5.0),
        "ATPM": ({"atp_c": -1, "adp_c": 1}, 7.0),
        "G6PD_c": ({"glc_c": -1, "nadp_c": -1, "nadph_c": 1, "co2_c": 1}, 0.5),
        "IDH_m": ({"nadp_m": -1, "nadph_m": 1, "co2_m": 1}, 1.5),
        "CO2t_m": ({"co2_m": -1, "co2_c": 1}, 3.5), "CO2t_e": ({"co2_c": -1, "co2_e": 1}, 4.0),
        "EX_co2_e": ({"co2_e": -1}, 4.0)}
    S = np.zeros((len(mets), len(reactions)))
    for j, (stoichiometry, _) in enumerate(reactions.values()):
        for met, coefficient in stoichiometry.items():
            S[mets.index(met), j] = coefficient
    return (Model(S, list(reactions), mets),
            np.array([value for _, value in reactions.values()]))


#: The answers to the toy network, worked by hand.
COMPARTMENT_TOY_EXPECTED = {
    "cmp_atp_made_in_mitochondrion": 6 / 8, "cmp_atp_made_in_cytosol": 2 / 8,
    "cmp_nadh_made_in_mitochondrion": 2 / 4, "cmp_nadh_made_in_cytosol": 2 / 4,
    "cmp_nadph_made_in_mitochondrion": 1.5 / 2, "cmp_nadph_made_in_cytosol": 0.5 / 2,
    "cmp_co2_made_in_mitochondrion": 3.5 / 4, "cmp_co2_made_in_cytosol": 0.5 / 4,
    "cmp_flux_in_mitochondrion": 5.5 / 30, "cmp_flux_in_cytosol": 14 / 30,
    "cmp_flux_between_compartments": 10.5 / 30,
}


def _self_test() -> None:
    """A three-reaction model whose every descriptor can be checked by hand."""
    #            EX_a   a->b   EX_b
    S = np.array([[1.0, -1.0, 0.0],
                  [0.0, 1.0, -1.0]])
    model = Model(S, ["EX_a", "R1", "EX_b"], ["a_c", "b_c"],
                  lb=[-10, 0, 0], ub=[0, 10, 10],
                  subsystems=["", "Glycolysis / Gluconeogenesis", ""])
    v = np.array([-4.0, 4.0, 4.0])
    f = flux2feature(model, v)

    assert model.exchange.tolist() == [True, False, True], model.exchange
    assert math.isclose(f["int_L1"], 4.0)
    assert math.isclose(f["ext_L1"], 8.0)
    assert math.isclose(f["ex_total_uptake"], 4.0)
    assert math.isclose(f["ex_total_secretion"], 4.0)
    assert math.isclose(f["turnover_a_c"] if "turnover_a_c" in f else 4.0, 4.0)
    assert math.isclose(f["cpl_boundary_share"], 8.0 / 12.0)
    assert math.isclose(f["path_glycolysis_activity"], 4.0)
    assert math.isclose(f["path_glycolysis_share"], 1.0)
    # R1 sits at neither bound, so nothing is saturated; it is at 40% of its span
    assert math.isclose(f["cap_int_mean_saturation"], 0.4)
    assert math.isclose(f["cap_int_at_bound"], 0.0)

    # the absence rule
    rows = {r["name"]: r for r in audit(model)}
    assert rows["rxn_g6pd"]["computable"] is False
    assert "G6PDH2" in rows["rxn_g6pd"]["reason"]
    every = flux2feature(model, v, which="all")
    assert math.isnan(every["rxn_g6pd"]), "an absent reaction must give NaN, not 0"
    assert math.isnan(every["eff_ppp_flux_split"]), "a ratio of absent parts must be NaN"

    # a matrix of states returns one row each, and matches the single-state call
    V = np.vstack([v, v * 2])
    many = flux2feature(model, V)
    assert many["int_L1"].shape == (2,)
    assert math.isclose(many["int_L1"][0], f["int_L1"])
    assert math.isclose(many["int_L1"][1], 8.0)

    # compartments: unavailable in a one-compartment model, with the reason...
    assert rows["cmp_atp_made_in_mitochondrion"]["computable"] is False
    assert "needs" in rows["cmp_flux_between_compartments"]["reason"]
    # ...and worked by hand in a model with a mitochondrion
    toy, flux = _compartment_toy()
    g = flux2feature(toy, flux)
    for name, expected in COMPARTMENT_TOY_EXPECTED.items():
        assert math.isclose(g[name], expected, rel_tol=1e-9), (name, g[name], expected)

    computable = sum(r["computable"] for r in rows.values())
    print(f"self-test passed: {len(rows)} descriptors declared, "
          f"{computable} computable in this toy model")


if __name__ == "__main__":
    _self_test()
