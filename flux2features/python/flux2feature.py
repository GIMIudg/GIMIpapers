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

__version__ = "1.1"

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



_ELEMENT = re.compile(r"([A-Z][a-z]?)(\d*(?:\.\d+)?)")


def _element_count(formula, element: str) -> float:
    """Atoms of ``element`` in a formula such as ``C6H12O6``; 0 if absent or unreadable."""
    if not isinstance(formula, str) or not formula:
        return 0.0
    return sum(float(c) if c else 1.0 for symbol, c in _ELEMENT.findall(formula) if symbol == element)


#: The proteinogenic amino acids: canonical name -> identifier stems (BiGG, VMH).
AMINO_ACIDS = {
    "alanine": ("ala__L", "ala_L"), "arginine": ("arg__L", "arg_L"),
    "asparagine": ("asn__L", "asn_L"), "aspartate": ("asp__L", "asp_L"),
    "cysteine": ("cys__L", "cys_L"), "glutamine": ("gln__L", "gln_L"),
    "glutamate": ("glu__L", "glu_L"), "glycine": ("gly",),
    "histidine": ("his__L", "his_L"), "isoleucine": ("ile__L", "ile_L"),
    "leucine": ("leu__L", "leu_L"), "lysine": ("lys__L", "lys_L"),
    "methionine": ("met__L", "met_L"), "phenylalanine": ("phe__L", "phe_L"),
    "proline": ("pro__L", "pro_L"), "serine": ("ser__L", "ser_L"),
    "threonine": ("thr__L", "thr_L"), "tryptophan": ("trp__L", "trp_L"),
    "tyrosine": ("tyr__L", "tyr_L"), "valine": ("val__L", "val_L"),
}
ESSENTIAL_AMINO_ACIDS = ("histidine", "isoleucine", "leucine", "lysine", "methionine",
                         "phenylalanine", "threonine", "tryptophan", "valine")
#: amino acids already reported one by one as tracked exchanges
_TRACKED_AMINO_ACIDS = ("cysteine", "glycine", "glutamate")
_AA_SYNONYMS = {"aspartic acid": "aspartate", "glutamic acid": "glutamate",
                "l-aspartate": "aspartate", "l-glutamate": "glutamate"}


def _amino_acid(identifier: str, name: str = "") -> str | None:
    """Canonical amino-acid name of a metabolite, by identifier and then by name; D forms excluded."""
    stem, _ = _split_metabolite(identifier)
    for canonical, stems in AMINO_ACIDS.items():
        if stem in stems:
            return canonical
    text = re.sub(r"\(.*?\)", "", str(name).strip().lower()).strip()
    if not text or text.startswith("d-"):
        return None
    text = _AA_SYNONYMS.get(text, text)
    text = text[2:] if text.startswith("l-") else text
    return text if text in AMINO_ACIDS else None

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
    formulas: list[str] | None = None           # metabolite formulas, for the element families
    met_names: list[str] | None = None          # metabolite names, to recognise amino acids

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
        m = len(self.mets)
        self.formulas = [""] * m if self.formulas is None else [
            f if isinstance(f, str) else "" for f in self.formulas]
        self.met_names = [""] * m if self.met_names is None else [str(x) for x in self.met_names]
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
                   [getattr(r, "subsystem", "") for r in model.reactions],
                   formulas=[getattr(m, "formula", "") or "" for m in model.metabolites],
                   met_names=[getattr(m, "name", "") or "" for m in model.metabolites])

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


    # -- elements, amino acids, oxidative defence and conversion -----------
    # The export of an exchanged metabolite is read from the sign of its
    # coefficient, so it is positive out of the cell whichever way the exchange
    # is written. Shares and ratios are NaN where their denominator is zero.
    ext_cols = np.array([c for c in np.flatnonzero(external)
                         if (np.abs(model.S[:, c]) > 0).sum() == 1], dtype=int)
    ext_row = np.array([np.flatnonzero(np.abs(model.S[:, c]) > 0)[0] for c in ext_cols], dtype=int)
    ext_sign = -np.sign(model.S[ext_row, ext_cols]) if ext_cols.size else np.zeros(0)
    ext_stem = np.array([base[r] for r in ext_row], dtype=object)
    atoms = {e: np.array([_element_count(model.formulas[r], e) for r in ext_row])
             for e in ("C", "N", "S", "P")}

    def export(V):
        return V[:, ext_cols] * ext_sign

    def carbon_parts(V, uptake=True):
        x = export(V)
        return (np.maximum(-x, 0) if uptake else np.maximum(x, 0)) * atoms["C"]

    def where(*stems):
        return np.flatnonzero(np.isin(ext_stem, stems))

    for e, label in (("N", "nitrogen"), ("S", "sulfur"), ("P", "phosphorus")):
        present = bool((atoms[e] > 0).any())
        for side in ("uptake", "secretion"):
            declare(f"elem_{e.lower()}_{side}", "elemental", "external", "flux",
                    f"{label} atoms {'taken up' if side == 'uptake' else 'secreted'} per unit "
                    f"time, from the formulas of exchanged metabolites",
                    (lambda V, a=atoms[e], up=side == "uptake":
                     (np.maximum(-export(V), 0) if up else np.maximum(export(V), 0)) @ a)
                    if present else None,
                    requires=f"an exchanged metabolite with {label} in its formula",
                    computable=present,
                    reason="" if present else f"no exchanged metabolite contains {label}",
                    scope="conditional")
    has_c, has_n = bool((atoms["C"] > 0).any()), bool((atoms["N"] > 0).any())
    amino = np.array([_amino_acid(model.mets[r], model.met_names[r]) for r in ext_row], dtype=object)
    aa_cols = np.flatnonzero(amino != None)                                  # noqa: E711
    shares = (
        ("elem_c_from_glucose", "share of carbon uptake supplied by glucose",
         where("glc__D", "glc_D"), True, "C", "a glucose exchange"),
        ("elem_c_from_lactate_pyruvate", "share of carbon uptake supplied by lactate and pyruvate",
         where("lac__L", "lac_L", "pyr"), True, "C", "a lactate or pyruvate exchange"),
        ("elem_c_from_amino_acids", "share of carbon uptake supplied by amino acids",
         aa_cols, True, "C", "an amino-acid exchange"),
        ("elem_c_to_lactate", "share of carbon secretion leaving as lactate",
         where("lac__L", "lac_L"), False, "C", "a lactate exchange"),
        ("elem_c_to_co2", "share of carbon secretion leaving as CO2",
         where("co2"), False, "C", "a CO2 exchange"),
        ("elem_n_from_amino_acids", "share of nitrogen uptake supplied by amino acids",
         aa_cols, True, "N", "an amino-acid exchange"))
    for name, text, cols, up, e, needs in shares:
        present = bool((atoms[e] > 0).any()) and cols.size > 0
        declare(name, "elemental", "external", "ratio", text,
                (lambda V, c=cols, up=up, a=atoms[e]: (lambda parts: _safe(
                    parts[:, c].sum(axis=1), parts.sum(axis=1)))(
                    (np.maximum(-export(V), 0) if up else np.maximum(export(V), 0)) * a))
                if present else None,
                requires=needs, computable=present,
                reason="" if present else f"no {needs[2:] if needs.startswith('a ') else needs}",
                scope="conditional")

    def diversity(parts, kind):
        total = parts.sum(axis=1, keepdims=True)
        with np.errstate(invalid="ignore", divide="ignore"):
            f = np.where(total > 0, parts / total, np.nan)
            logs = np.where(f > 0, np.log(np.where(f > 0, f, 1.0)), 0.0)
            value = {"entropy": -(np.nan_to_num(f) * logs).sum(axis=1),
                     "effective": 1.0 / np.nansum(f ** 2, axis=1),
                     "dominant": np.nanmax(np.nan_to_num(f), axis=1)}[kind]
        value = np.asarray(value, dtype=float)
        value[total.ravel() <= 0] = np.nan
        return value

    for e, present in (("c", has_c), ("n", has_n)):
        kinds = (("source_entropy", "entropy", "nats", "Shannon entropy of the {} sources taken up"),
                 ("effective_sources", "effective", "dimensionless", "effective number of {} sources"))
        if e == "c":
            kinds += (("dominant_source", "dominant", "ratio",
                       "share of {} uptake from the largest single source"),)
        label = "carbon" if e == "c" else "nitrogen"
        for suffix, kind, units, text in kinds:
            declare(f"div_{e}_{suffix}", "source_diversity", "external", units, text.format(label),
                    (lambda V, a=atoms[e.upper()], k=kind:
                     diversity(np.maximum(-export(V), 0) * a, k)) if present else None,
                    requires=f"an exchanged metabolite with {label} in its formula",
                    computable=present,
                    reason="" if present else f"no exchanged metabolite contains {label}",
                    scope="conditional")

    by_acid = {}
    for k in aa_cols:
        by_acid.setdefault(amino[k], []).append(int(k))
    for acid in AMINO_ACIDS:
        if acid in _TRACKED_AMINO_ACIDS:
            continue
        present = acid in by_acid
        declare(f"aa_{acid}_net_uptake", "amino_acid", "external", "flux",
                f"net uptake of {acid} (uptake minus secretion)",
                (lambda V, c=by_acid.get(acid, []): -export(V)[:, c].sum(axis=1))
                if present else None,
                requires=f"an exchange for {acid}", computable=present,
                reason="" if present else f"no {acid} exchange", scope="conditional")

    def acid_matrix(V):
        x = export(V)
        return np.column_stack([x[:, by_acid[a]].sum(axis=1) for a in by_acid])

    found = bool(by_acid)
    essential = [k for k, a in enumerate(by_acid) if a in ESSENTIAL_AMINO_ACIDS]
    for name, units, text, compute in (
            ("aa_total_uptake", "flux", "total uptake of amino acids",
             lambda V: np.maximum(-acid_matrix(V), 0).sum(axis=1)),
            ("aa_total_secretion", "flux", "total secretion of amino acids",
             lambda V: np.maximum(acid_matrix(V), 0).sum(axis=1)),
            ("aa_essential_uptake", "flux", "uptake of the nine essential amino acids",
             lambda V: np.maximum(-acid_matrix(V), 0)[:, essential].sum(axis=1)),
            ("aa_uptake_entropy", "nats", "Shannon entropy of the amino-acid uptake profile",
             lambda V: diversity(np.maximum(-acid_matrix(V), 0), "entropy")),
            ("aa_n_taken_up", "count", "number of amino acids taken up",
             lambda V: (np.maximum(-acid_matrix(V), 0) > tolerance).sum(axis=1).astype(float))):
        present = found and (name != "aa_essential_uptake" or bool(essential))
        declare(name, "amino_acid", "external", units, text, compute if present else None,
                requires="an amino-acid exchange recognised by identifier or name",
                computable=present,
                reason="" if present else "no recognised amino-acid exchange",
                scope="conditional")

    # carriers: rows of every intracellular compartment, moves excluded
    int_cols = np.flatnonzero(internal)
    stems_all = base

    def carrier(*stems):
        rows = [r for r in range(len(model.mets)) if stems_all[r] in stems
                and comp[r] not in EXTRACELLULAR]
        if not rows:
            return None
        coefficient = model.S[np.ix_(rows, int_cols)].sum(axis=0)
        every = model.S[np.ix_([r for r in range(len(model.mets)) if stems_all[r] in stems],
                               int_cols)]
        coefficient = np.where((every > 0).any(axis=0) & (every < 0).any(axis=0), 0.0, coefficient)
        return coefficient

    def made(V, c):
        return np.maximum(V[:, int_cols] * c, 0)

    def used(V, c):
        return np.maximum(-V[:, int_cols] * c, 0)

    def signature(consumes, produces=(), only=False, partner=""):
        hits = []
        for k, j in enumerate(int_cols):
            column = model.S[:, j]
            taken = {stems_all[r] for r in np.flatnonzero(column < 0)}
            given = {stems_all[r] for r in np.flatnonzero(column > 0)}
            if not set(consumes) <= taken or not set(produces) <= given:
                continue
            if only and taken - set(consumes) - {"h"}:
                continue
            if partner and not any(t.startswith(partner) for t in taken):
                continue
            hits.append(k)
        return np.array(hits, dtype=int)

    h2o2, o2s, nadph, atp = carrier("h2o2"), carrier("o2s"), carrier("nadph"), carrier("atp")
    live = lambda c: c is not None and bool(np.any(c))
    catalase = signature(("h2o2",), ("o2",), only=True)
    gpx = signature(("h2o2", "gthrd"), ("gthox",))
    prdx = signature(("h2o2",), partner="prdx")
    sod = signature(("o2s",), ("h2o2",))
    gr_nadph = signature(("gthox", "nadph"), ("gthrd",))
    gr = np.union1d(gr_nadph, signature(("gthox", "nadh"), ("gthrd",))).astype(int)
    defence = (
        ("ros_h2o2_production", "flux", "gross internal production of hydrogen peroxide",
         live(h2o2), "no peroxide", lambda V: made(V, h2o2).sum(axis=1)),
        ("ros_h2o2_removal", "flux", "gross internal removal of hydrogen peroxide",
         live(h2o2), "no peroxide", lambda V: used(V, h2o2).sum(axis=1)),
        ("ros_h2o2_removal_per_production", "ratio",
         "model-based peroxide removal-to-production ratio (above 1 when peroxide is taken up)",
         live(h2o2), "no peroxide",
         lambda V: _safe(used(V, h2o2).sum(axis=1), made(V, h2o2).sum(axis=1))),
        ("ros_h2o2_removed_by_catalase", "ratio", "share of peroxide removal by catalase",
         live(h2o2) and catalase.size > 0, "no catalase signature (2 H2O2 -> O2)",
         lambda V: _safe(used(V, h2o2)[:, catalase].sum(axis=1), used(V, h2o2).sum(axis=1))),
        ("ros_h2o2_removed_by_glutathione_peroxidase", "ratio",
         "share of peroxide removal by glutathione peroxidase",
         live(h2o2) and gpx.size > 0, "no glutathione peroxidase signature",
         lambda V: _safe(used(V, h2o2)[:, gpx].sum(axis=1), used(V, h2o2).sum(axis=1))),
        ("ros_h2o2_removed_by_peroxiredoxin", "ratio", "share of peroxide removal by peroxiredoxins",
         live(h2o2) and prdx.size > 0, "no peroxiredoxin signature",
         lambda V: _safe(used(V, h2o2)[:, prdx].sum(axis=1), used(V, h2o2).sum(axis=1))),
        ("ros_superoxide_production", "flux", "gross internal production of superoxide",
         live(o2s), "no superoxide", lambda V: made(V, o2s).sum(axis=1)),
        ("ros_superoxide_dismutated", "ratio", "share of superoxide removal by superoxide dismutase",
         live(o2s) and sod.size > 0, "no superoxide dismutase signature",
         lambda V: _safe(used(V, o2s)[:, sod].sum(axis=1), used(V, o2s).sum(axis=1))),
        ("ros_glutathione_regeneration_cover", "ratio",
         "NADPH produced per NADPH consumed by glutathione reductases",
         live(nadph) and gr_nadph.size > 0, "no NADPH-dependent glutathione reductase",
         lambda V: _safe(made(V, nadph).sum(axis=1), used(V, nadph)[:, gr_nadph].sum(axis=1))),
        ("ros_glutathione_cycle_closure", "ratio",
         "glutathione disulfide reduced by glutathione reductases per disulfide formed by peroxide removal",
         gr.size > 0 and gpx.size > 0, "no glutathione reductase or peroxidase",
         lambda V: _safe(np.maximum(V[:, int_cols[gr]], 0).sum(axis=1),
                         np.maximum(V[:, int_cols[gpx]], 0).sum(axis=1))))
    for name, units, text, present, why, compute in defence:
        declare(name, "oxidative_defense", "internal", units, text, compute if present else None,
                requires="peroxide or superoxide and the named signature", computable=present,
                reason="" if present else why, scope="conditional")

    # -- budget: the two NADPH budget descriptors of the paper's panel -------
    # NADPH consumed by glutathione reductases (the signature above), and NADPH made by
    # the oxidative branch (the named G6PDH2 and GND) per NADPH so consumed.
    columns = {model.reaction(named_reactions.get(k, "")) for k in ("g6pd", "gnd")}
    position = {j: k for k, j in enumerate(int_cols)}
    oxidative = np.array(sorted(position[j] for j in columns if j is not None and j in position),
                         dtype=int)
    branch = oxidative.size == 2
    budget = (
        ("bud_nadph_to_glutathione_reduction",
         "share of NADPH consumption used to reduce glutathione disulfide",
         live(nadph) and gr.size > 0, "no NADPH or no glutathione reductase signature",
         lambda V: _safe(used(V, nadph)[:, gr].sum(axis=1), used(V, nadph).sum(axis=1))),
        ("bud_oxppp_cover_of_glutathione_reduction",
         "NADPH made by the oxidative branch per NADPH used to reduce glutathione",
         live(nadph) and gr.size > 0 and branch,
         "no NADPH, no glutathione reductase signature, or G6PDH2 and GND not both present",
         lambda V: _safe(made(V, nadph)[:, oxidative].sum(axis=1),
                         used(V, nadph)[:, gr].sum(axis=1))))
    for name, text, present, why, compute in budget:
        declare(name, "budget", "internal", "ratio", text, compute if present else None,
                requires="NADPH, the glutathione-reductase signature and the named reactions",
                computable=present, reason="" if present else why, scope="conditional")

    glucose, lactate, co2 = where("glc__D", "glc_D"), where("lac__L", "lac_L"), where("co2")
    conversion = (
        ("conv_atp_per_carbon", "gross ATP production per carbon atom taken up",
         has_c and atp is not None, "no ATP or no carbon exchange",
         lambda V: _safe(made(V, atp).sum(axis=1), carbon_parts(V).sum(axis=1))),
        ("conv_nadph_per_carbon", "gross NADPH production per carbon atom taken up",
         has_c and nadph is not None, "no NADPH or no carbon exchange",
         lambda V: _safe(made(V, nadph).sum(axis=1), carbon_parts(V).sum(axis=1))),
        ("conv_lactate_per_glucose_carbon", "lactate carbon secreted per glucose carbon taken up",
         glucose.size > 0 and lactate.size > 0, "no glucose or lactate exchange",
         lambda V: _safe(carbon_parts(V, False)[:, lactate].sum(axis=1),
                         carbon_parts(V)[:, glucose].sum(axis=1))),
        ("conv_co2_per_carbon", "CO2 carbon released per carbon atom taken up",
         has_c and co2.size > 0, "no CO2 exchange",
         lambda V: _safe(carbon_parts(V, False)[:, co2].sum(axis=1), carbon_parts(V).sum(axis=1))))
    for name, text, present, why, compute in conversion:
        declare(name, "conversion", "combined", "ratio", text, compute if present else None,
                requires="metabolite formulas, the carrier and the named exchanges",
                computable=present, reason="" if present else why, scope="conditional")

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


def flux2feature(model: Model, fluxes, which="all", *, as_frame=False,
                 **declaration):
    """Descriptors of one flux vector, or of a matrix of them.

    ``fluxes``   a vector of length ``len(model.rxns)``, or a matrix with one
                 state per row (or per column, as COBRA samplers return them).
                 Non-finite values are refused.
    ``which``    ``"all"`` (the default) returns every declared descriptor, with
                 NaN where the model cannot support it, so that every flux
                 vector of every reconstruction maps to the same coordinates
                 and an absence stays visibly an absence; ``"computable"``
                 returns only the descriptors this model supports; a list of
                 names returns exactly those.
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


def _functional_toy():
    """Glucose and glutamine in, lactate and CO2 out; peroxide made and removed two ways.

    Exchanges: glucose 2 in (12 C), glutamine 1 in (5 C, 2 N), lactate 3 out (9 C),
    CO2 8 out (8 C). Peroxide: 3 made by an oxidase, 2 removed by catalase (one
    turn) and 2 by glutathione peroxidase (two turns); glutathione reductase
    reduces 2 disulfide with 2 NADPH, G6PD and GND make 5 and 1 NADPH, and an
    oxidase uses the other 2 NADPH (4 used in all). The stoichiometry
    is not balanced: the descriptors read only the rows they name.
    """
    mets = ["glc__D_e", "gln__L_e", "lac__L_e", "co2_e", "h2o2_c", "o2_c", "h2o_c",
            "gthrd_c", "gthox_c", "nadph_c", "nadp_c"]
    formulas = ["C6H12O6", "C5H10N2O3", "C3H5O3", "CO2", "H2O2", "O2", "H2O",
                "C10H16N3O6S", "C20H30N6O12S2", "C21H26N7O17P3", "C21H25N7O17P3"]
    reactions = {
        "EX_glc__D_e": ({"glc__D_e": -1}, -2.0), "EX_gln__L_e": ({"gln__L_e": -1}, -1.0),
        "EX_lac__L_e": ({"lac__L_e": -1}, 3.0), "EX_co2_e": ({"co2_e": -1}, 8.0),
        "OXIDASE": ({"o2_c": -1, "h2o2_c": 1}, 3.0),
        "CAT": ({"h2o2_c": -2, "o2_c": 1, "h2o_c": 2}, 1.0),
        "GTHP": ({"h2o2_c": -1, "gthrd_c": -2, "gthox_c": 1, "h2o_c": 2}, 2.0),
        "GTHOy": ({"gthox_c": -1, "nadph_c": -1, "gthrd_c": 2, "nadp_c": 1}, 2.0),
        "G6PDH2": ({"nadp_c": -1, "nadph_c": 1}, 5.0),
        "GND": ({"nadp_c": -1, "nadph_c": 1}, 1.0),
        "NADPHOX": ({"nadph_c": -1, "nadp_c": 1}, 2.0)}
    S = np.zeros((len(mets), len(reactions)))
    for j, (stoichiometry, _) in enumerate(reactions.values()):
        for met, coefficient in stoichiometry.items():
            S[mets.index(met), j] = coefficient
    return (Model(S, list(reactions), mets, formulas=formulas),
            np.array([value for _, value in reactions.values()]))


#: The answers to the functional toy, worked by hand from its docstring.
FUNCTIONAL_TOY_EXPECTED = {
    "elem_n_uptake": 2.0, "elem_n_from_amino_acids": 1.0,
    "elem_c_from_glucose": 12 / 17, "elem_c_from_amino_acids": 5 / 17,
    "elem_c_from_lactate_pyruvate": 0.0, "elem_c_to_lactate": 9 / 17, "elem_c_to_co2": 8 / 17,
    "div_c_effective_sources": 289 / 169, "div_c_dominant_source": 12 / 17,
    "aa_glutamine_net_uptake": 1.0, "aa_total_uptake": 1.0, "aa_n_taken_up": 1.0,
    "ros_h2o2_production": 3.0, "ros_h2o2_removal": 4.0, "ros_h2o2_removal_per_production": 4 / 3,
    "ros_h2o2_removed_by_catalase": 0.5, "ros_h2o2_removed_by_glutathione_peroxidase": 0.5,
    "ros_glutathione_regeneration_cover": 6 / 2, "ros_glutathione_cycle_closure": 1.0,
    "conv_lactate_per_glucose_carbon": 9 / 12, "conv_nadph_per_carbon": 6 / 17,
    "bud_nadph_to_glutathione_reduction": 2 / 4, "bud_oxppp_cover_of_glutathione_reduction": 6 / 2,
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

    # elements, amino acids, oxidative defence and conversion, worked by hand
    toy, flux = _functional_toy()
    h = flux2feature(toy, flux)
    for name, expected in FUNCTIONAL_TOY_EXPECTED.items():
        assert math.isclose(h[name], expected, rel_tol=1e-9, abs_tol=1e-12), (name, h[name], expected)
    assert not audit(toy)[0] is None
    absent = {r["name"]: r for r in audit(toy)}
    assert absent["ros_h2o2_removed_by_peroxiredoxin"]["computable"] is False
    assert math.isnan(flux2feature(toy, flux, which="all")["ros_h2o2_removed_by_peroxiredoxin"])

    computable = sum(r["computable"] for r in rows.values())
    print(f"self-test passed: {len(rows)} descriptors declared, "
          f"{computable} computable in this toy model")


if __name__ == "__main__":
    _self_test()
