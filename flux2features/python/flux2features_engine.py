"""flux2features engine: every declared descriptor of Phi_M(v).

GENERATED from notebook/metabolic_flux_phenotyping.ipynb by
python/make_engine_module.py. Do not edit: change the notebook and regenerate.

This is the engine the study ran: it computes all every declared
descriptors, each NaN when the reconstruction lacks what it requires, with
NumPy, SciPy and pandas. The standalone tool (flux2feature.py) is a smaller,
dependency-light subset of it with a MATLAB twin.

    import flux2features_engine as fe
    gem = fe.load_gem("RBC-GEM.mat")          # a COBRA model saved by MATLAB
    table = fe.phi(gem, V)                    # V: states x reactions
    fe.declared(gem)                          # every descriptor, defined or not, and why

The declarations of the study (PATHWAY_GROUPS, TURNOVER_METABOLITES,
KEY_REACTIONS, EXCHANGE_ADMISSIBILITY, ...) describe RBC-GEM. For another
reconstruction, replace them before calling phi: a descriptor whose
requirement the model does not meet is reported as NaN, never as zero.
"""

import pandas as pd

"""What the study fixes before any state is read: the design, the named
reactions and metabolites, the pathway groups, the tolerances and the seed."""
from dataclasses import dataclass

# --------------------------------------------------------------------------
# Experimental design
# --------------------------------------------------------------------------

#: Fractional upper bound applied to G6PDH2 (glucose-6-phosphate dehydrogenase).
G6PD_LEVELS = (0.10, 0.28, 0.50)

#: Fractional upper bound applied to EX_h2o2_e (hydrogen peroxide exchange).
H2O2_LEVELS = (0.50, 0.75, 1.00)

#: WHO-inspired ordinal severity labels for the G6PD axis.
G6PD_SEVERITY = {0.10: "severe", 0.28: "moderate", 0.50: "mild"}

#: Ordinal oxidative-load labels for the H2O2 axis. A *larger* permitted
#: exchange bound means the cell can move more peroxide across its boundary.
H2O2_LOAD = {0.50: "low", 0.75: "intermediate", 1.00: "high"}

N_CHAINS = 3

#: States kept per chain, by reconstruction. Every chain keeps the same
#: number, so no condition outweighs another; the stride is derived from each
#: chain's own length and written to its JSON, so a chain sampled with a
#: different number of points is thinned to the same count rather than
#: silently weighted differently.
KEEP_PER_CHAIN = {"original": 1000, "thermo": 1000}

#: Reactions per reconstruction; used to orient the sampling matrices.
N_REACTIONS = {"original": 3275, "thermo": 2580}


# --------------------------------------------------------------------------
# Key reactions and metabolites of the RBC-GEM
# --------------------------------------------------------------------------

G6PD_REACTION = "G6PDH2"
H2O2_EXCHANGE = "EX_h2o2_e"

#: Metabolites whose turnover defines the biochemical phenotype.
TURNOVER_METABOLITES = (
    "atp_c", "adp_c", "amp_c", "pi_c",
    "nadh_c", "nad_c", "nadph_c", "nadp_c",
    "gthrd_c", "gthox_c",
    "glc__D_c", "g6p_c", "f6p_c", "pyr_c", "lac__L_c",
    "h2o2_c", "o2_c", "o2s_c",
    "13dpg_c", "23dpg_c", "r5p_c", "ru5p__D_c", "6pgc_c",
    "gtp_c", "utp_c", "ctp_c", "fad_c", "fadh2_c",
    "accoa_c", "oaa_c", "akg_c", "glu__L_c", "gln__L_c",
)

#: Currency metabolites excluded from flux-weighted network statistics, where
#: their promiscuous participation would otherwise dominate every centrality.
CURRENCY_METABOLITES = (
    "h_c", "h_e", "h2o_c", "h2o_e", "pi_c", "pi_e", "ppi_c",
    "co2_c", "co2_e", "na1_c", "na1_e", "k_c", "k_e", "cl_c", "cl_e",
)

#: Curated pathway groupings. Keys are the labels used throughout the paper;
#: values are the RBC-GEM ``subSystems`` strings that map onto them. The
#: strings are verbatim from the model, so an unmatched entry is a bug rather
#: than a silent omission (:func:`mfp.annotation.reactions.check_pathway_map`).
PATHWAY_GROUPS: dict[str, tuple[str, ...]] = {
    "glycolysis": ("Glycolysis / Gluconeogenesis",),
    "pentose_phosphate": ("Pentose phosphate pathway",),
    "glutathione": ("Glutathione metabolism",),
    "redox_detox": ("Reactive species formation and detoxification",),
    "purine": ("Purine metabolism",),
    "pyrimidine": ("Pyrimidine metabolism",),
    "nucleotide": (
        "Nucleotide metabolism",
        "Amino sugar and nucleotide sugar metabolism",
    ),
    "pyruvate_tca": (
        "Pyruvate metabolism",
        "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism",
    ),
    "amino_acid": (
        "Alanine, aspartate and glutamate metabolism",
        "Glycine, serine and threonine metabolism",
        "Cysteine and methionine metabolism",
        "Arginine and proline metabolism",
        "Phenylalanine, tyrosine and tryptophan metabolism",
        "Valine, leucine, and isoleucine metabolism",
        "Histidine metabolism",
        "Lysine metabolism",
        "Selenocompound metabolism",
        "Sulfur metabolism",
        "Aminoacyl-tRNA biosynthesis",
        "Peptide metabolism",
    ),
    "lipid": (
        "Fatty acid oxidation",
        "Fatty acid activation (cytosolic)",
        "Fatty acid biosynthesis",
        "Fatty acid biosynthesis (even-chain)",
        "Fatty acid biosynthesis (unsaturated)",
        "Fatty acid elongation (even-chain)",
        "Fatty acid elongation (odd-chain)",
        "Fatty acid desaturation (even-chain)",
        "Acyl-CoA hydrolysis",
        "Carnitine shuttle (cytosolic)",
        "Glycerophospholipid metabolism",
        "Glycerolipid metabolism",
        "Sphingolipid metabolism",
        "Ether lipid metabolism",
        "Cholesterol metabolism",
        "Inositol phosphate metabolism",
        "Docosahexaenoic acid metabolism",
        "Eicosapentaenoic acid metabolism",
        "Arachidonic acid metabolism",
        "Linoleate metabolism",
        "Omega-3 fatty acid metabolism",
        "Omega-6 fatty acid metabolism",
        "Prostaglandin biosynthesis",
        "Leukotriene metabolism",
    ),
    "glycan": (
        "Glycosphingolipid biosynthesis - ganglio series",
        "Glycosphingolipid biosynthesis - lacto and neolacto series",
        "Glycosphingolipid biosynthesis - globo and isoglobo series",
        "N-glycan metabolism",
        "Starch and sucrose metabolism",
        "Galactose metabolism",
        "Fructose and mannose metabolism",
        "Pentose and glucuronate interconversions",
    ),
    "cofactor_vitamin": (
        "Nicotinate and nicotinamide metabolism",
        "Folate metabolism",
        "Riboflavin metabolism",
        "Thiamine metabolism",
        "Vitamin B6 metabolism",
        "Vitamin E metabolism",
        "Retinol metabolism",
        "Biopterin metabolism",
        "Pantothenate and CoA metabolism",
        "Ascorbate and aldarate metabolism",
    ),
    "porphyrin_hemoglobin": (
        "Porphyrin metabolism",
        "Hemoglobin binding and degradation",
    ),
    "xenobiotic_hormone": (
        "Estrogen metabolism",
        "Androgen metabolism",
        "Steroid metabolism",
        "Azathioprine and 6-mercaptopurine metabolism",
        "Ribavirin metabolism",
        "5-fluorouracil metabolism",
        "Renin-angiotensin system",
    ),
    "protein_turnover": (
        "Protein modification",
        "Ubiquitin mediated proteolysis",
    ),
    "transport": ("Transport, extracellular",),
    "pseudoreactions": ("Pseudoreactions", "Miscellaneous"),
}

#: The Rapoport-Luebering shunt is not annotated as a subsystem in RBC-GEM but
#: is central to erythrocyte physiology, so it is defined by reaction identity.
RAPOPORT_LUEBERING_RXNS = ("DPGM", "DPGase")

#: Reactions that define the oxidative branch of the pentose phosphate pathway.
OXIDATIVE_PPP_RXNS = ("G6PDH2", "PGL", "GND")


# --------------------------------------------------------------------------
# Numerical tolerances
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Tolerances:
    #: below this magnitude a flux is treated as numerically zero
    flux_zero: float = 1e-9
    #: below this range a reaction is treated as invariant across samples
    invariant_range: float = 1e-8
    #: relative tolerance used when checking mass balance of a sample
    mass_balance: float = 1e-6
    #: features with |Spearman rho| above this are considered redundant
    redundancy_rho: float = 0.95
    #: features with coefficient of variation below this are considered constant
    min_cv: float = 1e-6


TOL = Tolerances()


# --------------------------------------------------------------------------
# Analysis choices
# --------------------------------------------------------------------------

#: The one seed every stochastic step of the analysis derives its generator
#: from (subsampling, forests, permutations, k-means). The sampling seeds are
#: set separately, one per chain, by the MATLAB script (baseSeed + 100 k + c).
RANDOM_SEED = 20260812

#: A descriptor responds to the gradient if its chain medians correlate with
#: the imposed capacity at least this strongly (and survive FDR at 5 %).
RESPONSIVE_RHO = 0.7
#: Two responses of the same sign are "the same size" within this fold change.
MAGNITUDE_FOLD = 2.0


# The stages read these settings as ``config.X``, so the code reads the same as
# it would in a script.
import types

config = types.SimpleNamespace(
    **{k: v for k, v in globals().items() if k.isupper() or k == "ensure_dirs"})


from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.io as sio
import scipy.sparse as sp


def _as_str_array(value) -> np.ndarray:
    """Normalise the many shapes MATLAB cell arrays take after ``loadmat``."""
    out = []
    for item in np.atleast_1d(value):
        if isinstance(item, np.ndarray):
            flat = item.ravel()
            out.append(str(flat[0]) if flat.size else "")
        else:
            out.append(str(item))
    return np.array(out, dtype=object)


@dataclass
class GEM:
    """A minimal, analysis-oriented view of a COBRA model."""

    name: str
    S: sp.csr_matrix
    lb: np.ndarray
    ub: np.ndarray
    b: np.ndarray
    c: np.ndarray
    rxns: np.ndarray
    rxn_names: np.ndarray
    subsystems: np.ndarray
    mets: np.ndarray
    met_names: np.ndarray
    met_formulas: np.ndarray
    met_charges: np.ndarray
    genes: np.ndarray
    gene_names: np.ndarray
    rxn_gene_mat: sp.csr_matrix | None
    gr_rules: np.ndarray

    # -- basic properties -------------------------------------------------
    @property
    def n_reactions(self) -> int:
        return self.S.shape[1]

    @property
    def n_metabolites(self) -> int:
        return self.S.shape[0]

    def rxn_index(self, rxn_id: str) -> int:
        hits = np.flatnonzero(self.rxns == rxn_id)
        if hits.size == 0:
            raise KeyError(f"reaction {rxn_id!r} is not in the model")
        return int(hits[0])

    def met_index(self, met_id: str) -> int:
        hits = np.flatnonzero(self.mets == met_id)
        if hits.size == 0:
            raise KeyError(f"metabolite {met_id!r} is not in the model")
        return int(hits[0])

    def has_met(self, met_id: str) -> bool:
        return bool(np.any(self.mets == met_id))

    def has_rxn(self, rxn_id: str) -> bool:
        return bool(np.any(self.rxns == rxn_id))

    def with_bounds(self, changes: dict[str, tuple[float | None, float | None]]) -> "GEM":
        """Return a copy of the model with selected bounds replaced."""
        lb, ub = self.lb.copy(), self.ub.copy()
        for rxn_id, (low, high) in changes.items():
            idx = self.rxn_index(rxn_id)
            if low is not None:
                lb[idx] = low
            if high is not None:
                ub[idx] = high
        clone = GEM(**{**self.__dict__, "lb": lb, "ub": ub})
        return clone

    def carbon_count(self) -> np.ndarray:
        """Number of carbon atoms per metabolite, parsed from ``metFormulas``."""
        import re

        counts = np.zeros(self.n_metabolites)
        pattern = re.compile(r"C(\d*)(?![a-z])")
        for i, formula in enumerate(self.met_formulas):
            if not isinstance(formula, str) or not formula:
                continue
            match = pattern.search(formula)
            if match:
                counts[i] = int(match.group(1)) if match.group(1) else 1
        return counts


def load_gem(path: str | Path, name: str | None = None) -> GEM:
    """Read a COBRA model stored as a MATLAB v7 ``.mat`` file."""
    path = Path(path)
    raw = sio.loadmat(str(path), struct_as_record=False, squeeze_me=True)
    keys = [k for k in raw if not k.startswith("__")]
    if len(keys) != 1:
        raise ValueError(f"expected exactly one variable in {path.name}, found {keys}")
    model = raw[keys[0]]

    S = sp.csr_matrix(model.S)
    n_rxn = S.shape[1]

    def field(attr, default=None):
        return getattr(model, attr, default)

    rgm = field("rxnGeneMat")
    return GEM(
        name=name or keys[0],
        S=S,
        lb=np.asarray(model.lb, dtype=float).ravel(),
        ub=np.asarray(model.ub, dtype=float).ravel(),
        b=np.asarray(field("b", np.zeros(S.shape[0])), dtype=float).ravel(),
        c=np.asarray(field("c", np.zeros(n_rxn)), dtype=float).ravel(),
        rxns=_as_str_array(model.rxns),
        rxn_names=_as_str_array(field("rxnNames", model.rxns)),
        subsystems=_as_str_array(field("subSystems", np.array([""] * n_rxn))),
        mets=_as_str_array(model.mets),
        met_names=_as_str_array(field("metNames", model.mets)),
        met_formulas=_as_str_array(field("metFormulas", np.array([""] * S.shape[0]))),
        met_charges=np.asarray(field("metCharges", np.zeros(S.shape[0])), dtype=float).ravel(),
        genes=_as_str_array(field("genes", np.array([]))),
        gene_names=_as_str_array(field("geneNames", np.array([]))),
        rxn_gene_mat=sp.csr_matrix(rgm) if rgm is not None else None,
        gr_rules=_as_str_array(field("grRules", np.array([""] * n_rxn))),
    )


def load_flux_vector(path: str | Path) -> np.ndarray:
    """Read a plain numeric vector stored in a ``.mat`` file (e.g. ``v_or.mat``)."""
    raw = sio.loadmat(str(path), struct_as_record=False, squeeze_me=True)
    keys = [k for k in raw if not k.startswith("__")]
    if len(keys) != 1:
        raise ValueError(f"expected exactly one variable in {path}, found {keys}")
    return np.asarray(raw[keys[0]], dtype=float).ravel()


def load_fba_solution(path: str | Path) -> dict:
    """Read a COBRA ``optimizeCbModel`` solution structure."""
    raw = sio.loadmat(str(path), struct_as_record=False, squeeze_me=True)
    keys = [k for k in raw if not k.startswith("__")]
    sol = raw[keys[0]]
    out = {"objective": float(np.atleast_1d(sol.f)[0])}
    for attr in ("v", "x"):
        if hasattr(sol, attr):
            out["v"] = np.asarray(getattr(sol, attr), dtype=float).ravel()
            break
    out["status"] = int(np.atleast_1d(getattr(sol, "stat", 1))[0])
    return out


import re
from dataclasses import dataclass

import numpy as np
import pandas as pd


class ReactionClass:
    """Mutually exclusive structural roles a reaction can take."""

    EXCHANGE = "exchange"          # EX_*: mass crossing the system boundary
    SINK = "sink"                  # SK_*: intracellular pseudo-pools
    TRANSPORT = "transport"        # moves a metabolite between compartments
    METABOLIC = "metabolic"        # genuine intracellular conversion
    PSEUDO = "pseudoreaction"      # lumped/pseudo processes that are not sinks


# --------------------------------------------------------------------------
# Biological admissibility of boundary reactions
# --------------------------------------------------------------------------
#
# ``role`` answers "what is this compound to an erythrocyte in plasma?" and
# ``modifiable`` answers "could a clinician plausibly change its availability?".
# ``toxicity`` is an ordinal penalty weight in [0, 1] used by the intervention
# utility function; it encodes the harm of *increasing* exposure to the
# metabolite, and is what prevents mathematically attractive but clinically
# absurd interventions from being proposed.

@dataclass(frozen=True)
class Admissibility:
    role: str
    modifiable: bool
    toxicity: float
    note: str = ""


#: Curated admissibility for the boundary metabolites that carry the biology of
#: G6PD deficiency. Everything not listed falls back to :func:`_default_admissibility`.
EXCHANGE_ADMISSIBILITY: dict[str, Admissibility] = {
    # --- primary carbon and energy ------------------------------------
    "glc__D_e": Admissibility("nutrient", True, 0.05, "principal RBC carbon source"),
    "lac__L_e": Admissibility("waste", True, 0.10, "glycolytic end product"),
    "pyr_e": Admissibility("nutrient", True, 0.10, "alternative carbon source"),
    "fru_e": Admissibility("nutrient", True, 0.10, ""),
    "gal_e": Admissibility("nutrient", True, 0.10, ""),
    "man_e": Admissibility("nutrient", True, 0.10, ""),
    "rib__D_e": Admissibility("nutrient", True, 0.10, "feeds non-oxidative PPP"),
    "adn_e": Admissibility("nutrient", True, 0.15, "adenosine salvage"),
    "ade_e": Admissibility("nutrient", True, 0.15, "adenine salvage"),
    "ins_e": Admissibility("nutrient", True, 0.15, "inosine salvage"),
    "pi_e": Admissibility("nutrient", True, 0.05, ""),
    # --- redox and oxidative stress -----------------------------------
    "h2o2_e": Admissibility("toxicant", True, 0.95, "oxidative insult"),
    "o2_e": Admissibility("nutrient", True, 0.20, "substrate for ROS generation"),
    "o2s_e": Admissibility("toxicant", False, 1.00, "superoxide"),
    "no_e": Admissibility("signalling", False, 0.60, ""),
    "gthrd_e": Admissibility("therapeutic", True, 0.05, "reduced glutathione"),
    "gthox_e": Admissibility("waste", True, 0.30, "oxidised glutathione"),
    "cys__L_e": Admissibility("therapeutic", True, 0.10, "glutathione precursor"),
    "gly_e": Admissibility("therapeutic", True, 0.05, "glutathione precursor"),
    "glu__L_e": Admissibility("therapeutic", True, 0.10, "glutathione precursor"),
    "glycys_e": Admissibility("therapeutic", True, 0.10, "gamma-glutamylcysteine"),
    "ascb__L_e": Admissibility("therapeutic", True, 0.10, "vitamin C"),
    "dhdascb_e": Admissibility("waste", True, 0.25, "dehydroascorbate"),
    "nac_e": Admissibility("therapeutic", True, 0.10, "nicotinate, NAD precursor"),
    "ncam_e": Admissibility("therapeutic", True, 0.10, "nicotinamide, NAD precursor"),
    "ribflv_e": Admissibility("therapeutic", True, 0.10, "riboflavin, FAD precursor"),
    "pydxn_e": Admissibility("therapeutic", True, 0.10, "vitamin B6"),
    "thm_e": Admissibility("therapeutic", True, 0.10, "thiamine, transketolase cofactor"),
    "avite1_e": Admissibility("therapeutic", True, 0.10, "alpha-tocopherol"),
    "selni_e": Admissibility("therapeutic", True, 0.55, "selenium, narrow therapeutic index"),
    # --- known haemolytic triggers in G6PD deficiency -------------------
    "pchlor_e": Admissibility("toxicant", False, 1.00, "primaquine-type oxidant"),
    "meoh_e": Admissibility("toxicant", False, 0.90, ""),
    "cyan_e": Admissibility("toxicant", False, 1.00, "cyanide"),
    "co_e": Admissibility("toxicant", False, 0.95, "carbon monoxide"),
    "nh4_e": Admissibility("waste", True, 0.45, ""),
    "urea_e": Admissibility("waste", True, 0.10, ""),
    "co2_e": Admissibility("waste", True, 0.05, ""),
    "h_e": Admissibility("ion", True, 0.20, "pH"),
    "h2o_e": Admissibility("solvent", False, 0.00, ""),
}

#: Metabolite name fragments that mark a boundary compound as pharmacological
#: or xenobiotic rather than nutritional. Increasing exposure to these is
#: penalised because they are drugs or environmental agents, not nutrients.
_XENOBIOTIC_PATTERNS = (
    "azathioprine", "mercaptopurine", "ribavirin", "fluorouracil",
    "estr", "andro", "steroid", "prostaglandin", "leukotriene",
    "primaquine", "dapsone", "sulfa", "nitro",
)


def _default_admissibility(met_id: str, met_name: str, formula: str) -> Admissibility:
    lowered = f"{met_id} {met_name}".lower()
    if any(pat in lowered for pat in _XENOBIOTIC_PATTERNS):
        return Admissibility("xenobiotic", False, 0.75, "pharmacological agent")
    if re.search(r"\bion\b|^na1_|^k_|^cl_|^ca2_|^mg2_|^fe2_|^fe3_|^zn2_|^cu2_", met_id):
        return Admissibility("ion", True, 0.30, "")
    # An uncharacterised organic boundary metabolite is treated as a possible
    # nutrient but is not considered clinically modifiable, so it can never be
    # promoted into an intervention without explicit curation.
    if formula and "C" in formula:
        return Admissibility("metabolite", False, 0.40, "uncurated boundary metabolite")
    return Admissibility("metabolite", False, 0.40, "uncurated boundary metabolite")


# --------------------------------------------------------------------------
# Annotation
# --------------------------------------------------------------------------

def _pathway_lookup() -> dict[str, str]:
    out: dict[str, str] = {}
    for pathway, subsystems in PATHWAY_GROUPS.items():
        for sub in subsystems:
            out[sub] = pathway
    return out


def check_pathway_map(gem: GEM) -> pd.DataFrame:
    """Report model subsystems that no curated pathway group claims.

    Unmapped subsystems are not an error, but they must be visible: they are
    reported in the supplementary material so that the pathway-level feature
    family can be audited for coverage.
    """
    lookup = _pathway_lookup()
    unmapped = sorted({s for s in gem.subsystems if s and s not in lookup})
    counts = pd.Series(gem.subsystems).value_counts()
    return pd.DataFrame(
        {"subsystem": unmapped, "n_reactions": [int(counts.get(s, 0)) for s in unmapped]}
    ).sort_values("n_reactions", ascending=False, ignore_index=True)


def annotate_reactions(gem: GEM) -> pd.DataFrame:
    """Build the reaction metadata table.

    Returns one row per reaction with structural role, compartment span,
    pathway assignment, bound information and — for boundary reactions — the
    biological admissibility record used by the intervention framework.
    """
    S = gem.S.tocsc()
    n_rxn = gem.n_reactions
    lookup = _pathway_lookup()

    met_comp = np.array([m.rsplit("_", 1)[-1] for m in gem.mets], dtype=object)

    rows = []
    for j in range(n_rxn):
        rxn = gem.rxns[j]
        col = S[:, [j]].tocoo()
        met_idx = col.row
        coeffs = col.data
        comps = sorted(set(met_comp[met_idx].tolist()))
        subsystem = gem.subsystems[j]

        if rxn.startswith("EX_"):
            cls = ReactionClass.EXCHANGE
        elif rxn.startswith(("SK_", "sink_", "DM_")):
            cls = ReactionClass.SINK
        elif len(met_idx) == 1:
            cls = ReactionClass.SINK
        elif len(comps) > 1:
            cls = ReactionClass.TRANSPORT
        elif subsystem in ("Pseudoreactions", "Miscellaneous"):
            cls = ReactionClass.PSEUDO
        else:
            cls = ReactionClass.METABOLIC

        boundary_met = ""
        boundary_name = ""
        boundary_formula = ""
        if cls == ReactionClass.EXCHANGE and len(met_idx) == 1:
            boundary_met = gem.mets[met_idx[0]]
            boundary_name = gem.met_names[met_idx[0]]
            boundary_formula = gem.met_formulas[met_idx[0]]

        pathway = lookup.get(subsystem, "unassigned")
        if rxn in RAPOPORT_LUEBERING_RXNS:
            pathway = "rapoport_luebering"

        rows.append(
            {
                "rxn_id": rxn,
                "rxn_name": gem.rxn_names[j],
                "subsystem": subsystem,
                "pathway": pathway,
                "reaction_class": cls,
                "n_metabolites": int(len(met_idx)),
                "n_compartments": len(comps),
                "compartments": "|".join(comps),
                "lb": float(gem.lb[j]),
                "ub": float(gem.ub[j]),
                "reversible": bool(gem.lb[j] < 0 < gem.ub[j]),
                "objective_coefficient": float(gem.c[j]),
                "boundary_metabolite": boundary_met,
                "boundary_metabolite_name": boundary_name,
                "boundary_formula": boundary_formula,
                "stoich_sum_abs": float(np.abs(coeffs).sum()),
                "gene_rule": gem.gr_rules[j] if j < len(gem.gr_rules) else "",
                "is_oxidative_ppp": rxn in OXIDATIVE_PPP_RXNS,
            }
        )

    table = pd.DataFrame(rows)

    # -- admissibility for boundary reactions ---------------------------
    roles, modifiable, toxicity, notes = [], [], [], []
    for _, row in table.iterrows():
        if row["reaction_class"] != ReactionClass.EXCHANGE:
            roles.append("")
            modifiable.append(False)
            toxicity.append(np.nan)
            notes.append("")
            continue
        met = row["boundary_metabolite"]
        adm = EXCHANGE_ADMISSIBILITY.get(met)
        if adm is None:
            adm = _default_admissibility(
                met, row["boundary_metabolite_name"], row["boundary_formula"]
            )
        roles.append(adm.role)
        modifiable.append(adm.modifiable)
        toxicity.append(adm.toxicity)
        notes.append(adm.note)

    table["boundary_role"] = roles
    table["clinically_modifiable"] = modifiable
    table["toxicity_weight"] = toxicity
    table["admissibility_note"] = notes

    # -- the internal/external decomposition ----------------------------
    #
    # "External" is deliberately narrow: only reactions that exchange mass with
    # the extracellular environment qualify, because only those correspond to
    # something a clinician could in principle manipulate. Intracellular sink
    # pseudo-reactions are boundary reactions in the linear-algebraic sense but
    # are not environmental, so they are kept on the internal side and flagged.
    table["is_external"] = table["reaction_class"] == ReactionClass.EXCHANGE
    table["is_internal"] = ~table["is_external"]
    table["carbon_count"] = [
        _carbon_count(f) for f in table["boundary_formula"].fillna("")
    ]
    return table


def _carbon_count(formula: str) -> float:
    if not isinstance(formula, str) or not formula:
        return np.nan
    match = re.search(r"C(\d*)(?![a-z])", formula)
    if not match:
        return 0.0
    return float(match.group(1)) if match.group(1) else 1.0


def internal_mask(reactions: pd.DataFrame) -> np.ndarray:
    return reactions["is_internal"].to_numpy(dtype=bool)


def external_mask(reactions: pd.DataFrame) -> np.ndarray:
    return reactions["is_external"].to_numpy(dtype=bool)


from dataclasses import dataclass, field
from typing import Callable, Iterable

import numpy as np
import pandas as pd
import scipy.sparse as sp


@dataclass(frozen=True)
class FeatureSpec:
    """Declarative description of one scalar feature of a flux vector."""

    name: str
    family: str
    block: str          # "internal" | "external" | "combined"
    description: str
    formula: str        # LaTeX, without surrounding delimiters
    units: str          # "flux", "dimensionless", "nats", "count", "ratio"
    signed: bool        # computed on v (True) or on |v| (False)
    model_dependent: bool = False

    def as_record(self) -> dict:
        return {
            "feature": self.name,
            "family": self.family,
            "block": self.block,
            "description": self.description,
            "formula": self.formula,
            "units": self.units,
            "signed_input": self.signed,
            "model_dependent": self.model_dependent,
        }


@dataclass
class FeatureContext:
    """Everything a feature family needs to know about the model."""

    gem: GEM
    reactions: pd.DataFrame

    # derived, filled in __post_init__
    internal_idx: np.ndarray = field(init=False)
    external_idx: np.ndarray = field(init=False)
    S_int: sp.csr_matrix = field(init=False)
    S_ext: sp.csr_matrix = field(init=False)

    def __post_init__(self) -> None:
        self.internal_idx = np.flatnonzero(internal_mask(self.reactions))
        self.external_idx = np.flatnonzero(external_mask(self.reactions))
        S = self.gem.S.tocsc()
        self.S_int = sp.csr_matrix(S[:, self.internal_idx])
        self.S_ext = sp.csr_matrix(S[:, self.external_idx])

    # -- convenience ----------------------------------------------------
    def split(self, V: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Split a ``(n_samples, n_reactions)`` matrix into internal/external."""
        return V[:, self.internal_idx], V[:, self.external_idx]

    def internal_reactions(self) -> pd.DataFrame:
        return self.reactions.iloc[self.internal_idx].reset_index(drop=True)

    def external_reactions(self) -> pd.DataFrame:
        return self.reactions.iloc[self.external_idx].reset_index(drop=True)

    def pathway_groups(self, block: str = "internal") -> dict[str, np.ndarray]:
        """Map pathway label -> positional indices within the given block."""
        table = self.internal_reactions() if block == "internal" else self.external_reactions()
        return {
            pathway: np.asarray(idx, dtype=int)
            for pathway, idx in table.groupby("pathway").groups.items()
        }


class FeatureFamily:
    """Base class for a group of features computed together."""

    name: str = "unnamed"

    def __init__(self, ctx: FeatureContext):
        self.ctx = ctx

    def specs(self) -> list[FeatureSpec]:  # pragma: no cover - interface
        raise NotImplementedError

    def compute(self, V_int: np.ndarray, V_ext: np.ndarray) -> pd.DataFrame:  # pragma: no cover
        raise NotImplementedError


class FeatureRegistry:
    """Ordered collection of feature specifications."""

    def __init__(self) -> None:
        self._specs: dict[str, FeatureSpec] = {}

    def add(self, spec: FeatureSpec) -> None:
        if spec.name in self._specs:
            raise ValueError(f"duplicate feature name {spec.name!r}")
        self._specs[spec.name] = spec

    def extend(self, specs: Iterable[FeatureSpec]) -> None:
        for spec in specs:
            self.add(spec)

    def __contains__(self, name: str) -> bool:
        return name in self._specs

    def __len__(self) -> int:
        return len(self._specs)

    @property
    def names(self) -> list[str]:
        return list(self._specs)

    def get(self, name: str) -> FeatureSpec:
        return self._specs[name]

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame([s.as_record() for s in self._specs.values()])

    def subset(self, *, block: str | None = None, family: str | None = None) -> list[str]:
        out = []
        for name, spec in self._specs.items():
            if block is not None and spec.block != block:
                continue
            if family is not None and spec.family != family:
                continue
            out.append(name)
        return out


# --------------------------------------------------------------------------
# Small numerical helpers shared by several families
# --------------------------------------------------------------------------

def safe_divide(num: np.ndarray, den: np.ndarray, fill: float = np.nan) -> np.ndarray:
    """Element-wise division that returns ``fill`` where the denominator vanishes."""
    num = np.asarray(num, dtype=float)
    den = np.asarray(den, dtype=float)
    out = np.full(np.broadcast(num, den).shape, fill, dtype=float)
    ok = np.abs(den) > 0
    np.divide(num, den, out=out, where=ok)
    return out


def positive_part(V: np.ndarray) -> np.ndarray:
    return np.where(V > 0, V, 0.0)


def negative_part(V: np.ndarray) -> np.ndarray:
    """Magnitude of the negative entries (returned as a non-negative array)."""
    return np.where(V < 0, -V, 0.0)


def gini(A: np.ndarray) -> np.ndarray:
    """Gini coefficient of each row of a non-negative matrix."""
    A = np.sort(np.abs(A), axis=1)
    n = A.shape[1]
    index = np.arange(1, n + 1)
    total = A.sum(axis=1)
    numerator = ((2 * index - n - 1) * A).sum(axis=1)
    return safe_divide(numerator, n * total, fill=0.0)


def hoyer_sparsity(A: np.ndarray) -> np.ndarray:
    """Hoyer sparsity in [0, 1]; 0 = fully dense, 1 = a single non-zero."""
    n = A.shape[1]
    l1 = np.abs(A).sum(axis=1)
    l2 = np.sqrt((A ** 2).sum(axis=1))
    root_n = np.sqrt(n)
    return safe_divide(root_n - safe_divide(l1, l2, fill=0.0), root_n - 1.0, fill=0.0)


def shannon_entropy(A: np.ndarray) -> np.ndarray:
    """Shannon entropy (nats) of the flux-magnitude distribution of each row."""
    P = np.abs(A)
    total = P.sum(axis=1, keepdims=True)
    P = safe_divide(P, total, fill=0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = np.where(P > 0, P * np.log(P), 0.0)
    return -terms.sum(axis=1)


def renyi_entropy(A: np.ndarray, alpha: float = 2.0) -> np.ndarray:
    """Renyi entropy of order ``alpha`` (nats)."""
    P = np.abs(A)
    total = P.sum(axis=1, keepdims=True)
    P = safe_divide(P, total, fill=0.0)
    powered = np.power(P, alpha, where=P > 0, out=np.zeros_like(P)).sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.log(powered) / (1.0 - alpha)
    return np.where(powered > 0, out, 0.0)


def top_fraction(A: np.ndarray, fraction: float) -> np.ndarray:
    """Share of total flux magnitude carried by the top ``fraction`` reactions."""
    M = np.abs(A)
    n = M.shape[1]
    k = max(1, int(np.ceil(fraction * n)))
    part = np.partition(M, n - k, axis=1)[:, n - k:]
    return safe_divide(part.sum(axis=1), M.sum(axis=1), fill=0.0)


import numpy as np
import pandas as pd
from scipy import stats


_QUANTILES = (0.05, 0.25, 0.50, 0.75, 0.95)


def structural_block(V: np.ndarray, prefix: str, *, reduced: bool = False) -> pd.DataFrame:
    """Compute the structural feature block for one flux matrix.

    Parameters
    ----------
    V
        ``(n_samples, n_reactions)`` signed flux matrix.
    prefix
        Prepended to every column name, e.g. ``"int"`` or ``"ext"``.
    reduced
        If True, only the subset of features that remain interpretable for a
        boundary flux vector is returned.
    """
    A = np.abs(V)
    n = V.shape[1]
    total = A.sum(axis=1)
    out: dict[str, np.ndarray] = {}

    # -- magnitude / norms ---------------------------------------------
    out["L0"] = (A > TOL.flux_zero).sum(axis=1).astype(float)
    out["L1"] = total
    out["L2"] = np.sqrt((V ** 2).sum(axis=1))
    out["Linf"] = A.max(axis=1)
    out["RMS"] = np.sqrt((V ** 2).mean(axis=1))

    # -- distribution of the magnitude vector ---------------------------
    out["mean_abs"] = A.mean(axis=1)
    out["median_abs"] = np.median(A, axis=1)
    out["sd_abs"] = A.std(axis=1, ddof=1)
    out["max_abs"] = out["Linf"]
    out["iqr_abs"] = np.subtract(*np.percentile(A, [75, 25], axis=1))
    out["mad_abs"] = np.median(np.abs(A - np.median(A, axis=1, keepdims=True)), axis=1)
    out["cv_abs"] = safe_divide(out["sd_abs"], out["mean_abs"])
    for q in _QUANTILES:
        out[f"q{int(q * 100):02d}_abs"] = np.quantile(A, q, axis=1)

    # -- distribution of the signed vector ------------------------------
    out["mean_signed"] = V.mean(axis=1)
    out["median_signed"] = np.median(V, axis=1)
    out["sd_signed"] = V.std(axis=1, ddof=1)
    out["skewness_signed"] = stats.skew(V, axis=1, bias=False)
    out["kurtosis_signed"] = stats.kurtosis(V, axis=1, bias=False)

    # -- sparsity -------------------------------------------------------
    active = A > TOL.flux_zero
    out["fraction_zero"] = 1.0 - active.mean(axis=1)
    out["fraction_active"] = active.mean(axis=1)
    out["hoyer_sparsity"] = hoyer_sparsity(V)
    out["gini"] = gini(V)

    # -- entropy of the flux-magnitude distribution ---------------------
    H = shannon_entropy(V)
    out["shannon_entropy"] = H
    out["normalised_entropy"] = H / np.log(n) if n > 1 else np.zeros_like(H)
    out["renyi_entropy_2"] = renyi_entropy(V, alpha=2.0)

    # -- flux architecture ----------------------------------------------
    pos = positive_part(V)
    neg = negative_part(V)
    out["fraction_positive"] = (V > TOL.flux_zero).mean(axis=1)
    out["fraction_negative"] = (V < -TOL.flux_zero).mean(axis=1)
    out["positive_flux"] = pos.sum(axis=1)
    out["negative_flux"] = neg.sum(axis=1)
    out["directionality"] = safe_divide(pos.sum(axis=1) - neg.sum(axis=1), total, fill=0.0)
    out["dominant_share"] = safe_divide(out["Linf"], total, fill=0.0)
    out["top1pct_share"] = top_fraction(V, 0.01)
    out["top5pct_share"] = top_fraction(V, 0.05)
    out["top10pct_share"] = top_fraction(V, 0.10)
    P = safe_divide(A, total[:, None], fill=0.0)
    out["hhi_concentration"] = (P ** 2).sum(axis=1)
    out["effective_reactions"] = np.exp(H)

    if reduced:
        keep = [
            "L0", "L1", "L2", "Linf", "mean_abs", "median_abs", "sd_abs",
            "cv_abs", "fraction_active", "gini", "shannon_entropy",
            "normalised_entropy", "fraction_positive", "fraction_negative",
            "positive_flux", "negative_flux", "directionality",
            "dominant_share", "hhi_concentration", "effective_reactions",
        ]
        out = {k: out[k] for k in keep}

    return pd.DataFrame({f"{prefix}_{k}": v for k, v in out.items()})


_DEFS: list[tuple[str, str, str, str, bool]] = [
    # (suffix, description, latex, units, signed)
    ("L0", "number of active reactions", r"\lVert v\rVert_0=\#\{i:|v_i|>\tau\}", "count", True),
    ("L1", "total flux magnitude", r"\lVert v\rVert_1=\sum_i|v_i|", "flux", False),
    ("L2", "Euclidean flux norm", r"\lVert v\rVert_2=\sqrt{\sum_i v_i^2}", "flux", True),
    ("Linf", "largest single flux", r"\lVert v\rVert_\infty=\max_i|v_i|", "flux", False),
    ("RMS", "root-mean-square flux", r"\sqrt{n^{-1}\sum_i v_i^2}", "flux", True),
    ("mean_abs", "mean flux magnitude", r"n^{-1}\sum_i|v_i|", "flux", False),
    ("median_abs", "median flux magnitude", r"\mathrm{med}_i|v_i|", "flux", False),
    ("sd_abs", "dispersion of flux magnitudes", r"\mathrm{sd}_i(|v_i|)", "flux", False),
    ("max_abs", "maximum flux magnitude", r"\max_i|v_i|", "flux", False),
    ("iqr_abs", "interquartile range of magnitudes", r"Q_{75}(|v|)-Q_{25}(|v|)", "flux", False),
    ("mad_abs", "median absolute deviation of magnitudes", r"\mathrm{med}\,||v|-\mathrm{med}|v||", "flux", False),
    ("cv_abs", "coefficient of variation of magnitudes", r"\mathrm{sd}(|v|)/\overline{|v|}", "dimensionless", False),
    ("mean_signed", "mean signed flux", r"n^{-1}\sum_i v_i", "flux", True),
    ("median_signed", "median signed flux", r"\mathrm{med}_i v_i", "flux", True),
    ("sd_signed", "dispersion of signed fluxes", r"\mathrm{sd}_i(v_i)", "flux", True),
    ("skewness_signed", "skewness of the signed flux distribution", r"\gamma_1(v)", "dimensionless", True),
    ("kurtosis_signed", "excess kurtosis of the signed flux distribution", r"\gamma_2(v)", "dimensionless", True),
    ("fraction_zero", "fraction of inactive reactions", r"n^{-1}\#\{i:|v_i|\le\tau\}", "ratio", False),
    ("fraction_active", "fraction of active reactions", r"n^{-1}\#\{i:|v_i|>\tau\}", "ratio", False),
    ("hoyer_sparsity", "Hoyer sparsity", r"\frac{\sqrt n-\lVert v\rVert_1/\lVert v\rVert_2}{\sqrt n-1}", "dimensionless", True),
    ("gini", "Gini coefficient of flux magnitudes", r"G(|v|)", "dimensionless", False),
    ("shannon_entropy", "entropy of the flux-magnitude distribution", r"H=-\sum_i p_i\ln p_i,\ p_i=|v_i|/\lVert v\rVert_1", "nats", False),
    ("normalised_entropy", "entropy normalised by its maximum", r"H/\ln n", "dimensionless", False),
    ("renyi_entropy_2", "Renyi entropy of order two", r"-\ln\sum_i p_i^2", "nats", False),
    ("fraction_positive", "fraction of reactions in the forward direction", r"n^{-1}\#\{i:v_i>\tau\}", "ratio", True),
    ("fraction_negative", "fraction of reactions in the reverse direction", r"n^{-1}\#\{i:v_i<-\tau\}", "ratio", True),
    ("positive_flux", "total forward flux", r"\sum_i v_i^+", "flux", True),
    ("negative_flux", "total reverse flux", r"\sum_i v_i^-", "flux", True),
    ("directionality", "net directional bias", r"\frac{\sum_i v_i^+-\sum_i v_i^-}{\sum_i|v_i|}", "dimensionless", True),
    ("dominant_share", "share carried by the single largest reaction", r"\lVert v\rVert_\infty/\lVert v\rVert_1", "ratio", False),
    ("top1pct_share", "share carried by the top 1% of reactions", r"\sum_{i\in T_{1\%}}|v_i|/\lVert v\rVert_1", "ratio", False),
    ("top5pct_share", "share carried by the top 5% of reactions", r"\sum_{i\in T_{5\%}}|v_i|/\lVert v\rVert_1", "ratio", False),
    ("top10pct_share", "share carried by the top 10% of reactions", r"\sum_{i\in T_{10\%}}|v_i|/\lVert v\rVert_1", "ratio", False),
    ("hhi_concentration", "Herfindahl-Hirschman flux concentration", r"\sum_i p_i^2", "dimensionless", False),
    ("effective_reactions", "effective number of flux-carrying reactions", r"\exp(H)", "count", False),
]
_DEFS += [
    (f"q{int(q * 100):02d}_abs", f"{int(q * 100)}th percentile of flux magnitudes",
     rf"Q_{{{int(q * 100)}}}(|v|)", "flux", False)
    for q in _QUANTILES
]

_REDUCED_KEEP = {
    "L0", "L1", "L2", "Linf", "mean_abs", "median_abs", "sd_abs", "cv_abs",
    "fraction_active", "gini", "shannon_entropy", "normalised_entropy",
    "fraction_positive", "fraction_negative", "positive_flux",
    "negative_flux", "directionality", "dominant_share",
    "hhi_concentration", "effective_reactions",
}


class StructuralFamily(FeatureFamily):
    """Mathematical, statistical, sparsity, entropy and architecture features."""

    name = "structural"

    def __init__(self, ctx: FeatureContext, *, external_reduced: bool = True):
        super().__init__(ctx)
        self.external_reduced = external_reduced

    def specs(self) -> list[FeatureSpec]:
        out: list[FeatureSpec] = []
        for suffix, desc, formula, units, signed in _DEFS:
            family = _family_of(suffix)
            out.append(FeatureSpec(f"int_{suffix}", family, "internal",
                                   f"internal fluxes: {desc}", formula, units, signed))
        for suffix, desc, formula, units, signed in _DEFS:
            if self.external_reduced and suffix not in _REDUCED_KEEP:
                continue
            family = _family_of(suffix)
            out.append(FeatureSpec(f"ext_{suffix}", family, "external",
                                   f"exchange fluxes: {desc}", formula, units, signed))
        return out

    def compute(self, V_int: np.ndarray, V_ext: np.ndarray) -> pd.DataFrame:
        left = structural_block(V_int, "int", reduced=False)
        right = structural_block(V_ext, "ext", reduced=self.external_reduced)
        return pd.concat([left, right], axis=1)


def _family_of(suffix: str) -> str:
    if suffix in {"L0", "L1", "L2", "Linf", "RMS"}:
        return "magnitude"
    if suffix.startswith(("q0", "q9")) or suffix in {
        "mean_abs", "median_abs", "sd_abs", "max_abs", "iqr_abs", "mad_abs",
        "cv_abs", "mean_signed", "median_signed", "sd_signed",
        "skewness_signed", "kurtosis_signed",
    }:
        return "distribution"
    if suffix in {"fraction_zero", "fraction_active", "hoyer_sparsity", "gini"}:
        return "sparsity"
    if "entropy" in suffix:
        return "entropy"
    return "architecture"


import numpy as np
import pandas as pd
import scipy.sparse as sp


# --------------------------------------------------------------------------
# F3 - pathway activity
# --------------------------------------------------------------------------

class PathwayFamily(FeatureFamily):
    """Absolute activity, relative share, directionality and engagement per pathway."""

    name = "pathway"

    #: Pathways carrying too few reactions to support a directionality statistic
    #: are still reported for activity, but flagged in the registry.
    MIN_REACTIONS = 2

    def __init__(self, ctx: FeatureContext):
        super().__init__(ctx)
        self.groups = {
            name: idx
            for name, idx in ctx.pathway_groups("internal").items()
            if name not in ("unassigned",)
        }

    def specs(self) -> list[FeatureSpec]:
        out = []
        for pathway, idx in self.groups.items():
            n = len(idx)
            out.append(FeatureSpec(
                f"path_{pathway}_activity", "pathway", "internal",
                f"total flux magnitude carried by {pathway} ({n} reactions)",
                r"F_P=\sum_{i\in P}|v_i|", "flux", False))
            out.append(FeatureSpec(
                f"path_{pathway}_share", "pathway", "internal",
                f"share of internal flux carried by {pathway}",
                r"R_P=F_P/\sum_i|v_i|", "ratio", False))
            out.append(FeatureSpec(
                f"path_{pathway}_engagement", "pathway", "internal",
                f"fraction of {pathway} reactions carrying flux",
                r"A_P=|P|^{-1}\#\{i\in P:|v_i|>\tau\}", "ratio", False))
            if n >= self.MIN_REACTIONS:
                out.append(FeatureSpec(
                    f"path_{pathway}_directionality", "pathway", "internal",
                    f"net directional bias of {pathway}",
                    r"D_P=\frac{\sum_{i\in P}v_i^+-\sum_{i\in P}v_i^-}{\sum_{i\in P}|v_i|}",
                    "dimensionless", True))
        return out

    def compute(self, V_int: np.ndarray, V_ext: np.ndarray) -> pd.DataFrame:
        A = np.abs(V_int)
        total = A.sum(axis=1)
        out: dict[str, np.ndarray] = {}
        for pathway, idx in self.groups.items():
            block = V_int[:, idx]
            mag = np.abs(block)
            activity = mag.sum(axis=1)
            out[f"path_{pathway}_activity"] = activity
            out[f"path_{pathway}_share"] = safe_divide(activity, total, fill=0.0)
            out[f"path_{pathway}_engagement"] = (mag > TOL.flux_zero).mean(axis=1)
            if len(idx) >= self.MIN_REACTIONS:
                pos = positive_part(block).sum(axis=1)
                neg = negative_part(block).sum(axis=1)
                out[f"path_{pathway}_directionality"] = safe_divide(pos - neg, activity, fill=0.0)
        return pd.DataFrame(out)


# --------------------------------------------------------------------------
# F4-F6 - biochemical phenotype, turnover and efficiency
# --------------------------------------------------------------------------

#: Reactions that define erythrocyte redox and energy handling. Each entry maps
#: a label to the RBC-GEM reaction identifier; missing identifiers are skipped
#: and reported, so the family degrades gracefully on other models.
KEY_REACTIONS = {
    "g6pd": "G6PDH2",             # glucose-6-phosphate dehydrogenase
    "pgl": "PGL",                 # 6-phosphogluconolactonase
    "gnd": "GND",                 # 6-phosphogluconate dehydrogenase
    "catalase": "CAT",            # H2O2 -> H2O + O2
    "gpx": "GTHP",                # glutathione peroxidase
    "sod": "SPODM",               # superoxide dismutase
    "hexokinase": "HEX1",
    "pfk": "PFK",
    "pyruvate_kinase": "PYK",
    "ldh": "LDH_L",
    "dpg_mutase": "DPGM",
    "dpg_phosphatase": "DPGase",
    "adenylate_kinase": "ADK1",
    "na_k_atpase": "NaKt",
}

#: Efficiency ratios, expressed as (label, numerator, denominator). Numerators
#: and denominators are feature names produced earlier in this family, so every
#: ratio is traceable to explicit reaction sets rather than being an arbitrary
#: quotient of two convenient numbers.
EFFICIENCY_RATIOS = [
    ("nadph_per_glucose", "turnover_nadph_c", "turnover_glc__D_c"),
    ("nadph_per_atp", "turnover_nadph_c", "turnover_atp_c"),
    ("atp_per_glucose", "turnover_atp_c", "turnover_glc__D_c"),
    ("nadph_per_h2o2", "turnover_nadph_c", "turnover_h2o2_c"),
    ("gsh_per_h2o2", "turnover_gthrd_c", "turnover_h2o2_c"),
    ("redox_buffer_ratio", "turnover_gthrd_c", "turnover_gthox_c"),
    ("nadph_per_nadh", "turnover_nadph_c", "turnover_nadh_c"),
    ("oxidative_ppp_per_glycolysis", "path_pentose_phosphate_activity", "path_glycolysis_activity"),
    ("ppp_flux_split", "rxn_g6pd", "rxn_hexokinase"),
    ("peroxide_clearance_split", "rxn_catalase", "rxn_gpx"),
]


class BiochemicalFamily(FeatureFamily):
    """Metabolite turnover, key-reaction fluxes, redox budget and efficiencies."""

    name = "biochemical"

    def __init__(self, ctx: FeatureContext):
        super().__init__(ctx)
        gem = ctx.gem
        self.metabolites = [m for m in TURNOVER_METABOLITES if gem.has_met(m)]
        self.missing_metabolites = [m for m in TURNOVER_METABOLITES if not gem.has_met(m)]
        met_rows = [gem.met_index(m) for m in self.metabolites]
        # |S| restricted to the internal block and the metabolites of interest
        self.absS = sp.csr_matrix(np.abs(ctx.S_int[met_rows, :].toarray()))

        self.key_rxns: dict[str, int] = {}
        self.missing_reactions: list[str] = []
        internal_positions = {gem.rxns[j]: p for p, j in enumerate(ctx.internal_idx)}
        for label, rxn in KEY_REACTIONS.items():
            if rxn in internal_positions:
                self.key_rxns[label] = internal_positions[rxn]
            else:
                self.missing_reactions.append(rxn)

        # positions of the oxidative PPP reactions inside the internal block
        self.ox_ppp = [internal_positions[r] for r in OXIDATIVE_PPP_RXNS
                       if r in internal_positions]

    # -- registry -------------------------------------------------------
    def specs(self) -> list[FeatureSpec]:
        out = []
        for met in self.metabolites:
            out.append(FeatureSpec(
                f"turnover_{met}", "turnover", "internal",
                f"gross intracellular turnover of {met}",
                r"T_m=\sum_i|S_{mi}v_i|", "flux", False, model_dependent=True))
        for label in self.key_rxns:
            out.append(FeatureSpec(
                f"rxn_{label}", "biochemical", "internal",
                f"flux through {KEY_REACTIONS[label]} ({label})",
                r"v_{" + KEY_REACTIONS[label] + r"}", "flux", True, model_dependent=True))
        out.append(FeatureSpec(
            "oxidative_ppp_flux", "biochemical", "internal",
            "mean flux through the oxidative branch of the pentose phosphate pathway",
            r"\tfrac13(v_{G6PDH2}+v_{PGL}+v_{GND})", "flux", True, model_dependent=True))
        out.append(FeatureSpec(
            "nadph_from_oxidative_ppp", "biochemical", "internal",
            "share of NADPH turnover attributable to the oxidative PPP",
            r"\frac{\sum_{i\in\mathrm{oxPPP}}|S_{\mathrm{nadph},i}v_i|}{T_{\mathrm{nadph}}}",
            "ratio", False, model_dependent=True))
        out.append(FeatureSpec(
            "redox_burden", "biochemical", "internal",
            "peroxide turnover relative to reducing-equivalent turnover",
            r"T_{\mathrm{h2o2}}/(T_{\mathrm{nadph}}+T_{\mathrm{gthrd}})",
            "dimensionless", False, model_dependent=True))
        for label, num, den in EFFICIENCY_RATIOS:
            if not self._ratio_available(num, den):
                continue
            out.append(FeatureSpec(
                f"eff_{label}", "efficiency", "internal",
                f"{num} divided by {den}",
                rf"\mathrm{{{num}}}/\mathrm{{{den}}}", "dimensionless",
                False, model_dependent=True))
        return out

    def _ratio_available(self, num: str, den: str) -> bool:
        known = (
            {f"turnover_{m}" for m in self.metabolites}
            | {f"rxn_{k}" for k in self.key_rxns}
            | {f"path_{p}_activity" for p in self.ctx.pathway_groups("internal")}
        )
        return num in known and den in known

    # -- computation ----------------------------------------------------
    def compute(self, V_int: np.ndarray, V_ext: np.ndarray) -> pd.DataFrame:
        out: dict[str, np.ndarray] = {}

        turnover = np.abs(V_int) @ self.absS.T.toarray()   # (n_samples, n_metabolites)
        for k, met in enumerate(self.metabolites):
            out[f"turnover_{met}"] = turnover[:, k]

        for label, pos in self.key_rxns.items():
            out[f"rxn_{label}"] = V_int[:, pos]

        if self.ox_ppp:
            out["oxidative_ppp_flux"] = V_int[:, self.ox_ppp].mean(axis=1)

        # NADPH sourced from the oxidative PPP, as a share of NADPH turnover
        if "nadph_c" in self.metabolites and self.ox_ppp:
            row = self.metabolites.index("nadph_c")
            coeffs = np.abs(self.absS[row].toarray()).ravel()
            contribution = (np.abs(V_int[:, self.ox_ppp]) * coeffs[self.ox_ppp]).sum(axis=1)
            out["nadph_from_oxidative_ppp"] = safe_divide(
                contribution, turnover[:, row], fill=np.nan
            )

        if all(m in self.metabolites for m in ("h2o2_c", "nadph_c", "gthrd_c")):
            i_h = self.metabolites.index("h2o2_c")
            i_n = self.metabolites.index("nadph_c")
            i_g = self.metabolites.index("gthrd_c")
            out["redox_burden"] = safe_divide(
                turnover[:, i_h], turnover[:, i_n] + turnover[:, i_g], fill=np.nan
            )

        frame = pd.DataFrame(out)
        return frame

    def add_efficiencies(self, features: pd.DataFrame) -> pd.DataFrame:
        """Append efficiency ratios once pathway features are also available."""
        extra = {}
        for label, num, den in EFFICIENCY_RATIOS:
            if num not in features.columns or den not in features.columns:
                continue
            extra[f"eff_{label}"] = safe_divide(
                features[num].to_numpy(), features[den].to_numpy(), fill=np.nan
            )
        return pd.DataFrame(extra, index=features.index)


# --------------------------------------------------------------------------
# F7 - flux-weighted network descriptors
# --------------------------------------------------------------------------

class NetworkFamily(FeatureFamily):
    """Descriptors of how metabolic activity is spread across the network.

    Purely topological quantities are identical for every sample and therefore
    carry no phenotypic information; only *flux-weighted* versions are
    computed. Currency metabolites are excluded, because their promiscuous
    participation would otherwise dominate every weighting.
    """

    name = "network"

    def __init__(self, ctx: FeatureContext):
        super().__init__(ctx)
        gem = ctx.gem
        keep = np.array([m not in CURRENCY_METABOLITES for m in gem.mets])
        S_int = ctx.S_int
        self.absS = sp.csr_matrix(np.abs(S_int[keep, :].toarray()))
        # structural degree of each internal reaction (metabolites touched)
        self.reaction_degree = np.asarray((S_int != 0).sum(axis=0)).ravel().astype(float)
        # structural degree of each retained metabolite (reactions touching it)
        self.metabolite_degree = np.asarray((self.absS != 0).sum(axis=1)).ravel().astype(float)

    def specs(self) -> list[FeatureSpec]:
        return [
            FeatureSpec("net_weighted_reaction_degree", "network", "internal",
                        "flux-weighted mean number of metabolites per active reaction",
                        r"\sum_i p_i d_i,\ p_i=|v_i|/\lVert v\rVert_1", "count", False),
            FeatureSpec("net_metabolite_activity_gini", "network", "internal",
                        "inequality of metabolite turnover across the network",
                        r"G(T),\ T_m=\sum_i|S_{mi}v_i|", "dimensionless", False),
            FeatureSpec("net_metabolite_activity_entropy", "network", "internal",
                        "entropy of the metabolite turnover distribution",
                        r"-\sum_m q_m\ln q_m,\ q_m=T_m/\sum_{m'}T_{m'}", "nats", False),
            FeatureSpec("net_active_metabolites", "network", "internal",
                        "number of metabolites with non-zero turnover",
                        r"\#\{m:T_m>\tau\}", "count", False),
            FeatureSpec("net_hub_turnover_share", "network", "internal",
                        "share of total turnover carried by the ten most connected metabolites",
                        r"\sum_{m\in H_{10}}T_m/\sum_m T_m", "ratio", False),
            FeatureSpec("net_weighted_metabolite_degree", "network", "internal",
                        "turnover-weighted mean metabolite connectivity",
                        r"\sum_m q_m k_m", "count", False),
        ]

    def compute(self, V_int: np.ndarray, V_ext: np.ndarray) -> pd.DataFrame:
        A = np.abs(V_int)
        total = A.sum(axis=1, keepdims=True)
        P = safe_divide(A, total, fill=0.0)
        T = A @ self.absS.T.toarray()      # (n_samples, n_metabolites)
        T_total = T.sum(axis=1, keepdims=True)
        Q = safe_divide(T, T_total, fill=0.0)

        hubs = np.argsort(-self.metabolite_degree)[:10]
        with np.errstate(divide="ignore", invalid="ignore"):
            entropy = -np.where(Q > 0, Q * np.log(Q), 0.0).sum(axis=1)

        return pd.DataFrame(
            {
                "net_weighted_reaction_degree": P @ self.reaction_degree,
                "net_metabolite_activity_gini": gini(T),
                "net_metabolite_activity_entropy": entropy,
                "net_active_metabolites": (T > TOL.flux_zero).sum(axis=1).astype(float),
                "net_hub_turnover_share": safe_divide(
                    T[:, hubs].sum(axis=1), T_total.ravel(), fill=0.0
                ),
                "net_weighted_metabolite_degree": Q @ self.metabolite_degree,
            }
        )


import numpy as np
import pandas as pd


#: Boundary metabolites reported individually because they carry the biology of
#: erythrocyte energy metabolism and of G6PD deficiency. Everything else
#: contributes only through the aggregate role-level features, which keeps the
#: external phenotype interpretable instead of 480-dimensional.
TRACKED_EXCHANGES = {
    "glc__D_e": "glucose",
    "lac__L_e": "lactate",
    "pyr_e": "pyruvate",
    "o2_e": "oxygen",
    "co2_e": "carbon_dioxide",
    "h2o2_e": "hydrogen_peroxide",
    "o2s_e": "superoxide",
    "gthrd_e": "reduced_glutathione",
    "gthox_e": "oxidised_glutathione",
    "cys__L_e": "cysteine",
    "gly_e": "glycine",
    "glu__L_e": "glutamate",
    "ascb__L_e": "ascorbate",
    "dhdascb_e": "dehydroascorbate",
    "adn_e": "adenosine",
    "ade_e": "adenine",
    "ins_e": "inosine",
    "rib__D_e": "ribose",
    "pi_e": "phosphate",
    "nac_e": "nicotinate",
    "ncam_e": "nicotinamide",
    "ribflv_e": "riboflavin",
    "thm_e": "thiamine",
    "h_e": "proton",
    "urea_e": "urea",
    "nh4_e": "ammonium",
}

_ROLES = ("nutrient", "waste", "toxicant", "therapeutic", "ion", "xenobiotic", "metabolite")


class ExternalFamily(FeatureFamily):
    """Nutrient, waste, oxidative, redox and boundary-efficiency features."""

    name = "external"

    def __init__(self, ctx: FeatureContext):
        super().__init__(ctx)
        table = ctx.external_reactions()
        self.table = table
        self.met_ids = table["boundary_metabolite"].to_numpy()
        self.roles = table["boundary_role"].to_numpy()
        self.carbon = table["carbon_count"].to_numpy(dtype=float)
        self.carbon = np.nan_to_num(self.carbon, nan=0.0)
        self.toxicity = np.nan_to_num(table["toxicity_weight"].to_numpy(dtype=float), nan=0.4)

        self.tracked: dict[str, int] = {}
        for position, met in enumerate(self.met_ids):
            if met in TRACKED_EXCHANGES:
                self.tracked[TRACKED_EXCHANGES[met]] = position
        self.role_positions = {
            role: np.flatnonzero(self.roles == role) for role in _ROLES
        }

    # -- registry -------------------------------------------------------
    def specs(self) -> list[FeatureSpec]:
        out: list[FeatureSpec] = []
        for label in self.tracked:
            out.append(FeatureSpec(
                f"ex_{label}_net", "boundary_metabolite", "external",
                f"net exchange flux of {label} (negative = uptake)",
                r"v_{\mathrm{EX}}", "flux", True, model_dependent=True))
            out.append(FeatureSpec(
                f"ex_{label}_uptake", "boundary_metabolite", "external",
                f"uptake flux of {label}, reported as a positive magnitude",
                r"\max(-v_{\mathrm{EX}},0)", "flux", True, model_dependent=True))
        for role in _ROLES:
            if self.role_positions[role].size == 0:
                continue
            out.append(FeatureSpec(
                f"exrole_{role}_uptake", "boundary_role", "external",
                f"total uptake across boundary metabolites classified as {role}",
                r"\sum_{i\in R}\max(-v_i,0)", "flux", True))
            out.append(FeatureSpec(
                f"exrole_{role}_secretion", "boundary_role", "external",
                f"total secretion across boundary metabolites classified as {role}",
                r"\sum_{i\in R}\max(v_i,0)", "flux", True))
        out += [
            FeatureSpec("ex_total_uptake", "exchange_burden", "external",
                        "total uptake flux across all boundary reactions",
                        r"U=\sum_i\max(-v_i,0)", "flux", True),
            FeatureSpec("ex_total_secretion", "exchange_burden", "external",
                        "total secretion flux across all boundary reactions",
                        r"E=\sum_i\max(v_i,0)", "flux", True),
            FeatureSpec("ex_net_balance", "exchange_burden", "external",
                        "net mass flow across the boundary",
                        r"E-U", "flux", True),
            FeatureSpec("ex_turnover_ratio", "exchange_burden", "external",
                        "secretion relative to uptake",
                        r"E/U", "dimensionless", True),
            FeatureSpec("ex_n_uptaken", "exchange_burden", "external",
                        "number of metabolites taken up",
                        r"\#\{i:v_i<-\tau\}", "count", True),
            FeatureSpec("ex_n_secreted", "exchange_burden", "external",
                        "number of metabolites secreted",
                        r"\#\{i:v_i>\tau\}", "count", True),
            FeatureSpec("ex_uptake_entropy", "exchange_burden", "external",
                        "entropy of the uptake distribution across metabolites",
                        r"-\sum_i u_i\ln u_i,\ u_i=\max(-v_i,0)/U", "nats", True),
            FeatureSpec("ex_carbon_uptake", "boundary_efficiency", "external",
                        "carbon atoms entering the system per unit time",
                        r"C_U=\sum_i c_i\max(-v_i,0)", "flux", True,
                        model_dependent=True),
            FeatureSpec("ex_carbon_secretion", "boundary_efficiency", "external",
                        "carbon atoms leaving the system per unit time",
                        r"C_E=\sum_i c_i\max(v_i,0)", "flux", True,
                        model_dependent=True),
            FeatureSpec("ex_carbon_retention", "boundary_efficiency", "external",
                        "fraction of incoming carbon not immediately re-exported",
                        r"(C_U-C_E)/C_U", "ratio", True, model_dependent=True),
            FeatureSpec("ex_oxidative_exposure", "oxidative_environment", "external",
                        "combined uptake of peroxide, superoxide and oxygen",
                        r"\sum_{i\in\mathrm{ROS}}\max(-v_i,0)", "flux", True,
                        model_dependent=True),
            FeatureSpec("ex_reductant_supply", "redox_environment", "external",
                        "combined uptake of glutathione, cysteine and ascorbate",
                        r"\sum_{i\in\mathrm{RED}}\max(-v_i,0)", "flux", True,
                        model_dependent=True),
            FeatureSpec("ex_redox_supply_balance", "redox_environment", "external",
                        "reductant supply relative to oxidative exposure",
                        r"\mathrm{RED}/(\mathrm{RED}+\mathrm{ROS})", "ratio", True,
                        model_dependent=True),
            FeatureSpec("ex_toxicity_load", "admissibility", "external",
                        "toxicity-weighted uptake, used only as a descriptor",
                        r"\sum_i w_i\max(-v_i,0)", "flux", True, model_dependent=True),
        ]
        return out

    # -- computation ----------------------------------------------------
    def compute(self, V_int: np.ndarray, V_ext: np.ndarray) -> pd.DataFrame:
        uptake = negative_part(V_ext)     # positive magnitudes
        secretion = positive_part(V_ext)
        out: dict[str, np.ndarray] = {}

        for label, position in self.tracked.items():
            out[f"ex_{label}_net"] = V_ext[:, position]
            out[f"ex_{label}_uptake"] = uptake[:, position]

        for role, positions in self.role_positions.items():
            if positions.size == 0:
                continue
            out[f"exrole_{role}_uptake"] = uptake[:, positions].sum(axis=1)
            out[f"exrole_{role}_secretion"] = secretion[:, positions].sum(axis=1)

        U = uptake.sum(axis=1)
        E = secretion.sum(axis=1)
        out["ex_total_uptake"] = U
        out["ex_total_secretion"] = E
        out["ex_net_balance"] = E - U
        out["ex_turnover_ratio"] = safe_divide(E, U, fill=np.nan)
        out["ex_n_uptaken"] = (V_ext < -TOL.flux_zero).sum(axis=1).astype(float)
        out["ex_n_secreted"] = (V_ext > TOL.flux_zero).sum(axis=1).astype(float)
        out["ex_uptake_entropy"] = shannon_entropy(uptake)

        C_U = uptake @ self.carbon
        C_E = secretion @ self.carbon
        out["ex_carbon_uptake"] = C_U
        out["ex_carbon_secretion"] = C_E
        out["ex_carbon_retention"] = safe_divide(C_U - C_E, C_U, fill=np.nan)

        ros = [self.tracked[k] for k in ("hydrogen_peroxide", "superoxide", "oxygen")
               if k in self.tracked]
        red = [self.tracked[k] for k in ("reduced_glutathione", "cysteine", "ascorbate")
               if k in self.tracked]
        # absent exchanges make the descriptor missing, never zero
        ros_load = uptake[:, ros].sum(axis=1) if ros else np.full(len(V_ext), np.nan)
        red_load = uptake[:, red].sum(axis=1) if red else np.full(len(V_ext), np.nan)
        out["ex_oxidative_exposure"] = ros_load
        out["ex_reductant_supply"] = red_load
        out["ex_redox_supply_balance"] = safe_divide(
            red_load, red_load + ros_load, fill=np.nan
        )
        out["ex_toxicity_load"] = uptake @ self.toxicity

        return pd.DataFrame(out)


"""Three families the first version of the framework did not have.

They close the gaps the design review named: cofactor and redox allocation,
the relation between a state and the model's own capacity limits, and
descriptors that live *between* the internal and boundary blocks.

Two rules are enforced here and nowhere else in the code:

* a feature whose reactions or metabolites are missing from the model returns
  ``NaN``, never ``0``. A zero would be read downstream as a real biological
  value and would make two models look different for the wrong reason;
* no feature is a function of the reaction whose bound defines the
  perturbation. The saturation of G6PDH2 would recover the severity axis by
  construction, which is why the capacity family keeps only aggregate
  statistics over the whole vector.
"""
import numpy as np
import pandas as pd


#: Cofactor pairs whose turnover allocation the redox family reports.
COFACTOR_PAIRS = {
    "nadph_nadp": ("nadph_c", "nadp_c"),
    "nadh_nad": ("nadh_c", "nad_c"),
    "glutathione": ("gthrd_c", "gthox_c"),
}

#: Adenylates, in the order the energy-charge proxy expects.
ADENYLATES = ("atp_c", "adp_c", "amp_c")


def _metabolite_turnover(S, V_block, idx) -> np.ndarray:
    """Gross stoichiometric turnover of one metabolite over a block of states."""
    row = S[idx, :]
    return np.abs(V_block) @ np.abs(np.asarray(row.todense()).ravel())


class _TurnoverMixin:
    def _index(self, met: str):
        try:
            return int(self.ctx.gem.mets.index(met))
        except (ValueError, AttributeError):
            mets = list(self.ctx.gem.mets)
            return int(mets.index(met)) if met in mets else None

    def _turnover(self, V_int, met: str):
        idx = self._index(met)
        if idx is None:
            return None
        return _metabolite_turnover(self.ctx.S_int, V_int, idx)


class CofactorFamily(_TurnoverMixin, FeatureFamily):
    """How a state allocates reducing equivalents and adenylate turnover.

    Every feature is a *share*: a turnover divided by the turnover of the pool
    it belongs to. Shares are used rather than absolute turnovers because the
    sampled polytope is not in physiological units, so only relative
    allocation is comparable between states and between models.
    """

    name = "cofactor"

    def specs(self) -> list[FeatureSpec]:
        out = []
        for key, (reduced, oxidised) in COFACTOR_PAIRS.items():
            out.append(FeatureSpec(
                name=f"cof_{key}_share", family=self.name, block="internal",
                description=f"turnover share of {reduced} within the {reduced}/{oxidised} pool",
                formula=rf"T_{{{reduced}}} / (T_{{{reduced}}} + T_{{{oxidised}}})",
                units="ratio", signed=False, model_dependent=True,
            ))
        out.append(FeatureSpec(
            name="cof_nadph_over_nadh", family=self.name, block="internal",
            description="share of reducing-equivalent turnover carried by NADPH rather than NADH",
            formula=r"T_{nadph} / (T_{nadph} + T_{nadh})",
            units="ratio", signed=False, model_dependent=True,
        ))
        out.append(FeatureSpec(
            name="cof_energy_charge_proxy", family=self.name, block="internal",
            description=("adenylate energy charge computed on turnover rather than on "
                         "concentration; a flux-level proxy, not a measured charge"),
            formula=r"(T_{atp} + 0.5\,T_{adp}) / (T_{atp} + T_{adp} + T_{amp})",
            units="ratio", signed=False, model_dependent=True,
        ))
        out.append(FeatureSpec(
            name="cof_redox_per_energy", family=self.name, block="internal",
            description="redox turnover relative to adenylate turnover",
            formula=r"(T_{nadph} + T_{gthrd}) / (T_{atp} + T_{adp} + T_{amp})",
            units="ratio", signed=False, model_dependent=True,
        ))
        return out

    def compute(self, V_int: np.ndarray, V_ext: np.ndarray) -> pd.DataFrame:
        n = V_int.shape[0]
        nan = np.full(n, np.nan)
        turnover = {}
        for met in {m for pair in COFACTOR_PAIRS.values() for m in pair} | set(ADENYLATES):
            turnover[met] = self._turnover(V_int, met)

        def share(a, b):
            ta, tb = turnover.get(a), turnover.get(b)
            if ta is None or tb is None:
                return nan
            total = ta + tb
            return np.where(total > 0, ta / np.where(total > 0, total, 1.0), np.nan)

        columns = {f"cof_{key}_share": share(red, ox) for key, (red, ox) in COFACTOR_PAIRS.items()}
        columns["cof_nadph_over_nadh"] = share("nadph_c", "nadh_c")

        atp, adp, amp = (turnover.get(m) for m in ADENYLATES)
        if atp is None or adp is None or amp is None:
            columns["cof_energy_charge_proxy"] = nan
            columns["cof_redox_per_energy"] = nan
        else:
            pool = atp + adp + amp
            columns["cof_energy_charge_proxy"] = np.where(
                pool > 0, (atp + 0.5 * adp) / np.where(pool > 0, pool, 1.0), np.nan)
            nadph, gthrd = turnover.get("nadph_c"), turnover.get("gthrd_c")
            if nadph is None or gthrd is None:
                columns["cof_redox_per_energy"] = nan
            else:
                columns["cof_redox_per_energy"] = np.where(
                    pool > 0, (nadph + gthrd) / np.where(pool > 0, pool, 1.0), np.nan)
        return pd.DataFrame(columns)


class CapacityFamily(FeatureFamily):
    """Where a state sits relative to the bounds the model imposes.

    These are the only features that depend on the constraint set rather than
    on the flux distribution alone, which is exactly why they are reported
    separately: a difference between two models in this family is a statement
    about the models, not about metabolism.
    """

    name = "capacity"

    def __init__(self, ctx):
        super().__init__(ctx)
        lb, ub = np.asarray(ctx.gem.lb, float), np.asarray(ctx.gem.ub, float)
        self.lb_int, self.ub_int = lb[ctx.internal_idx], ub[ctx.internal_idx]
        self.lb_ext, self.ub_ext = lb[ctx.external_idx], ub[ctx.external_idx]
        self._tol = 1e-6

    def specs(self) -> list[FeatureSpec]:
        common = dict(family=self.name, units="ratio", signed=True, model_dependent=True)
        return [
            FeatureSpec(name="cap_at_upper_internal", block="internal",
                        description="fraction of internal reactions within tolerance of their upper bound",
                        formula=r"|\{j : v_j \geq u_j - \epsilon\}| / n_{int}", **common),
            FeatureSpec(name="cap_at_lower_internal", block="internal",
                        description="fraction of internal reactions within tolerance of their lower bound",
                        formula=r"|\{j : v_j \leq l_j + \epsilon\}| / n_{int}", **common),
            FeatureSpec(name="cap_mean_saturation_internal", block="internal",
                        description="mean position within the admissible interval, over bounded internal reactions",
                        formula=r"\mathrm{mean}_j\ (v_j - l_j) / (u_j - l_j)", **common),
            FeatureSpec(name="cap_at_bound_boundary", block="external",
                        description="fraction of exchange reactions at either of their bounds",
                        formula=r"|\{j : v_j \leq l_j + \epsilon \lor v_j \geq u_j - \epsilon\}| / n_{bnd}",
                        **common),
            FeatureSpec(name="cap_mean_saturation_boundary", block="external",
                        description="mean position within the admissible interval, over bounded exchanges",
                        formula=r"\mathrm{mean}_j\ (v_j - l_j) / (u_j - l_j)", **common),
        ]

    @staticmethod
    def _saturation(V, lb, ub):
        span = ub - lb
        usable = np.isfinite(span) & (span > 0)
        if not usable.any():
            return np.full(V.shape[0], np.nan)
        position = (V[:, usable] - lb[usable]) / span[usable]
        return position.mean(axis=1)

    def compute(self, V_int: np.ndarray, V_ext: np.ndarray) -> pd.DataFrame:
        tol = self._tol
        finite_i = np.isfinite(self.ub_int) & np.isfinite(self.lb_int)
        finite_e = np.isfinite(self.ub_ext) & np.isfinite(self.lb_ext)
        n_i = max(int(finite_i.sum()), 1)
        n_e = max(int(finite_e.sum()), 1)
        at_upper = (V_int[:, finite_i] >= self.ub_int[finite_i] - tol).sum(axis=1) / n_i
        at_lower = (V_int[:, finite_i] <= self.lb_int[finite_i] + tol).sum(axis=1) / n_i
        at_bound_e = (
            (V_ext[:, finite_e] >= self.ub_ext[finite_e] - tol)
            | (V_ext[:, finite_e] <= self.lb_ext[finite_e] + tol)
        ).sum(axis=1) / n_e
        return pd.DataFrame({
            "cap_at_upper_internal": at_upper,
            "cap_at_lower_internal": at_lower,
            "cap_mean_saturation_internal": self._saturation(V_int, self.lb_int, self.ub_int),
            "cap_at_bound_boundary": at_bound_e,
            "cap_mean_saturation_boundary": self._saturation(V_ext, self.lb_ext, self.ub_ext),
        })


class CouplingFamily(FeatureFamily):
    """Descriptors that only exist because the two blocks are kept apart.

    They are state-level quantities, one number per sampled state, and they are
    not the same thing as the statistical coupling between blocks measured
    later by canonical correlation: these say how much of a state's activity
    crosses the boundary, those say how predictable one block is from the other.
    """

    name = "coupling"

    def __init__(self, ctx):
        super().__init__(ctx)
        groups = ctx.pathway_groups("internal")
        self.ppp_idx = groups.get("pentose_phosphate", np.array([], dtype=int))
        self.glycolysis_idx = groups.get("glycolysis", np.array([], dtype=int))
        external = ctx.external_reactions()
        self.glucose_idx = self._find(external, "EX_glc__D_e")
        self.oxygen_idx = self._find(external, "EX_o2_e")

    @staticmethod
    def _find(table: pd.DataFrame, rxn_id: str):
        hit = np.flatnonzero((table["rxn_id"] == rxn_id).to_numpy())
        return int(hit[0]) if hit.size else None

    def specs(self) -> list[FeatureSpec]:
        common = dict(family=self.name, block="combined", units="ratio",
                      signed=False, model_dependent=True)
        return [
            FeatureSpec(name="cpl_boundary_share", description=(
                "share of the total flux magnitude of a state that crosses the system boundary"),
                formula=r"\|v_{bnd}\|_1 / \|v\|_1", **common),
            FeatureSpec(name="cpl_secretion_over_uptake", description=(
                "secreted magnitude relative to taken-up magnitude at the boundary"),
                formula=r"\sum_j \max(v_j, 0) / \sum_j \max(-v_j, 0)", **common),
            FeatureSpec(name="cpl_internal_per_boundary", description=(
                "internal activity supported per unit of boundary exchange"),
                formula=r"\|v_{int}\|_1 / \|v_{bnd}\|_1", **common),
            FeatureSpec(name="cpl_ppp_per_glucose_uptake", description=(
                "pentose-phosphate activity per unit of glucose taken up; NaN when the "
                "model has no glucose exchange"),
                formula=r"\sum_{i \in PPP} |v_i| / |v_{EX\_glc}|", **common),
            FeatureSpec(name="cpl_ppp_per_oxygen_uptake", description=(
                "pentose-phosphate activity per unit of oxygen taken up; NaN when the "
                "model has no oxygen exchange"),
                formula=r"\sum_{i \in PPP} |v_i| / |v_{EX\_o2}|", **common),
        ]

    def compute(self, V_int: np.ndarray, V_ext: np.ndarray) -> pd.DataFrame:
        n = V_int.shape[0]
        nan = np.full(n, np.nan)
        internal_l1 = np.abs(V_int).sum(axis=1)
        boundary_l1 = np.abs(V_ext).sum(axis=1)
        total = internal_l1 + boundary_l1
        uptake = np.maximum(-V_ext, 0).sum(axis=1)
        secretion = np.maximum(V_ext, 0).sum(axis=1)
        ppp = np.abs(V_int[:, self.ppp_idx]).sum(axis=1) if self.ppp_idx.size else nan

        def ratio(numerator, denominator):
            if numerator is nan or denominator is None:
                return nan
            safe = np.where(np.abs(denominator) > 0, denominator, np.nan)
            return numerator / safe

        glucose = np.abs(V_ext[:, self.glucose_idx]) if self.glucose_idx is not None else None
        oxygen = np.abs(V_ext[:, self.oxygen_idx]) if self.oxygen_idx is not None else None

        return pd.DataFrame({
            "cpl_boundary_share": ratio(boundary_l1, total),
            "cpl_secretion_over_uptake": ratio(secretion, uptake),
            "cpl_internal_per_boundary": ratio(internal_l1, boundary_l1),
            "cpl_ppp_per_glucose_uptake": ratio(ppp, glucose),
            "cpl_ppp_per_oxygen_uptake": ratio(ppp, oxygen),
        })


r"""The cofactor and carbon budget: where a currency metabolite comes from and where it goes.

Turnover adds production and consumption together. A budget keeps them apart
and asks which reactions make a metabolite and which use it, as shares of its
total production or consumption over the internal block of one state:

$$P_m = \sum_i \max(S_{mi} v_i, 0), \qquad C_m = \sum_i \max(-S_{mi} v_i, 0),
\qquad \text{share}_{m,R} = \frac{\sum_{i \in R} \max(\pm S_{mi} v_i, 0)}{P_m \text{ or } C_m}.$$

A reaction counts as a producer or a consumer according to the sign of its
flux in that state, so a reaction running backwards is booked on the other side.
The sources and sinks are named reactions (glycolytic ATP from PGK and PYK, NADPH
from the oxidative branch, NADH from GAPDH) or a stoichiometric signature (a
glutathione reductase is any reaction that consumes GSSG and a nicotinamide
cofactor and produces GSH), so that the same descriptor can be resolved in a
reconstruction with other reaction names. When a named reaction or metabolite
is absent the descriptor is missing, never zero.

Residence times and half-lives are not descriptors: with a fixed pool size
they are the reciprocal of production and carry no information that the
production does not (they are reported as a supplementary table in stage 21).
"""

#: name -> (metabolite, side, numerator, reading). ``side`` is "production" or
#: "consumption"; the numerator is a tuple of reaction identifiers or a signature
#: key resolved by ``BudgetFamily._signature``. The denominator is the metabolite's
#: total on that side, except where a second entry names another numerator.
BUDGET_DESCRIPTORS = {
    "bud_atp_from_glycolysis": (
        "atp_c", "production", ("PGK", "PYK"),
        "share of ATP production made by the two glycolytic kinases"),
    "bud_atp_to_na_k_pump": (
        "atp_c", "consumption", ("NaKt",),
        "share of ATP consumption spent by the Na+/K+ ATPase"),
    "bud_atp_to_glycolytic_priming": (
        "atp_c", "consumption", ("HEX1", "PFK"),
        "share of ATP consumption invested by hexokinase and phosphofructokinase"),
    "bud_nadph_from_oxidative_ppp": (
        "nadph_c", "production", ("G6PDH2", "GND"),
        "share of NADPH production made by the oxidative pentose-phosphate branch"),
    "bud_nadph_to_glutathione_reduction": (
        "nadph_c", "consumption", "glutathione_reductase",
        "share of NADPH consumption used to reduce glutathione disulfide"),
    "bud_nadh_from_gapdh": (
        "nadh_c", "production", ("GAPD",),
        "share of NADH production made by glyceraldehyde-3-phosphate dehydrogenase"),
    "bud_nadh_to_lactate": (
        "nadh_c", "consumption", ("LDH_L",),
        "share of NADH consumption used to reduce pyruvate to lactate"),
    "bud_nadh_to_cytochrome_b5": (
        "nadh_c", "consumption", ("FCYTOB5OXR",),
        "share of NADH consumption by cytochrome b5 reductase, the methaemoglobin-reducing route"),
    "bud_bpg_through_shunt": (
        "13dpg_c", "consumption", ("DPGM",),
        "share of 1,3-bisphosphoglycerate diverted into the Rapoport–Luebering shunt, "
        "each unit forgoing the ATP of phosphoglycerate kinase"),
    "bud_pyruvate_to_lactate": (
        "pyr_c", "consumption", ("LDH_L",),
        "share of pyruvate consumption reduced to lactate"),
    "bud_co2_from_oxidative_ppp": (
        "co2_c", "production", ("GND",),
        "share of internal CO2 production released by 6-phosphogluconate dehydrogenase, "
        "the carbon of glucose lost through the oxidative branch"),
    "bud_gsh_from_reduction": (
        "gthrd_c", "production", "glutathione_reductase",
        "share of reduced glutathione regenerated from its disulfide rather than made anew"),
}

#: The one budget descriptor that is a ratio of two numerators, not a share.
BUDGET_COVER = ("bud_oxppp_cover_of_glutathione_reduction", "nadph_c",
                "NADPH made by the oxidative branch per NADPH used to reduce glutathione")


class BudgetFamily(FeatureFamily):
    """Production and consumption shares of ATP, NADPH, NADH, 1,3-BPG, pyruvate, CO2 and GSH."""

    name = "budget"

    def __init__(self, ctx):
        super().__init__(ctx)
        mets = list(ctx.gem.mets)
        rxns = list(np.asarray(ctx.gem.rxns)[ctx.internal_idx])
        self._position = {r: k for k, r in enumerate(rxns)}
        self._rows = {m: np.asarray(ctx.S_int[mets.index(m)].todense()).ravel()
                      for m in {spec[0] for spec in BUDGET_DESCRIPTORS.values()} | {"gthox_c"}
                      if m in mets}

    def _signature(self, key: str) -> np.ndarray | None:
        """Internal reactions matching a stoichiometric signature, as positions."""
        if key != "glutathione_reductase":
            raise KeyError(key)
        if not {"gthrd_c", "gthox_c"} <= set(self._rows):
            return None
        cofactor = np.zeros_like(self._rows["gthrd_c"], dtype=bool)
        for met in ("nadph_c", "nadh_c"):
            row = np.asarray(self.ctx.S_int[list(self.ctx.gem.mets).index(met)].todense()).ravel() \
                if met in self.ctx.gem.mets else np.zeros_like(cofactor, dtype=float)
            cofactor |= row < 0
        match = (self._rows["gthox_c"] < 0) & (self._rows["gthrd_c"] > 0) & cofactor
        return np.flatnonzero(match)

    def _members(self, numerator) -> np.ndarray | None:
        if isinstance(numerator, str):
            return self._signature(numerator)
        if not all(r in self._position for r in numerator):
            return None
        return np.array([self._position[r] for r in numerator], dtype=int)

    def specs(self) -> list[FeatureSpec]:
        out = [FeatureSpec(name=name, family=self.name, block="internal", description=reading,
                           formula=(rf"\sum_{{i \in R}} \max({'+' if side == 'production' else '-'}"
                                    rf"S_{{{met},i}} v_i, 0) / {'P' if side == 'production' else 'C'}"
                                    rf"_{{{met}}}"),
                           units="ratio", signed=True, model_dependent=True)
               for name, (met, side, _, reading) in BUDGET_DESCRIPTORS.items()]
        name, met, reading = BUDGET_COVER
        out.append(FeatureSpec(name=name, family=self.name, block="internal", description=reading,
                               formula=r"P_{nadph}^{G6PDH2,GND} / C_{nadph}^{GSSG\ reductases}",
                               units="ratio", signed=True, model_dependent=True))
        return out

    @staticmethod
    def _sides(V: np.ndarray, row: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Per state and reaction: the amount produced and the amount consumed."""
        flow = V * row                                   # S_mi v_i, states x reactions
        return np.maximum(flow, 0.0), np.maximum(-flow, 0.0)

    def compute(self, V_int: np.ndarray, V_ext: np.ndarray) -> pd.DataFrame:
        n = V_int.shape[0]
        out, parts = {}, {}
        sides = {m: self._sides(V_int, row) for m, row in self._rows.items()}
        for name, (met, side, numerator, _) in BUDGET_DESCRIPTORS.items():
            members = self._members(numerator)
            if met not in sides or members is None:
                out[name] = np.full(n, np.nan)
                continue
            produced, consumed = sides[met]
            flows = produced if side == "production" else consumed
            part = flows[:, members].sum(axis=1)
            parts[name] = part
            out[name] = safe_divide(part, flows.sum(axis=1))
        cover, _, _ = BUDGET_COVER
        if {"bud_nadph_from_oxidative_ppp", "bud_nadph_to_glutathione_reduction"} <= set(parts):
            out[cover] = safe_divide(parts["bud_nadph_from_oxidative_ppp"],
                                     parts["bud_nadph_to_glutathione_reduction"])
        else:
            out[cover] = np.full(n, np.nan)
        return pd.DataFrame(out)


def budget_support(feature: str, gem) -> set[str] | None:
    """The reactions a budget descriptor reads: every reaction touching its metabolite."""
    if feature == BUDGET_COVER[0]:
        feature, extra = "bud_nadph_from_oxidative_ppp", {"gthox_c", "gthrd_c"}
    elif feature in BUDGET_DESCRIPTORS:
        extra = set()
    else:
        return None
    met = BUDGET_DESCRIPTORS[feature][0]
    S = gem.S.tocsr()
    mets, rxns = list(gem.mets), np.asarray(gem.rxns)
    out = set()
    for m in {met} | extra:
        if m in mets:
            out |= set(rxns[S[mets.index(m)].indices])
    return out


r"""Compartments: where in the cell a currency metabolite is made, and where the flux runs.

A reconstruction with organelles makes the same metabolite in several places:
ATP in the cytosol by glycolysis and in the mitochondrion by oxidative
phosphorylation, NADH in both, NADPH in both. These descriptors ask how a
state divides that production between the mitochondrion and the cytosol, and
how its internal flux divides between them and the transport that joins them.

For a metabolite $x$ in compartment $K$ (identifier ``x_K`` or ``x[K]``),
gross production is booked by the sign of each flux, as in the budget family, but only
over reactions that make or use $x$ chemically. A reaction that carries $x$
between two compartments moves it and makes nothing, so it is left out:

$$P_{x,K} = \sum_{i \notin T_x} \max(S_{x_K,i}\, v_i, 0), \qquad
\text{share}_{x,K} = \frac{P_{x,K}}{\sum_{K'} P_{x,K'}},$$

with $T_x$ the reactions in which $x$ appears in two compartments and $K'$ every
intracellular compartment. These are gross stoichiometric production shares, not
concentrations, availability or net balances. The flux of a compartment is the share of internal
absolute flux carried by reactions all of whose intracellular metabolites lie in it; the flux
between compartments, the share carried by reactions that join two
intracellular compartments.

Requirements are structural. A share needs the metabolite in at least two
intracellular compartments, one of them the compartment named; the flux shares
need at least two intracellular compartments. A reconstruction without
organelles (RBC-GEM has only the cytosol and the extracellular space) cannot
compute any of them: they are declared, listed as structurally unavailable for
that model, and absent from its table of states, never zero. A share whose
total production is zero in a state is a missing value of that state.
"""
import re

#: The compartments the shares name; any other intracellular compartment still
#: enters the denominator.
COMPARTMENT_NAMES = {"m": "mitochondrion", "c": "cytosol"}
#: Metabolites whose production is divided between compartments.
COMPARTMENT_METABOLITES = ("atp", "nadh", "nadph", "co2")
#: Compartments that are not inside the cell.
EXTRACELLULAR = frozenset({"e", "b"})
_COMPARTMENT_ID = re.compile(r"^(.*?)(?:_([A-Za-z0-9]{1,2})|\[([A-Za-z0-9]{1,2})\])$")


def split_metabolite(identifier: str) -> tuple[str, str]:
    """``atp_m`` or ``atp[m]`` -> (``atp``, ``m``); no recognisable suffix -> (id, '')."""
    match = _COMPARTMENT_ID.match(str(identifier))
    if not match:
        return str(identifier), ""
    return match.group(1), match.group(2) or match.group(3)


def compartment_layout(mets, S_int) -> dict:
    """Where each metabolite lives, and which internal reactions join compartments."""
    parts = [split_metabolite(m) for m in mets]
    base = np.array([p[0] for p in parts], dtype=object)
    comp = np.array([p[1] for p in parts], dtype=object)
    inside = sorted({c for c in comp if c and c not in EXTRACELLULAR})
    A = abs(sp.csc_matrix(S_int)) > 0
    n = A.shape[1]
    single, between = np.full(n, "", dtype=object), np.zeros(n, dtype=bool)
    for j in range(n):
        rows = A.indices[A.indptr[j]:A.indptr[j + 1]]
        places = {comp[r] for r in rows if comp[r] not in EXTRACELLULAR}
        if len(places) == 1:
            single[j] = next(iter(places))
        between[j] = len(places) >= 2
    return {"base": base, "comp": comp, "inside": inside, "single": single,
            "between": between}


def compartment_specs() -> list[tuple[str, str, str, str]]:
    """(name, kind, metabolite or compartment, reading) for every declared descriptor."""
    out = []
    for met in COMPARTMENT_METABOLITES:
        for code, place in COMPARTMENT_NAMES.items():
            out.append((f"cmp_{met}_made_in_{place}", "share", f"{met}:{code}",
                        f"share of the cell's gross {met.upper()} production made in the {place} (stoichiometric, not net)"))
    for code, place in COMPARTMENT_NAMES.items():
        out.append((f"cmp_flux_in_{place}", "flux", code,
                    f"share of internal flux carried by reactions inside the {place}"))
    out.append(("cmp_flux_between_compartments", "flux", "between",
                "share of internal flux carried by reactions joining two intracellular "
                "compartments"))
    return out


def compartment_availability(layout: dict) -> dict[str, tuple[bool, str]]:
    """Per descriptor: can this reconstruction compute it, and if not, why."""
    inside = layout["inside"]
    out = {}
    for name, kind, target, _ in compartment_specs():
        if kind == "share":
            met, code = target.split(":")
            places = sorted({c for b, c in zip(layout["base"], layout["comp"])
                             if b == met and c and c not in EXTRACELLULAR})
            ok = code in places and len(places) >= 2
            why = "" if ok else (f"{met} is in {len(places)} intracellular compartment(s) "
                                 f"({', '.join(places) or 'none'}); needs {code} and one other")
        else:
            ok = len(inside) >= 2 and (target == "between" or target in inside)
            why = "" if ok else (f"{len(inside)} intracellular compartment(s) "
                                 f"({', '.join(inside) or 'none'}); needs two")
        out[name] = (ok, why)
    return out


class CompartmentFamily(FeatureFamily):
    """Mitochondrial and cytosolic shares of ATP, NADH, NADPH and CO2 production, and of flux."""

    name = "compartment"

    def __init__(self, ctx):
        super().__init__(ctx)
        self.layout = compartment_layout(list(ctx.gem.mets), ctx.S_int)
        self.availability = compartment_availability(self.layout)
        self.unavailable = {k for k, (ok, _) in self.availability.items() if not ok}

    def specs(self) -> list[FeatureSpec]:
        return [FeatureSpec(name=name, family=self.name, block="internal", description=reading,
                            formula=(r"P_{x,K} / \sum_{K'} P_{x,K'}" if kind == "share" else
                                     r"\sum_{i \in K} |v_i| / \sum_i |v_i|"),
                            units="ratio", signed=kind == "share", model_dependent=True)
                for name, kind, _, reading in compartment_specs()]

    def _production(self, V: np.ndarray, met: str) -> dict[str, np.ndarray]:
        """Chemical production of ``met`` in every intracellular compartment."""
        rows = [r for r, (b, c) in enumerate(zip(self.layout["base"], self.layout["comp"]))
                if b == met and c and c not in EXTRACELLULAR]
        S = self.ctx.S_int
        moving = np.asarray((abs(S[rows]) > 0).sum(axis=0)).ravel() >= 2       # met in two compartments
        out = {}
        for r in rows:
            row = np.asarray(S[r].todense()).ravel() * ~moving
            out[self.layout["comp"][r]] = np.maximum(V * row, 0.0).sum(axis=1)
        return out

    def compute(self, V_int: np.ndarray, V_ext: np.ndarray) -> pd.DataFrame:
        n = V_int.shape[0]
        out, made = {}, {}
        total = np.abs(V_int).sum(axis=1)
        for name, kind, target, _ in compartment_specs():
            if name in self.unavailable:
                out[name] = np.full(n, np.nan)
                continue
            if kind == "share":
                met, code = target.split(":")
                if met not in made:
                    made[met] = self._production(V_int, met)
                everywhere = sum(made[met].values())
                out[name] = safe_divide(made[met][code], everywhere)
            elif target == "between":
                out[name] = safe_divide(np.abs(V_int[:, self.layout["between"]]).sum(axis=1), total)
            else:
                out[name] = safe_divide(
                    np.abs(V_int[:, self.layout["single"] == target]).sum(axis=1), total)
        return pd.DataFrame(out)


def compartment_support(feature: str, gem) -> set[str] | None:
    """The reactions a compartment descriptor reads."""
    specs = {name: (kind, target) for name, kind, target, _ in compartment_specs()}
    if feature not in specs:
        return None
    kind, target = specs[feature]
    S = gem.S.tocsr()
    rxns = np.asarray(gem.rxns)
    if kind == "share":
        met = target.split(":")[0]
        rows = [r for r, m in enumerate(gem.mets) if split_metabolite(m)[0] == met]
        return set(rxns[np.unique(S[rows].indices)]) if rows else set()
    return set(rxns)


r"""Functional families: elements, amino acids, cofactor sources, oxidative defence, allocation.

These families describe what a state does with its matter and its cofactors,
computed only from the stoichiometric matrix, the metabolite formulas and
identifiers, and the pathway groups of the annotation. Recognising the
metabolites and reactions involved is kept apart from computing the
descriptors, so that the same definitions resolve in other reconstructions.

*Elements.* For an exchanged metabolite $i$ with $n_{i,E}$ atoms of element $E$
in its formula, and $x_i$ its export (positive out of the cell, read from the
sign of its stoichiometric coefficient),

$$U_E = \sum_i n_{i,E}\max(-x_i,0), \qquad X_E = \sum_i n_{i,E}\max(x_i,0),$$

and the share of carbon uptake that each source supplies is
$f_i = n_{i,C}\max(-x_i,0)/U_C$. Its Shannon entropy and its effective number
$1/\sum_i f_i^2$ describe how many sources the state draws carbon from; they
describe source use, not metabolic flexibility, which would need responses to
changing supply.

*Amino acids.* The twenty proteinogenic amino acids are recognised by
identifier (BiGG or VMH spelling, L forms and glycine) and, failing that, by a
normalised name; D forms are excluded. Uptake and secretion are read from
their exchange reactions.

*Cofactor sources and sinks.* For a carrier (ATP, GTP, NADPH, NADH, FADH2,
reduced glutathione, reduced thioredoxin), gross production and consumption
are booked by the sign of each flux over internal reactions, leaving out
reactions that only move the carrier between compartments, and divided among
pathway groups. At steady state the production and the consumption of an
internal carrier are equal, so their balance carries no information; how
many routes produce and consume it does. The entropy and the dominant share
of its sources and sinks say whether a carrier depends on one route or on many.

*Oxidative defence.* Peroxide production and removal are booked the same way,
and removal is divided among catalase, glutathione peroxidase and
peroxiredoxins, each recognised by a stoichiometric signature. These are
model-based removal-to-production ratios and shares, not oxidative stress: a
flux vector contains no concentration of reactive species and no damage.

*Allocation.* Gross ATP consumption, gross NADPH consumption and internal
absolute flux are divided among functions (energy metabolism, redox
maintenance, biosynthesis, transport), each a union of pathway groups.

*Conversion.* ATP and NADPH made per carbon atom taken up, lactate carbon
secreted per glucose carbon taken up, and CO2 carbon released per carbon taken
up relate the interior to the boundary.

Every descriptor whose requirement the reconstruction lacks (no exchanged
metabolite with that element, no recognised amino acid, no carrier, no
signature) is declared, structurally unavailable for that model, and never
zero. A share or ratio whose denominator is zero in a state is missing there.
"""
import re

import numpy as np
import pandas as pd
import scipy.sparse as sp

_ELEMENT = re.compile(r"([A-Z][a-z]?)(\d*(?:\.\d+)?)")


def element_count(formula, element: str) -> float:
    """Atoms of ``element`` in a formula such as ``C6H12O6``; 0 if absent or unreadable."""
    if not isinstance(formula, str) or not formula:
        return 0.0
    total = 0.0
    for symbol, count in _ELEMENT.findall(formula):
        if symbol == element:
            total += float(count) if count else 1.0
    return total


# ---------------------------------------------------------------- amino acids
#: canonical name -> (identifier stems, BiGG and VMH, without compartment)
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
#: names that also appear in the metabolite names of reconstructions
_AA_SYNONYMS = {"aspartic acid": "aspartate", "glutamic acid": "glutamate",
                "l-aspartate": "aspartate", "l-glutamate": "glutamate"}
ESSENTIAL_AMINO_ACIDS = ("histidine", "isoleucine", "leucine", "lysine", "methionine",
                         "phenylalanine", "threonine", "tryptophan", "valine")
#: amino acids the boundary family already reports one by one
_TRACKED_AMINO_ACIDS = ("cysteine", "glycine", "glutamate")


def _normalise_name(name: str) -> str:
    text = str(name).strip().lower()
    text = re.sub(r"\(.*?\)", "", text).strip()
    return text


def recognise_amino_acid(identifier: str, name: str = "") -> str | None:
    """Canonical amino-acid name for a metabolite, or None. D forms are not amino acids here."""
    stem, _ = split_metabolite(identifier)
    for canonical, stems in AMINO_ACIDS.items():
        if stem in stems:
            return canonical
    text = _normalise_name(name)
    if not text or text.startswith("d-"):
        return None
    text = _AA_SYNONYMS.get(text, text)
    if text.startswith("l-"):
        text = text[2:]
    return text if text in AMINO_ACIDS else None


# ------------------------------------------------------------ cofactor carriers
#: carrier label -> identifier stems of its reduced or charged form
CARRIERS = {
    "atp": ("atp",), "gtp": ("gtp",), "nadph": ("nadph",), "nadh": ("nadh",),
    "fadh2": ("fadh2",), "glutathione": ("gthrd",), "thioredoxin": ("trdrd", "trdxrd"),
}

#: function -> pathway groups of the annotation
FUNCTIONS = {
    "energy": ("glycolysis", "rapoport_luebering", "pyruvate_tca", "oxidative_phosphorylation"),
    "redox": ("pentose_phosphate", "glutathione", "redox_detox"),
    "biosynthesis": ("amino_acid", "nucleotide", "purine", "pyrimidine", "lipid", "glycan",
                     "cofactor_vitamin", "porphyrin_hemoglobin"),
    "transport": ("transport",),
}


def _fspec(name, family, block, description, formula, units):
    return FeatureSpec(name=name, family=family, block=block, description=description,
                       formula=formula, units=units, signed=True, model_dependent=True)


class _Resolved:
    """The metabolites, rows and reaction sets every functional family reads, found once."""

    def __init__(self, ctx):
        self.ctx = ctx
        gem = ctx.gem
        self.mets = np.asarray(gem.mets)
        self.stems = np.array([split_metabolite(m)[0] for m in self.mets], dtype=object)
        self.comps = np.array([split_metabolite(m)[1] for m in self.mets], dtype=object)
        self.S_int = sp.csr_matrix(ctx.S_int)
        internal = ctx.internal_reactions()
        self.int_rxns = internal["rxn_id"].to_numpy()
        self.pathway = internal["pathway"].fillna("unassigned").to_numpy()
        # exchanges: the exchanged metabolite and the sign that makes its export positive
        S_ext = sp.csc_matrix(ctx.S_ext)
        external = ctx.external_reactions()
        self.ext_rxns = external["rxn_id"].to_numpy()
        n_ext = S_ext.shape[1]
        self.ext_met = np.full(n_ext, -1)
        self.ext_sign = np.zeros(n_ext)
        for j in range(n_ext):
            rows = S_ext.indices[S_ext.indptr[j]:S_ext.indptr[j + 1]]
            data = S_ext.data[S_ext.indptr[j]:S_ext.indptr[j + 1]]
            if rows.size == 1:
                self.ext_met[j] = rows[0]
                self.ext_sign[j] = -np.sign(data[0])      # x_e <=>: export = v
        formulas = np.asarray(gem.met_formulas, dtype=object)
        names = np.asarray(gem.met_names, dtype=object)
        self.ext_formula = np.array([formulas[m] if m >= 0 else "" for m in self.ext_met], dtype=object)
        self.ext_stem = np.array([self.stems[m] if m >= 0 else "" for m in self.ext_met], dtype=object)
        self.ext_amino = np.array([recognise_amino_acid(self.mets[m], names[m]) if m >= 0 else None
                                   for m in self.ext_met], dtype=object)

    def export(self, V_ext: np.ndarray) -> np.ndarray:
        """Export of each exchanged metabolite, positive out of the cell."""
        return V_ext * self.ext_sign

    def ext_where(self, stems) -> np.ndarray:
        return np.flatnonzero(np.isin(self.ext_stem, list(stems)))

    def carrier_rows(self, stems) -> np.ndarray:
        """Metabolite rows of a carrier, every intracellular compartment."""
        return np.flatnonzero(np.isin(self.stems, list(stems)) & ~np.isin(self.comps, ["e", "b"]))

    def carrier_flow(self, stems) -> np.ndarray | None:
        """Stoichiometric coefficient of a carrier in each internal reaction that makes or uses it.

        Rows of every intracellular compartment are added; a reaction in which the
        carrier appears in two compartments, the extracellular one included, only
        moves it and gets zero.
        """
        rows = self.carrier_rows(stems)
        if rows.size == 0:
            return None
        coefficient = self.S_int[rows].toarray().sum(axis=0)
        # every compartment, the extracellular one included, decides what is a move
        everywhere = self.S_int[np.flatnonzero(np.isin(self.stems, list(stems)))].toarray()
        moves = (everywhere > 0).any(axis=0) & (everywhere < 0).any(axis=0)
        coefficient[moves] = 0.0
        return coefficient

    def reactions_of(self, positions) -> set[str]:
        return set(self.int_rxns[np.asarray(positions, dtype=int)])

    def signature(self, consumes: tuple, produces: tuple = (), only: bool = False,
                  partner_prefix: str = "") -> np.ndarray:
        """Internal reactions that consume (and produce) the given stems, as positions.

        ``only`` asks that nothing else be consumed; ``partner_prefix`` asks that a
        metabolite whose identifier starts with it also be consumed.
        """
        A = sp.csc_matrix(self.S_int)
        hits = []
        for j in range(A.shape[1]):
            rows = A.indices[A.indptr[j]:A.indptr[j + 1]]
            data = A.data[A.indptr[j]:A.indptr[j + 1]]
            used = {self.stems[r] for r, d in zip(rows, data) if d < 0}
            made = {self.stems[r] for r, d in zip(rows, data) if d > 0}
            if not set(consumes) <= used or not set(produces) <= made:
                continue
            if only and used - set(consumes) - {"h"}:
                continue
            if partner_prefix and not any(s.startswith(partner_prefix) for s in used):
                continue
            hits.append(j)
        return np.asarray(hits, dtype=int)


def _resolved(ctx) -> _Resolved:
    if getattr(ctx, "_functional", None) is None:
        ctx._functional = _Resolved(ctx)
    return ctx._functional


def _sides(V: np.ndarray, coefficient: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    flow = V * coefficient
    return np.maximum(flow, 0.0), np.maximum(-flow, 0.0)


def _entropy_and_top(parts: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Shannon entropy, effective number and largest share of rows of non-negative parts."""
    total = parts.sum(axis=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        f = np.where(total > 0, parts / total, np.nan)
        logs = np.where(f > 0, np.log(np.where(f > 0, f, 1.0)), 0.0)
        H = -(np.nan_to_num(f) * logs).sum(axis=1)
        D = 1.0 / np.nansum(f ** 2, axis=1)
    empty = total.ravel() <= 0
    H[empty], D[empty] = np.nan, np.nan
    return H, D, np.where(empty, np.nan, np.nanmax(np.nan_to_num(f), axis=1))


# ================================================================= elements
class ElementalFamily(FeatureFamily):
    """Uptake and secretion of N, S and P; carbon sources, sinks and their diversity."""

    name = "elemental"
    SPECS = [
        *[(f"elem_{e.lower()}_{side}", "elemental",
           f"{label} atoms {'taken up' if side == 'uptake' else 'secreted'} per unit time, from the formulas of exchanged metabolites",
           rf"\sum_i n_{{i,{e}}}\max({'-' if side == 'uptake' else ''}x_i,0)", "flux")
          for e, label in (("N", "nitrogen"), ("S", "sulfur"), ("P", "phosphorus"))
          for side in ("uptake", "secretion")],
        ("elem_c_from_glucose", "elemental", "share of carbon uptake supplied by glucose",
         r"U_C^{glc}/U_C", "ratio"),
        ("elem_c_from_lactate_pyruvate", "elemental", "share of carbon uptake supplied by lactate and pyruvate",
         r"(U_C^{lac}+U_C^{pyr})/U_C", "ratio"),
        ("elem_c_from_amino_acids", "elemental", "share of carbon uptake supplied by amino acids",
         r"U_C^{AA}/U_C", "ratio"),
        ("elem_c_to_lactate", "elemental", "share of carbon secretion leaving as lactate",
         r"X_C^{lac}/X_C", "ratio"),
        ("elem_c_to_co2", "elemental", "share of carbon secretion leaving as CO2",
         r"X_C^{CO_2}/X_C", "ratio"),
        ("elem_n_from_amino_acids", "elemental", "share of nitrogen uptake supplied by amino acids",
         r"U_N^{AA}/U_N", "ratio"),
        ("div_c_source_entropy", "source_diversity", "Shannon entropy of the carbon sources taken up",
         r"-\sum_i f_i\ln f_i,\ f_i=n_{i,C}\max(-x_i,0)/U_C", "nats"),
        ("div_c_effective_sources", "source_diversity", "effective number of carbon sources",
         r"1/\sum_i f_i^2", "dimensionless"),
        ("div_c_dominant_source", "source_diversity", "share of carbon uptake from the largest single source",
         r"\max_i f_i", "ratio"),
        ("div_n_source_entropy", "source_diversity", "Shannon entropy of the nitrogen sources taken up",
         r"-\sum_i g_i\ln g_i,\ g_i=n_{i,N}\max(-x_i,0)/U_N", "nats"),
        ("div_n_effective_sources", "source_diversity", "effective number of nitrogen sources",
         r"1/\sum_i g_i^2", "dimensionless"),
    ]

    def __init__(self, ctx):
        super().__init__(ctx)
        self.r = _resolved(ctx)
        self.atoms = {e: np.array([element_count(f, e) for f in self.r.ext_formula])
                      for e in ("C", "N", "S", "P")}
        self.glucose = self.r.ext_where(("glc__D", "glc_D"))
        self.lactate_pyruvate = self.r.ext_where(("lac__L", "lac_L", "pyr"))
        self.lactate = self.r.ext_where(("lac__L", "lac_L"))
        self.co2 = self.r.ext_where(("co2",))
        self.amino = np.flatnonzero(np.array([a is not None for a in self.r.ext_amino]))
        has = {e: bool((self.atoms[e] > 0).any()) for e in self.atoms}
        need = {
            "elem_n": (has["N"], "no exchanged metabolite contains nitrogen"),
            "elem_s": (has["S"], "no exchanged metabolite contains sulfur"),
            "elem_p": (has["P"], "no exchanged metabolite contains phosphorus"),
            "elem_c_from_glucose": (has["C"] and self.glucose.size > 0, "no glucose exchange"),
            "elem_c_from_lactate_pyruvate": (has["C"] and self.lactate_pyruvate.size > 0,
                                             "no lactate or pyruvate exchange"),
            "elem_c_from_amino_acids": (has["C"] and self.amino.size > 0, "no recognised amino-acid exchange"),
            "elem_c_to_lactate": (has["C"] and self.lactate.size > 0, "no lactate exchange"),
            "elem_c_to_co2": (has["C"] and self.co2.size > 0, "no CO2 exchange"),
            "elem_n_from_amino_acids": (has["N"] and self.amino.size > 0, "no recognised amino-acid exchange"),
            "div_c": (has["C"], "no exchanged metabolite contains carbon"),
            "div_n": (has["N"], "no exchanged metabolite contains nitrogen"),
        }
        self.availability = {}
        for name, *_ in self.SPECS:
            key = next(k for k in need if name.startswith(k))
            ok, why = need[key]
            self.availability[name] = (ok, "" if ok else why)

    def specs(self):
        return [_fspec(n, fam, "external", d, f, u) for n, fam, d, f, u in self.SPECS]

    def compute(self, V_int, V_ext):
        x = self.r.export(V_ext)
        up, out = np.maximum(-x, 0.0), np.maximum(x, 0.0)
        res = {}
        for e in ("N", "S", "P"):
            res[f"elem_{e.lower()}_uptake"] = up @ self.atoms[e]
            res[f"elem_{e.lower()}_secretion"] = out @ self.atoms[e]
        cu, cx = up * self.atoms["C"], out * self.atoms["C"]
        UC, XC = cu.sum(axis=1), cx.sum(axis=1)
        nu = up * self.atoms["N"]
        UN = nu.sum(axis=1)
        res["elem_c_from_glucose"] = safe_divide(cu[:, self.glucose].sum(axis=1), UC)
        res["elem_c_from_lactate_pyruvate"] = safe_divide(cu[:, self.lactate_pyruvate].sum(axis=1), UC)
        res["elem_c_from_amino_acids"] = safe_divide(cu[:, self.amino].sum(axis=1), UC)
        res["elem_c_to_lactate"] = safe_divide(cx[:, self.lactate].sum(axis=1), XC)
        res["elem_c_to_co2"] = safe_divide(cx[:, self.co2].sum(axis=1), XC)
        res["elem_n_from_amino_acids"] = safe_divide(nu[:, self.amino].sum(axis=1), UN)
        H, D, top = _entropy_and_top(cu)
        res["div_c_source_entropy"], res["div_c_effective_sources"], res["div_c_dominant_source"] = H, D, top
        H, D, _ = _entropy_and_top(nu)
        res["div_n_source_entropy"], res["div_n_effective_sources"] = H, D
        return pd.DataFrame(res)

    def support(self, feature: str) -> set[str]:
        r = self.r
        pick = {
            "elem_c_from_glucose": self.atoms["C"] > 0, "elem_c_from_lactate_pyruvate": self.atoms["C"] > 0,
            "elem_c_from_amino_acids": self.atoms["C"] > 0, "elem_c_to_lactate": self.atoms["C"] > 0,
            "elem_c_to_co2": self.atoms["C"] > 0, "elem_n_from_amino_acids": self.atoms["N"] > 0,
        }
        if feature in pick:
            return set(r.ext_rxns[pick[feature]])
        element = {"elem_n": "N", "elem_s": "S", "elem_p": "P", "div_c": "C", "div_n": "N"}
        key = next(k for k in element if feature.startswith(k))
        return set(r.ext_rxns[self.atoms[element[key]] > 0])


# =============================================================== amino acids
class AminoAcidFamily(FeatureFamily):
    """Uptake and secretion of the proteinogenic amino acids, one by one and in aggregate."""

    name = "amino_acid"

    def __init__(self, ctx):
        super().__init__(ctx)
        self.r = _resolved(ctx)
        self.position = {}
        for j, aa in enumerate(self.r.ext_amino):
            if aa is not None:
                self.position.setdefault(aa, []).append(j)
        self.single = [aa for aa in AMINO_ACIDS if aa not in _TRACKED_AMINO_ACIDS]
        found = bool(self.position)
        self.availability = {f"aa_{aa}_net_uptake": (aa in self.position,
                                                     "" if aa in self.position else f"no {aa} exchange")
                             for aa in self.single}
        for name in ("aa_total_uptake", "aa_total_secretion", "aa_uptake_entropy", "aa_n_taken_up"):
            self.availability[name] = (found, "" if found else "no recognised amino-acid exchange")
        essential = any(aa in self.position for aa in ESSENTIAL_AMINO_ACIDS)
        self.availability["aa_essential_uptake"] = (essential, "" if essential else
                                                    "no essential amino-acid exchange")

    def specs(self):
        out = [_fspec(f"aa_{aa}_net_uptake", "amino_acid", "external",
                     f"net uptake of {aa} (uptake minus secretion)", rf"-x_{{{aa}}}", "flux")
               for aa in self.single]
        out += [
            _fspec("aa_total_uptake", "amino_acid", "external", "total uptake of amino acids",
                  r"\sum_{a}\max(-x_a,0)", "flux"),
            _fspec("aa_total_secretion", "amino_acid", "external", "total secretion of amino acids",
                  r"\sum_{a}\max(x_a,0)", "flux"),
            _fspec("aa_essential_uptake", "amino_acid", "external",
                  "uptake of the nine essential amino acids", r"\sum_{a\in E}\max(-x_a,0)", "flux"),
            _fspec("aa_uptake_entropy", "amino_acid", "external",
                  "Shannon entropy of the amino-acid uptake profile",
                  r"-\sum_a u_a\ln u_a,\ u_a=\max(-x_a,0)/\sum_b\max(-x_b,0)", "nats"),
            _fspec("aa_n_taken_up", "amino_acid", "external", "number of amino acids taken up",
                  r"\#\{a:-x_a>\tau\}", "count"),
        ]
        return out

    def compute(self, V_int, V_ext):
        x = self.r.export(V_ext)
        n = V_ext.shape[0]
        per = {aa: x[:, idx].sum(axis=1) for aa, idx in self.position.items()}
        res = {f"aa_{aa}_net_uptake": (-per[aa] if aa in per else np.full(n, np.nan)) for aa in self.single}
        if per:
            X = np.column_stack([per[aa] for aa in per])
            up, out = np.maximum(-X, 0.0), np.maximum(X, 0.0)
            res["aa_total_uptake"] = up.sum(axis=1)
            res["aa_total_secretion"] = out.sum(axis=1)
            essential = [k for k, aa in enumerate(per) if aa in ESSENTIAL_AMINO_ACIDS]
            res["aa_essential_uptake"] = up[:, essential].sum(axis=1) if essential else np.full(n, np.nan)
            res["aa_uptake_entropy"] = _entropy_and_top(up)[0]
            res["aa_n_taken_up"] = (up > TOL.flux_zero).sum(axis=1).astype(float)
        else:
            for name in ("aa_total_uptake", "aa_total_secretion", "aa_essential_uptake",
                         "aa_uptake_entropy", "aa_n_taken_up"):
                res[name] = np.full(n, np.nan)
        return pd.DataFrame(res)

    def support(self, feature: str) -> set[str]:
        for aa in self.single:
            if feature == f"aa_{aa}_net_uptake":
                return set(self.r.ext_rxns[self.position.get(aa, [])])
        if feature == "aa_essential_uptake":
            return {self.r.ext_rxns[j] for aa in ESSENTIAL_AMINO_ACIDS for j in self.position.get(aa, [])}
        return {self.r.ext_rxns[j] for idx in self.position.values() for j in idx}


# ========================================================= cofactor sources
class CofactorSourceFamily(FeatureFamily):
    """How many pathway groups produce and consume each carrier, and how unevenly."""

    name = "cofactor_sources"
    KINDS = (("source_entropy", "Shannon entropy of its production across pathway groups", "nats"),
             ("dominant_source", "share of its production made by the largest pathway group", "ratio"),
             ("sink_entropy", "Shannon entropy of its consumption across pathway groups", "nats"),
             ("dominant_sink", "share of its consumption used by the largest pathway group", "ratio"))

    def __init__(self, ctx):
        super().__init__(ctx)
        self.r = _resolved(ctx)
        self.groups = sorted(set(self.r.pathway))
        member = np.array([[p == g for p in self.r.pathway] for g in self.groups], dtype=float)
        self.membership = member.T                        # reactions x groups
        self.coefficient = {c: self.r.carrier_flow(stems) for c, stems in CARRIERS.items()}
        self.availability = {}
        for carrier in CARRIERS:
            ok = self.coefficient[carrier] is not None and bool(np.any(self.coefficient[carrier]))
            for kind, *_ in self.KINDS:
                self.availability[f"src_{carrier}_{kind}"] = (ok, "" if ok else f"no {carrier} carrier")

    def specs(self):
        out = []
        for carrier in CARRIERS:
            for kind, reading, units in self.KINDS:
                side = "P" if "source" in kind else "C"
                formula = (rf"-\sum_g s_g\ln s_g,\ s_g={side}_{{{carrier},g}}/{side}_{{{carrier}}}"
                           if "entropy" in kind else rf"\max_g {side}_{{{carrier},g}}/{side}_{{{carrier}}}")
                out.append(_fspec(f"src_{carrier}_{kind}", "cofactor_sources", "internal",
                                 f"{carrier}: {reading}", formula, units))
        return out

    def compute(self, V_int, V_ext):
        n = V_int.shape[0]
        res = {}
        for carrier in CARRIERS:
            coefficient = self.coefficient[carrier]
            if coefficient is None or not np.any(coefficient):
                for kind, *_ in self.KINDS:
                    res[f"src_{carrier}_{kind}"] = np.full(n, np.nan)
                continue
            made, used = _sides(V_int, coefficient)
            for prefix, parts in (("source", made @ self.membership), ("sink", used @ self.membership)):
                H, _, top = _entropy_and_top(parts)
                res[f"src_{carrier}_{prefix}_entropy"] = H
                res[f"src_{carrier}_dominant_{prefix}"] = top
        return pd.DataFrame(res)

    def support(self, feature: str) -> set[str]:
        carrier = next(c for c in CARRIERS if feature.startswith(f"src_{c}_"))
        coefficient = self.coefficient[carrier]
        return set() if coefficient is None else self.r.reactions_of(np.flatnonzero(coefficient))


# ========================================================== oxidative defence
class OxidativeDefenseFamily(FeatureFamily):
    """Peroxide and superoxide production and removal, and the routes of removal."""

    name = "oxidative_defense"
    SPECS = [
        ("ros_h2o2_production", "gross internal production of hydrogen peroxide",
         r"P_{H_2O_2}=\sum_i\max(S_{H_2O_2,i}v_i,0)", "flux"),
        ("ros_h2o2_removal", "gross internal removal of hydrogen peroxide",
         r"D_{H_2O_2}=\sum_i\max(-S_{H_2O_2,i}v_i,0)", "flux"),
        ("ros_h2o2_removal_per_production",
         "model-based peroxide removal-to-production ratio (above 1 when peroxide is taken up)",
         r"D_{H_2O_2}/P_{H_2O_2}", "ratio"),
        ("ros_h2o2_removed_by_catalase", "share of peroxide removal by catalase",
         r"D_{H_2O_2}^{CAT}/D_{H_2O_2}", "ratio"),
        ("ros_h2o2_removed_by_glutathione_peroxidase",
         "share of peroxide removal by glutathione peroxidase",
         r"D_{H_2O_2}^{GPX}/D_{H_2O_2}", "ratio"),
        ("ros_h2o2_removed_by_peroxiredoxin", "share of peroxide removal by peroxiredoxins",
         r"D_{H_2O_2}^{PRDX}/D_{H_2O_2}", "ratio"),
        ("ros_superoxide_production", "gross internal production of superoxide",
         r"P_{O_2^-}=\sum_i\max(S_{O_2^-,i}v_i,0)", "flux"),
        ("ros_superoxide_dismutated", "share of superoxide removal by superoxide dismutase",
         r"D_{O_2^-}^{SOD}/D_{O_2^-}", "ratio"),
        ("ros_glutathione_regeneration_cover",
         "NADPH produced per NADPH consumed by glutathione reductases",
         r"P_{NADPH}/C_{NADPH}^{GR}", "ratio"),
        ("ros_glutathione_cycle_closure",
         "glutathione disulfide reduced by glutathione reductases per disulfide formed by peroxide removal",
         r"\sum_{GR}v/\sum_{GPX}v", "ratio"),
    ]

    def __init__(self, ctx):
        super().__init__(ctx)
        r = self.r = _resolved(ctx)
        self.h2o2 = r.carrier_flow(("h2o2",))
        self.o2s = r.carrier_flow(("o2s",))
        self.nadph = r.carrier_flow(("nadph",))
        self.catalase = r.signature(("h2o2",), ("o2",), only=True)
        self.gpx = r.signature(("h2o2", "gthrd"), ("gthox",))
        self.prdx = r.signature(("h2o2",), partner_prefix="prdx")
        self.sod = r.signature(("o2s",), ("h2o2",))
        gr_nadph = r.signature(("gthox", "nadph"), ("gthrd",))
        gr_nadh = r.signature(("gthox", "nadh"), ("gthrd",))
        self.gr = np.union1d(gr_nadph, gr_nadh).astype(int)
        self.gr_nadph = gr_nadph
        has = lambda x: x is not None and bool(np.any(x))
        need = {
            "ros_h2o2_production": has(self.h2o2), "ros_h2o2_removal": has(self.h2o2),
            "ros_h2o2_removal_per_production": has(self.h2o2),
            "ros_h2o2_removed_by_catalase": has(self.h2o2) and self.catalase.size > 0,
            "ros_h2o2_removed_by_glutathione_peroxidase": has(self.h2o2) and self.gpx.size > 0,
            "ros_h2o2_removed_by_peroxiredoxin": has(self.h2o2) and self.prdx.size > 0,
            "ros_superoxide_production": has(self.o2s),
            "ros_superoxide_dismutated": has(self.o2s) and self.sod.size > 0,
            "ros_glutathione_regeneration_cover": has(self.nadph) and gr_nadph.size > 0,
            "ros_glutathione_cycle_closure": self.gr.size > 0 and self.gpx.size > 0,
        }
        why = {"ros_h2o2_removed_by_catalase": "no catalase signature (2 H2O2 -> O2)",
               "ros_h2o2_removed_by_glutathione_peroxidase": "no glutathione peroxidase signature",
               "ros_h2o2_removed_by_peroxiredoxin": "no peroxiredoxin signature",
               "ros_superoxide_dismutated": "no superoxide dismutase signature",
               "ros_glutathione_regeneration_cover": "no NADPH-dependent glutathione reductase",
               "ros_glutathione_cycle_closure": "no glutathione reductase or peroxidase"}
        self.availability = {n: (ok, "" if ok else why.get(n, "no peroxide or superoxide"))
                             for n, ok in need.items()}

    def specs(self):
        return [_fspec(n, "oxidative_defense", "internal", d, f, u) for n, d, f, u in self.SPECS]

    def compute(self, V_int, V_ext):
        n = V_int.shape[0]
        nan = np.full(n, np.nan)
        res = {name: nan for name, *_ in self.SPECS}
        if self.h2o2 is not None and np.any(self.h2o2):
            made, used = _sides(V_int, self.h2o2)
            P, D = made.sum(axis=1), used.sum(axis=1)
            res["ros_h2o2_production"], res["ros_h2o2_removal"] = P, D
            res["ros_h2o2_removal_per_production"] = safe_divide(D, P)
            for name, idx in (("ros_h2o2_removed_by_catalase", self.catalase),
                              ("ros_h2o2_removed_by_glutathione_peroxidase", self.gpx),
                              ("ros_h2o2_removed_by_peroxiredoxin", self.prdx)):
                if idx.size:
                    res[name] = safe_divide(used[:, idx].sum(axis=1), D)
        if self.o2s is not None and np.any(self.o2s):
            made, used = _sides(V_int, self.o2s)
            res["ros_superoxide_production"] = made.sum(axis=1)
            if self.sod.size:
                res["ros_superoxide_dismutated"] = safe_divide(used[:, self.sod].sum(axis=1),
                                                               used.sum(axis=1))
        if self.nadph is not None and self.gr_nadph.size:
            made, used = _sides(V_int, self.nadph)
            res["ros_glutathione_regeneration_cover"] = safe_divide(
                made.sum(axis=1), used[:, self.gr_nadph].sum(axis=1))
        if self.gr.size and self.gpx.size:
            res["ros_glutathione_cycle_closure"] = safe_divide(
                np.maximum(V_int[:, self.gr], 0).sum(axis=1), np.maximum(V_int[:, self.gpx], 0).sum(axis=1))
        return pd.DataFrame(res)

    def support(self, feature: str) -> set[str]:
        r = self.r
        nz = lambda c: set() if c is None else r.reactions_of(np.flatnonzero(c))
        if feature.startswith("ros_h2o2"):
            return nz(self.h2o2)
        if feature.startswith("ros_superoxide"):
            return nz(self.o2s)
        if feature == "ros_glutathione_regeneration_cover":
            return nz(self.nadph)
        return r.reactions_of(np.union1d(self.gr, self.gpx))


# ================================================================ allocation
class AllocationFamily(FeatureFamily):
    """How ATP, NADPH and internal flux are divided among cellular functions."""

    name = "allocation"
    ATP_TO = ("transport", "biosynthesis", "energy", "redox")
    NADPH_TO = ("redox", "biosynthesis")
    FLUX_TO = ("energy", "redox", "biosynthesis", "transport")

    def __init__(self, ctx):
        super().__init__(ctx)
        r = self.r = _resolved(ctx)
        internal = ctx.internal_reactions()
        transport = (internal["reaction_class"].astype(str).str.lower().str.contains("transport").to_numpy()
                     | (r.pathway == "transport"))
        self.function = {f: np.isin(r.pathway, groups) for f, groups in FUNCTIONS.items()}
        self.function["transport"] = transport
        for f in ("energy", "redox", "biosynthesis"):
            self.function[f] = self.function[f] & ~transport
        self.atp = r.carrier_flow(("atp",))
        self.nadph = r.carrier_flow(("nadph",))
        has = lambda c: c is not None and bool(np.any(c))
        self.availability = {}
        for f in self.ATP_TO:
            ok = has(self.atp) and bool(self.function[f].any())
            self.availability[f"alloc_atp_to_{f}"] = (ok, "" if ok else f"no ATP or no {f} reactions")
        for f in self.NADPH_TO:
            ok = has(self.nadph) and bool(self.function[f].any())
            self.availability[f"alloc_nadph_to_{f}"] = (ok, "" if ok else f"no NADPH or no {f} reactions")
        for f in self.FLUX_TO:
            ok = bool(self.function[f].any())
            self.availability[f"alloc_flux_to_{f}"] = (ok, "" if ok else f"no {f} reactions")

    def specs(self):
        out = [_fspec(f"alloc_atp_to_{f}", "allocation", "internal",
                     f"share of gross ATP consumption by {f} reactions",
                     rf"C_{{ATP}}^{{{f}}}/C_{{ATP}}", "ratio") for f in self.ATP_TO]
        out += [_fspec(f"alloc_nadph_to_{f}", "allocation", "internal",
                      f"share of gross NADPH consumption by {f} reactions",
                      rf"C_{{NADPH}}^{{{f}}}/C_{{NADPH}}", "ratio") for f in self.NADPH_TO]
        out += [_fspec(f"alloc_flux_to_{f}", "allocation", "internal",
                      f"share of internal absolute flux carried by {f} reactions",
                      rf"\sum_{{i\in {f}}}|v_i|/\|v\|_1", "ratio") for f in self.FLUX_TO]
        return out

    def compute(self, V_int, V_ext):
        n = V_int.shape[0]
        res = {}
        for carrier, coefficient, targets in (("atp", self.atp, self.ATP_TO), ("nadph", self.nadph, self.NADPH_TO)):
            if coefficient is None or not np.any(coefficient):
                for f in targets:
                    res[f"alloc_{carrier}_to_{f}"] = np.full(n, np.nan)
                continue
            _, used = _sides(V_int, coefficient)
            total = used.sum(axis=1)
            for f in targets:
                res[f"alloc_{carrier}_to_{f}"] = safe_divide(used[:, self.function[f]].sum(axis=1), total)
        A = np.abs(V_int)
        total = A.sum(axis=1)
        for f in self.FLUX_TO:
            res[f"alloc_flux_to_{f}"] = safe_divide(A[:, self.function[f]].sum(axis=1), total)
        return pd.DataFrame(res)

    def support(self, feature: str) -> set[str]:
        if feature.startswith("alloc_flux"):
            return set(self.r.int_rxns)
        coefficient = self.atp if feature.startswith("alloc_atp") else self.nadph
        return set() if coefficient is None else self.r.reactions_of(np.flatnonzero(coefficient))


# ================================================================ conversion
class ConversionFamily(FeatureFamily):
    """What the interior makes per atom the boundary supplies."""

    name = "conversion"
    SPECS = [
        ("conv_atp_per_carbon", "gross ATP production per carbon atom taken up", r"P_{ATP}/U_C", "ratio"),
        ("conv_nadph_per_carbon", "gross NADPH production per carbon atom taken up", r"P_{NADPH}/U_C", "ratio"),
        ("conv_lactate_per_glucose_carbon", "lactate carbon secreted per glucose carbon taken up",
         r"X_C^{lac}/U_C^{glc}", "ratio"),
        ("conv_co2_per_carbon", "CO2 carbon released per carbon atom taken up", r"X_C^{CO_2}/U_C", "ratio"),
    ]

    def __init__(self, ctx):
        super().__init__(ctx)
        r = self.r = _resolved(ctx)
        self.carbon = np.array([element_count(f, "C") for f in r.ext_formula])
        self.atp, self.nadph = r.carrier_flow(("atp",)), r.carrier_flow(("nadph",))
        self.glucose = r.ext_where(("glc__D", "glc_D"))
        self.lactate = r.ext_where(("lac__L", "lac_L"))
        self.co2 = r.ext_where(("co2",))
        has_c = bool((self.carbon > 0).any())
        need = {"conv_atp_per_carbon": (has_c and self.atp is not None, "no ATP or no carbon exchange"),
                "conv_nadph_per_carbon": (has_c and self.nadph is not None, "no NADPH or no carbon exchange"),
                "conv_lactate_per_glucose_carbon": (self.glucose.size > 0 and self.lactate.size > 0,
                                                    "no glucose or lactate exchange"),
                "conv_co2_per_carbon": (has_c and self.co2.size > 0, "no CO2 exchange")}
        self.availability = {k: (ok, "" if ok else why) for k, (ok, why) in need.items()}

    def specs(self):
        return [_fspec(n, "conversion", "combined", d, f, u) for n, d, f, u in self.SPECS]

    def compute(self, V_int, V_ext):
        n = V_int.shape[0]
        x = self.r.export(V_ext)
        up, out = np.maximum(-x, 0.0) * self.carbon, np.maximum(x, 0.0) * self.carbon
        UC = up.sum(axis=1)
        nan = np.full(n, np.nan)
        res = {}
        res["conv_atp_per_carbon"] = (safe_divide(_sides(V_int, self.atp)[0].sum(axis=1), UC)
                                      if self.atp is not None else nan)
        res["conv_nadph_per_carbon"] = (safe_divide(_sides(V_int, self.nadph)[0].sum(axis=1), UC)
                                        if self.nadph is not None else nan)
        res["conv_lactate_per_glucose_carbon"] = (
            safe_divide(out[:, self.lactate].sum(axis=1), up[:, self.glucose].sum(axis=1))
            if self.glucose.size and self.lactate.size else nan)
        res["conv_co2_per_carbon"] = (safe_divide(out[:, self.co2].sum(axis=1), UC)
                                      if self.co2.size else nan)
        return pd.DataFrame(res)

    def support(self, feature: str) -> set[str]:
        r = self.r
        carbon = set(r.ext_rxns[self.carbon > 0])
        if feature == "conv_atp_per_carbon" and self.atp is not None:
            return carbon | r.reactions_of(np.flatnonzero(self.atp))
        if feature == "conv_nadph_per_carbon" and self.nadph is not None:
            return carbon | r.reactions_of(np.flatnonzero(self.nadph))
        if feature == "conv_lactate_per_glucose_carbon":
            return set(r.ext_rxns[np.r_[self.glucose, self.lactate].astype(int)])
        return carbon


FUNCTIONAL_FAMILIES = (ElementalFamily, AminoAcidFamily, CofactorSourceFamily,
                       OxidativeDefenseFamily, AllocationFamily, ConversionFamily)
FUNCTIONAL_PREFIXES = ("elem_", "div_", "aa_", "src_", "ros_", "alloc_", "conv_")

_FUNCTIONAL_CACHE: dict = {}


def functional_families(gem, reactions: pd.DataFrame) -> list:
    """The functional families resolved on one reconstruction (cached by identity)."""
    key = (id(gem), len(reactions))
    if key not in _FUNCTIONAL_CACHE:
        ctx = FeatureContext(gem=gem, reactions=reactions)
        _FUNCTIONAL_CACHE[key] = [family(ctx) for family in FUNCTIONAL_FAMILIES]
    return _FUNCTIONAL_CACHE[key]


def functional_support(feature: str, gem, reactions: pd.DataFrame) -> set[str] | None:
    """The reactions a functional descriptor reads, or None if it is not one."""
    if not feature.startswith(FUNCTIONAL_PREFIXES):
        return None
    for family in functional_families(gem, reactions):
        if feature in family.availability:
            return family.support(feature)
    return None


def functional_availability(gem, reactions: pd.DataFrame) -> dict[str, tuple[bool, str]]:
    out = {}
    for family in functional_families(gem, reactions):
        out.update(family.availability)
    return out


import numpy as np
import pandas as pd


class FeatureEngine:
    """Compute the full candidate feature space for a set of flux vectors."""

    def __init__(self, gem: GEM, reactions: pd.DataFrame):
        self.ctx = FeatureContext(gem=gem, reactions=reactions)
        self.structural = StructuralFamily(self.ctx)
        self.pathway = PathwayFamily(self.ctx)
        self.biochemical = BiochemicalFamily(self.ctx)
        self.network = NetworkFamily(self.ctx)
        self.external = ExternalFamily(self.ctx)
        self.cofactor = CofactorFamily(self.ctx)
        self.capacity = CapacityFamily(self.ctx)
        self.coupling = CouplingFamily(self.ctx)
        self.budget = BudgetFamily(self.ctx)
        self.compartment = CompartmentFamily(self.ctx)
        self.functional = [family(self.ctx) for family in FUNCTIONAL_FAMILIES]
        self.families = [
            self.structural, self.pathway, self.biochemical,
            self.network, self.external,
            self.cofactor, self.budget, self.compartment, *self.functional,
            self.capacity, self.coupling,
        ]
        #: Declared descriptors this reconstruction cannot compute, by structure:
        #: they stay in the registry and are absent from the table of states.
        self.unavailable = dict(sorted(
            (name, why)
            for family in (self.compartment, *self.functional)
            for name, (ok, why) in family.availability.items() if not ok))
        self.registry = FeatureRegistry()
        for family in self.families:
            self.registry.extend(family.specs())

    # ------------------------------------------------------------------
    def transform(self, V: np.ndarray, *, block_size: int = 4000) -> pd.DataFrame:
        """Map ``(n_samples, n_reactions)`` fluxes to ``(n_samples, n_features)``.

        Computation is chunked over samples so that memory stays bounded for
        the full 300,000-sample dataset.
        """
        chunks = []
        for start in range(0, V.shape[0], block_size):
            block = V[start:start + block_size]
            V_int, V_ext = self.ctx.split(block)
            parts = [family.compute(V_int, V_ext) for family in self.families]
            frame = pd.concat(parts, axis=1)
            frame = pd.concat([frame, self.biochemical.add_efficiencies(frame)], axis=1)
            chunks.append(frame)
        out = pd.concat(chunks, ignore_index=True)
        out = out.loc[:, ~out.columns.duplicated()]
        return out.drop(columns=[c for c in self.unavailable if c in out.columns])

    def transform_vector(self, v: np.ndarray) -> pd.Series:
        """Phenotype of a single flux vector (e.g. an FBA reference point)."""
        return self.transform(np.asarray(v, dtype=float).reshape(1, -1)).iloc[0]

    # ------------------------------------------------------------------
    def registry_frame(self) -> pd.DataFrame:
        """Every declared descriptor, with whether this reconstruction can compute it."""
        frame = self.registry.to_frame()
        frame["computable"] = ~frame["feature"].isin(list(self.unavailable))
        frame["unavailable_because"] = frame["feature"].map(self.unavailable).fillna("")
        return frame

    def coverage_report(self) -> pd.DataFrame:
        """Model entities the engine expected but could not find."""
        rows = [
            {"kind": "metabolite", "identifier": m, "consequence": "turnover feature omitted"}
            for m in self.biochemical.missing_metabolites
        ] + [
            {"kind": "reaction", "identifier": r, "consequence": "key-reaction feature omitted"}
            for r in self.biochemical.missing_reactions
        ]
        return pd.DataFrame(rows, columns=["kind", "identifier", "consequence"])


# --------------------------------------------------------------------------
# Interface
# --------------------------------------------------------------------------

def engine(gem: "GEM") -> "FeatureEngine":
    """The descriptor engine for one reconstruction."""
    return FeatureEngine(gem, annotate_reactions(gem))


def declared(gem: "GEM") -> pd.DataFrame:
    """Every declared descriptor, with whether this reconstruction defines it and why not."""
    built = engine(gem)
    table = built.registry_frame()
    table["defined"] = ~table["feature"].isin(built.unavailable)
    table["why_not"] = table["feature"].map(built.unavailable).fillna("")
    return table


def phi(gem: "GEM", V, *, built: "FeatureEngine | None" = None) -> pd.DataFrame:
    """Phi_M(v) of every state in V (states x reactions): every declared column,
    NaN where the reconstruction cannot define the descriptor."""
    built = built or engine(gem)
    V = np.atleast_2d(np.asarray(V, dtype=float))
    if V.shape[1] != len(gem.rxns) and V.shape[0] == len(gem.rxns):
        V = V.T
    values = built.transform(V)
    columns = built.registry_frame()["feature"].tolist()
    return values.reindex(columns=columns)
