"""FluxQuery - model validation (run after step 2 and, for test 3, after step 3).

Tests of iDopaNeuro, iCholiNeuro and iGlutaNeuro against knowledge that was not used to build them:
  1. Specificity: flux variability of the marker reactions of each neurotransmitter system and of
     enzymes absent from neurons. A model passes when only its own transmitter pathway can carry flux.
  2. Exometabolome (neurons with measured rates): each measured exchange rate is removed in turn, its exchange is
     opened (uptake only for medium components) and the flux range allowed by the rest of the model is
     compared with the measured direction and mean.
  3. Energy metabolism (needs samplingSummary.csv from step 3): oxygen/glucose and lactate/glucose
     uptake ratios of the sampled states.
  4. Unlimited sources: uptakes, demands or sinks that can supply metabolites up to the generic flux
     bound, and the ranges of glucose, H2O2 and phosphatidylserine-demand fluxes.

Usage:  python modelValidationFluxQuery.py <project folder>
Output: <project>/results/python/validation/*.csv and validationSummary.txt
Requirements: numpy, scipy (HiGHS), pandas, openpyxl
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
from scipy.optimize import linprog

PROJECT = Path(sys.argv[1] if len(sys.argv) > 1 else '.')
OUT = PROJECT / 'results' / 'python' / 'validation'
MODELS = {'dopaminergic': 'iDopaNeuro', 'cholinergic': 'iCholiNeuro', 'glutamatergic': 'iGlutaNeuro'}
# reaction, process, neurons in which it may carry flux, neurons in which it must carry flux
MARKERS = [('CHAT', 'acetylcholine synthesis', {'cholinergic'}, {'cholinergic'}),
           ('ACHVESSEC', 'vesicular acetylcholine release', {'cholinergic'}, {'cholinergic'}),
           # glutamate co-release by cholinergic and dopaminergic neurons
           ('GLUVESSEC', 'vesicular glutamate release', {'glutamatergic', 'cholinergic', 'dopaminergic'}, {'glutamatergic'}),
           ('TYR3MO2', 'tyrosine hydroxylase (catecholamines)', {'dopaminergic'}, {'dopaminergic'}),
           ('DOPAVESSEC', 'vesicular dopamine release', {'dopaminergic'}, {'dopaminergic'}),
           ('HISDC', 'histidine decarboxylase (histamine)', set(), set()),
           ('GLUDC', 'glutamate decarboxylase (GABA)', set(), set()),
           ('GLNS', 'glutamine synthetase (astrocytes)', set(), set()),
           # pyruvate carboxylase flux was measured in the iDopaNeuro cultures (rxns2constrain of iDopaNeuro)
           ('PCm', 'pyruvate carboxylase (astrocytes)', {'dopaminergic'}, set()),
           ('PGCD', 'phosphoglycerate dehydrogenase (astrocytes)', set(), set()),
           ('NAAGS', 'NAAG synthesis', {'glutamatergic', 'cholinergic', 'dopaminergic'}, set()),
           ('GABAVESSEC', 'vesicular GABA release (VGAT)', set(), set()),
           ('GLYVESSEC', 'vesicular glycine release (VGAT)', set(), set()),
           ('DOPABMO', 'dopamine beta-hydroxylase (noradrenaline)', set(), set()),
           ('CITL', 'ATP-independent citrate lyase (not human)', set(), set()),
           ('HMGCOAtm', 'mitochondrial HMG-CoA transport (ketogenesis bypass)', set(), set())]
TOL = 1e-6


def load(neuron):
    path = PROJECT / 'results' / 'matlab' / 'models' / neuron / f'{MODELS[neuron]}.mat'
    M = sio.loadmat(path, squeeze_me=True, struct_as_record=False)['model']
    m = dict(rxns=list(map(str, M.rxns)), S=sp.csr_matrix(M.S), lb=M.lb.astype(float), ub=M.ub.astype(float))
    m['A'] = m['bA'] = None
    if getattr(M, 'C', np.empty(0)).size:
        sense = np.array(list(''.join(np.atleast_1d(M.dsense))))
        sign = np.where(sense == 'G', -1.0, 1.0)            # C v >= d  ->  -C v <= -d
        m['A'] = sp.diags(sign) @ sp.csr_matrix(M.C)
        m['bA'] = sign * np.atleast_1d(M.d).astype(float)
    return m


def flux_range(m, rxn, lb=None, ub=None):
    lb = m['lb'] if lb is None else lb
    ub = m['ub'] if ub is None else ub
    j = m['rxns'].index(rxn)
    out = []
    for sense in (1, -1):
        c = np.zeros(len(m['rxns'])); c[j] = sense
        r = linprog(c, A_ub=m['A'], b_ub=m['bA'], A_eq=m['S'], b_eq=np.zeros(m['S'].shape[0]),
                    bounds=list(zip(lb, ub)), method='highs')
        out.append(sense * r.fun if r.status == 0 else np.nan)
    return out[0], out[1]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    lines = []
    models = {n: load(n) for n in MODELS if (PROJECT / 'results' / 'matlab' / 'models' / n / f'{MODELS[n]}.mat').exists()}

    # 1. specificity
    rows = []
    for n, m in models.items():
        for rxn, process, expected, required in MARKERS:
            lo, hi = flux_range(m, rxn) if rxn in m['rxns'] else (np.nan, np.nan)
            active = rxn in m['rxns'] and (hi > TOL or lo < -TOL)
            allowed = n in expected
            rows.append(dict(neuron=n, reaction=rxn, process=process, present=rxn in m['rxns'], minFlux=lo, maxFlux=hi,
                             canCarryFlux=active, allowed=allowed, required=n in required,
                             passed=(active or n not in required) and (allowed or not active)))
    spec = pd.DataFrame(rows)
    spec.to_csv(OUT / 'specificity.csv', index=False)
    for n in models:
        s = spec[spec.neuron == n]
        lines.append(f'{n}: specificity {s.passed.sum()}/{len(s)} passed; failed: {", ".join(s.loc[~s.passed, "reaction"]) or "none"}')

    # 2. leave-one-out prediction of the measured exometabolome (neurons with exometabolomic data)
    for n, m in models.items():
        exoFile = PROJECT / 'data' / 'omics' / n / f'exometabolomicData_{n}Neuron.txt'
        if not exoFile.exists():
            continue
        exo = pd.read_csv(exoFile)
        medium = set(pd.read_excel(PROJECT / 'results' / 'bibliomic' / f'bibliomicData_{n}Neuron.xlsx',
                                   sheet_name='mediaData').rxns.astype(str))
        rows = []
        for _, r in exo.iterrows():
            name = r.get('metName', r.rxns)
            measured = 'uptake' if r['mean'] + r.SD < 0 else 'secretion' if r['mean'] - r.SD > 0 else 'not significant'
            if r.rxns not in m['rxns']:
                rows.append(dict(rxns=r.rxns, metabolite=name, mean=r['mean'], SD=r.SD, measured=measured,
                                 predicted='not in model'))
                continue
            j = m['rxns'].index(r.rxns)
            lb, ub = m['lb'].copy(), m['ub'].copy()
            gb = max(np.max(m['ub']), -np.min(m['lb']))
            lb[j] = -gb if r.rxns in medium else 0
            ub[j] = gb
            lo, hi = flux_range(m, r.rxns, lb, ub)
            predicted = 'uptake' if hi < -TOL else 'secretion' if lo > TOL else 'either'
            rows.append(dict(rxns=r.rxns, metabolite=name, mean=r['mean'], SD=r.SD, measured=measured,
                             predictedMin=lo, predictedMax=hi, predicted=predicted,
                             meanInRange=bool(lo - TOL <= r['mean'] <= hi + TOL)))
        loo = pd.DataFrame(rows)
        loo.to_csv(OUT / f'exometabolome_leaveOneOut_{n}.csv', index=False)
        sig = loo[loo.measured != 'not significant']
        decided = sig[sig.predicted.isin(['uptake', 'secretion'])]
        lines.append(f'{n}: {len(sig)} significant exchange rates; direction predicted for {len(decided)}, '
                     f'correct {(decided.predicted == decided.measured).sum()}; measured mean inside the predicted range '
                     f'for {int(loo.meanInRange.fillna(False).sum())}/{len(loo)}')

    # 3. energy metabolism of the sampled states
    for n in models:
        f = PROJECT / 'results' / 'matlab' / 'sampling' / n / 'samplingSummary.csv'
        if not f.exists():
            continue
        s = pd.read_csv(f).set_index('rxns')['mean']
        glc, o2, lac = -s.get('EX_glc_D[e]', np.nan), -s.get('EX_o2[e]', np.nan), s.get('EX_lac_L[e]', np.nan)
        lines.append(f'{n}: mean glucose uptake {glc:.0f}, O2 uptake {o2:.0f}, lactate release {lac:.0f} umol/gDW/h; '
                     f'O2/glucose {o2 / glc:.2f} (complete oxidation 6), lactate/glucose {lac / glc:.2f} (maximum 2)')

    # 4. sources at the generic bound: uptakes or reversed demands/sinks that are only limited by the
    #    generic flux bound (water and protons excluded) supply unlimited carbon, electrons or lipid
    for n, m in models.items():
        gb = max(np.max(m['ub']), -np.min(m['lb']))          # generic bound of this model
        open_ = [r for r, l in zip(m['rxns'], m['lb']) if r.startswith(('EX_', 'DM_', 'sink_')) and l <= -0.99 * gb
                 and r not in ('EX_h2o[e]', 'EX_h[e]')]
        lines.append(f'{n}: {len(open_)} sources limited only by the generic bound: {", ".join(open_) or "none"}')
        for r in ('EX_glc_D[e]', 'EX_h2o2[e]', 'DM_ps_hs[c]'):
            if r in m['rxns']:
                lo, hi = flux_range(m, r)
                lines.append(f'   {r}: [{lo:.1f}, {hi:.1f}]')

    (OUT / 'validationSummary.txt').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
