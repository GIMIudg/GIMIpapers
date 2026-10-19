#!/usr/bin/env python3
"""
Convierte los CSV transcriptómicos de islotes pancreáticos
(MCBCI_DiaMED/datos/DE_DONDE_SALIO_ISLOTESPANCREATICOS/*.csv) a archivos .txt
en el formato que espera XomicsToModel (specificData.transcriptomicData):

    genes,expVal
    7105,1.50988667390527
    64102,0.0495192289370227
    ...

- genes : Entrez Gene ID (entero, sin sufijo de transcrito)
- expVal: expresión en escala LINEAL y no negativa (TPM/FPKM). XomicsToModel
          aplica log() internamente y compara contra param.transcriptomicThreshold.

Controles que se aplican a cada archivo:
  * columnas requeridas (genes, expVal); acepta alias comunes (ENTREZ_GENE_ID, Entrez, value...)
  * elimina filas con gen o valor faltante / no numérico
  * Entrez IDs enteros positivos
  * valores negativos -> 0 (se reportan)
  * genes duplicados -> se colapsan con la media (se reportan)
  * detecta datos en escala log (max < 30) y lo advierte

Uso:
    python3 csv2txt_transcriptomics.py                       # rutas por defecto
    python3 csv2txt_transcriptomics.py -i <dir_csv> -o <dir_txt>
"""
import argparse
import glob
import os
import sys

import pandas as pd

DEFAULT_ROOT = os.path.expanduser('~/work/Drive/GIMI/01_posgrado/MCBCI_DiaMED')
DEFAULT_IN = os.path.join(DEFAULT_ROOT, 'datos', 'DE_DONDE_SALIO_ISLOTESPANCREATICOS')
DEFAULT_OUT = os.path.join(DEFAULT_ROOT, 'resultados', 'transcriptomicData')

GENE_ALIASES = ['genes', 'gene', 'entrez', 'entrezid', 'entrez_id', 'entrez_gene_id', 'geneid']
VAL_ALIASES = ['expval', 'value', 'expression', 'tpm', 'fpkm', 'rpkm', 'counts']


def find_col(cols, aliases):
    low = {c.strip().lower(): c for c in cols}
    for a in aliases:
        if a in low:
            return low[a]
    return None


def convert(csv_path, out_dir):
    name = os.path.splitext(os.path.basename(csv_path))[0]
    df = pd.read_csv(csv_path)
    gcol = find_col(df.columns, GENE_ALIASES)
    vcol = find_col(df.columns, VAL_ALIASES)
    if gcol is None or vcol is None:
        if df.shape[1] == 2:
            gcol, vcol = df.columns
        else:
            raise ValueError(f'{name}: columnas no reconocidas {list(df.columns)}')

    out = pd.DataFrame({
        'genes': pd.to_numeric(df[gcol].astype(str).str.replace(r'\.\d+$', '', regex=True), errors='coerce'),
        'expVal': pd.to_numeric(df[vcol], errors='coerce'),
    })
    n0 = len(out)
    out = out.dropna()
    out = out[(out.genes > 0) & (out.genes == out.genes.round())]
    nDropped = n0 - len(out)
    out['genes'] = out['genes'].astype('int64')

    nNeg = int((out.expVal < 0).sum())
    out.loc[out.expVal < 0, 'expVal'] = 0.0

    nDup = int(out.genes.duplicated().sum())
    if nDup:
        out = out.groupby('genes', as_index=False, sort=False)['expVal'].mean()

    logWarning = out.expVal.max() < 30

    txt_path = os.path.join(out_dir, name + '.txt')
    out.to_csv(txt_path, index=False, header=True, float_format='%.15g', lineterminator='\n')
    return dict(sample=name, rows=len(out), dropped=nDropped, negatives=nNeg,
                duplicates=nDup, expressed_ge1=int((out.expVal >= 1).sum()),
                maxExp=out.expVal.max(), possibleLogScale=logWarning)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-i', '--input', default=DEFAULT_IN, help='carpeta con los .csv')
    ap.add_argument('-o', '--output', default=DEFAULT_OUT, help='carpeta destino de los .txt')
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.input, '*.csv')))
    if not files:
        sys.exit(f'No se encontraron .csv en {args.input}')
    os.makedirs(args.output, exist_ok=True)

    report = pd.DataFrame([convert(f, args.output) for f in files])
    report.to_csv(os.path.join(args.output, '_conversionReport.csv'), index=False)

    print(f'{len(report)} archivos convertidos -> {args.output}')
    print(report[['rows', 'dropped', 'negatives', 'duplicates', 'expressed_ge1']].describe().loc[['min', 'mean', 'max']])
    if report.possibleLogScale.any():
        print('ADVERTENCIA: posibles datos en escala log:', list(report.sample[report.possibleLogScale]))


if __name__ == '__main__':
    main()
