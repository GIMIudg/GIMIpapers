"""Funciones de las Fases 6 y 7: capas del InChI, tipos de inconsistencia, variables para
clustering y recomendación del InChI. Pruebas en pruebas_inchi_calidad.py.

Idea central: el InChI estándar separa la estructura en capas. Comparar las fuentes capa por
capa dice QUÉ difiere (fórmula, conectividad, hidrógenos, estereoquímica o solo protonación),
no solo SI difiere.
"""
import re
from collections import Counter
from functools import lru_cache

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")

FUENTES = ["INCHI_MODELO", "INCHI_SMILES", "INCHI_PUBCHEM", "INCHI_CHEBI", "INCHI_KEGG", "INCHI_HMDB"]
DEL_MODELO = ["INCHI_MODELO", "INCHI_SMILES"]          # salen del mismo GEM: cuentan como un votante
EXTERNAS = ["INCHI_PUBCHEM", "INCHI_CHEBI", "INCHI_KEGG", "INCHI_HMDB"]   # HMDB solo si hay .mol locales

# Capas en orden de "más grave" a "menos grave". Cada nivel agrega capas al anterior.
NIVELES = [
    ("formula", ["formula"]),
    ("conectividad", ["c"]),
    ("hidrogenos", ["h"]),
    ("estereoquimica", ["b", "t", "m", "s", "i"]),
    ("carga", ["q", "p"]),
]


# ---------------------------------------------------------------- lectura de capas
@lru_cache(maxsize=None)
def a_estandar(inchi):
    """Devuelve el InChI estándar (InChI=1S/...). Los no estándar (InChI=1/) se convierten con RDKit."""
    if not isinstance(inchi, str) or not inchi.startswith("InChI="):
        return None
    inchi = inchi.strip()
    if inchi.startswith("InChI=1S/"):
        return inchi
    mol = Chem.MolFromInchi(inchi, sanitize=True)
    if mol is None:
        return None
    out = Chem.MolToInchi(mol) or None
    # /s2 = estereo relativa y /s3 = racémica: no dicen la configuración absoluta. RDKit la inventa al
    # convertir (/m0/s1), así que se quitan /t /m /s y el centro queda "no definido". Los E/Z (/b) se conservan.
    if out and re.search(r"/s[23](/|$)", inchi):
        partes = out.split("/")
        out = "/".join(partes[:2] + [x for x in partes[2:] if x[:1] not in ("t", "m", "s")])
    return out


def capas(inchi):
    """Diccionario capa -> texto. 'formula' es la capa principal sin prefijo."""
    if not isinstance(inchi, str) or not inchi.startswith("InChI="):
        return None
    partes = inchi.split("/")
    # El protón (InChI=1S/p+1) no tiene capa de fórmula: la segunda parte ya es una capa con letra.
    resto = partes[2:] if len(partes) > 1 and partes[1][:1].isupper() else partes[1:]
    d = {"formula": partes[1] if len(partes) > 1 and partes[1][:1].isupper() else ""}
    for p in resto:
        if p:
            d.setdefault(p[0], p[1:])      # primera aparición (las de /f fixed-H no se usan)
    return d


def clave(inchi, hasta):
    """Texto con las capas desde 'formula' hasta el nivel 'hasta' (inclusive)."""
    c = capas(inchi)
    if c is None:
        return None
    partes = []
    for nombre, letras in NIVELES:
        partes += [f"{l}={c.get(l, '')}" for l in letras]
        if nombre == hasta:
            break
    return "|".join(partes)


def n_capas(inchi):
    c = capas(inchi)
    return None if c is None else len(c)


# ---------------------------------------------------------------- carga y fórmula
def suma_capa_carga(valor):
    """Suma una capa /q o /p. Acepta varios componentes separados por ';' y la notación 'N*carga'
    (N componentes iguales con esa carga; p. ej. 'q2*-1;+2' = -2 + 2 = 0)."""
    total = 0
    for parte in (valor or "").split(";"):
        m = re.fullmatch(r"\s*(?:(\d+)\*)?([+-]?\d+)\s*", parte)
        if m:
            total += int(m.group(1) or 1) * int(m.group(2))
    return total


def carga_neta(inchi):
    """Carga neta = suma de /q (por componente) + protones de /p."""
    c = capas(inchi)
    if c is None:
        return None
    return suma_capa_carga(c.get("q")) + suma_capa_carga(c.get("p"))


_ELEM = re.compile(r"([A-Z][a-z]?)(\d*)")


def _elementos_validos():
    pt = Chem.GetPeriodicTable()
    out = set()
    for z in range(1, 119):
        out.add(pt.GetElementSymbol(z))
    return out


ELEMENTOS = _elementos_validos()


def parsear_formula(f):
    """'C2H3O2' -> {'C':2,'H':3,'O':2}. Fórmulas con varios componentes ('2C2H4O2.Na') se suman.
    Devuelve None si trae símbolos que no son elementos (R, X, FULLR...)."""
    if not isinstance(f, str) or not f.strip():
        return None
    total = Counter()
    for comp in f.strip().split("."):
        m = re.match(r"^(\d+)(.*)$", comp)
        mult, comp = (int(m.group(1)), m.group(2)) if m else (1, comp)
        if re.sub(r"[A-Z][a-z]?\d*", "", comp):
            return None
        for sim, n in _ELEM.findall(comp):
            if sim not in ELEMENTOS:
                return None
            total[sim] += mult * (int(n) if n else 1)
    return dict(total)


def hill(d):
    orden = (["C", "H"] if "C" in d else []) + sorted(k for k in d if not ("C" in d and k in ("C", "H")))
    return "".join(f"{k}{d[k] if d[k] != 1 else ''}" for k in orden if d.get(k, 0) > 0)


def _formula_de_capas(c):
    """Fórmula del InChI como diccionario. Vacía = cero átomos (el protón), no 'desconocida'."""
    return {} if c["formula"] == "" else parsear_formula(c["formula"])


def formula_corregida(inchi):
    """Fórmula del InChI ajustando los hidrógenos por la capa /p (la que reporta el GEM)."""
    c = capas(inchi)
    if c is None:
        return None
    d = _formula_de_capas(c)
    if d is None:
        return None
    d["H"] = d.get("H", 0) + suma_capa_carga(c.get("p"))
    return hill(d)


def coincide_formula_gem(inchi, formula_gem):
    """True/False si la fórmula del InChI (con su protonación) es la del GEM; None si no se puede comparar."""
    g = parsear_formula(formula_gem)
    f = formula_corregida(inchi)
    if g is None or f is None:
        return None
    return parsear_formula(f) == g


def coincide_pesados_gem(inchi, formula_gem):
    """True/False si coinciden los átomos distintos de H (no depende de la protonación); None si no se puede."""
    g = parsear_formula(formula_gem)
    c = capas(inchi)
    f = _formula_de_capas(c) if c else None
    if g is None or f is None:
        return None
    sin_h = lambda d: {k: v for k, v in d.items() if k != "H"}
    return sin_h(f) == sin_h(g)


def formula_generica(formula_gem):
    """True si la fórmula del GEM no es química (trae R, X, FULLR...): clases de lípidos, polímeros genéricos."""
    return isinstance(formula_gem, str) and bool(formula_gem.strip()) and parsear_formula(formula_gem) is None


# ---------------------------------------------------------------- estructura neutralizada y estereo
from rdkit.Chem import rdCIPLabeler
from rdkit.Chem.MolStandardize import rdMolStandardize

_UNCHARGER = rdMolStandardize.Uncharger()


def sin_carga(inchi):
    """El InChI sin la capa /p. La /q se conserva: es carga fija (p. ej. el N cuaternario de la colina)."""
    partes = inchi.split("/")
    return "/".join(partes[:2] + [p for p in partes[2:] if p[:1] != "p"])


@lru_cache(maxsize=None)
def _mol_neutra(inchi):
    """Molécula neutralizada con el Uncharger de RDKit (agrega o quita H donde se puede).
    Hace falta porque el InChI estándar no siempre marca la protonación con /p (p. ej. arsenito AsO3 q-3
    frente a ácido arsenioso AsH3O3)."""
    m = Chem.MolFromInchi(inchi)
    if m is None:
        return None
    try:
        return _UNCHARGER.uncharge(m)
    except Exception:
        return m


@lru_cache(maxsize=None)
def neutro(inchi):
    """InChI estándar de la estructura neutralizada."""
    m = _mol_neutra(inchi)
    return (Chem.MolToInchi(m) if m is not None else None) or sin_carga(inchi)


@lru_cache(maxsize=None)
def _estereo(inchi):
    """Configuraciones definidas: {('a', átomo): 'R'/'S', ('b', i, j): 'E'/'Z'} con el etiquetador CIP de RDKit.
    Los índices siguen la numeración canónica del InChI, así que son comparables entre fuentes que
    comparten la misma capa principal."""
    m = _mol_neutra(inchi)
    if m is None:
        return ()
    m = Chem.Mol(m)
    try:
        rdCIPLabeler.AssignCIPLabels(m)
    except Exception:
        return ()
    d = [(("a", a.GetIdx()), a.GetProp("_CIPCode")) for a in m.GetAtoms() if a.HasProp("_CIPCode")]
    d += [(("b", b.GetBeginAtomIdx(), b.GetEndAtomIdx()), b.GetProp("_CIPCode"))
          for b in m.GetBonds() if b.HasProp("_CIPCode")]
    return tuple(d)


def estereo(inchi):
    return dict(_estereo(inchi))


def principal(inchi):
    """Fórmula + conectividad + hidrógenos de la estructura neutralizada."""
    return clave(neutro(inchi), "hidrogenos")


def contradice(a, b):
    """Misma constitución y al menos un centro (o doble enlace) con configuración opuesta."""
    if principal(a) != principal(b):
        return False
    ea, eb = estereo(a), estereo(b)
    return any(ea[k] != eb[k] for k in set(ea) & set(eb))


def compatibles(a, b):
    """Misma constitución y ninguna configuración opuesta (una puede definir menos centros que la otra)."""
    return principal(a) == principal(b) and not contradice(a, b)


def especificidad(inchi):
    return len(_estereo(inchi))


def _estado_carga(inchi):
    c = capas(inchi)
    return (c.get("formula"), c.get("q", ""), c.get("p", "")) if c else None


# ---------------------------------------------------------------- inconsistencias
TIPOS = ["formula", "conectividad", "anillo_cadena", "hidrogenos", "estereo_contradictoria", "estereo_incompleta",
         "carga"]


@lru_cache(maxsize=None)
def n_anillos(inchi):
    m = _mol_neutra(inchi)
    return None if m is None else m.GetRingInfo().NumRings()


def _presentes(inchis):
    out = {}
    for f in FUENTES:
        v = a_estandar(inchis.get(f))
        if v:
            out[f] = v
    return out


def tipos_inconsistencia(inchis, normalizar=True):
    """Tipos presentes, en orden de gravedad. [] = todas coinciden. None = menos de 2 fuentes.
    Fórmula, conectividad e hidrógenos se comparan sobre la estructura neutralizada (capas anidadas: un tipo
    está si su capa separa fuentes que la anterior juntaba). La estereoquímica se compara con R/S y E/Z
    entre fuentes de igual constitución. 'carga' = misma constitución, distinto estado de protonación/carga."""
    p = _presentes(inchis) if normalizar else {k: v for k, v in inchis.items() if v}
    if len(p) < 2:
        return None
    tipos, previo = [], 1
    for nivel in ["formula", "conectividad", "hidrogenos"]:
        n = len({clave(neutro(v), nivel) for v in p.values()})
        if n > previo:
            tipos.append(nivel)
        previo = n
    if "conectividad" in tipos:
        # Misma fórmula, distinta conectividad: ¿es forma abierta vs cíclica (azúcares) o un isómero distinto?
        mismos = [(a, b) for i, a in enumerate(p.values()) for b in list(p.values())[i + 1:]
                  if clave(neutro(a), "formula") == clave(neutro(b), "formula")
                  and clave(neutro(a), "conectividad") != clave(neutro(b), "conectividad")]
        anillo = [n_anillos(a) != n_anillos(b) for a, b in mismos]
        if any(anillo):
            tipos.insert(tipos.index("conectividad") + 1, "anillo_cadena")
            if all(anillo):
                tipos.remove("conectividad")
    pares = [(a, b) for i, a in enumerate(p.values()) for b in list(p.values())[i + 1:] if principal(a) == principal(b)]
    if any(contradice(a, b) for a, b in pares):
        tipos.append("estereo_contradictoria")
    if any(not contradice(a, b) and set(estereo(a)) != set(estereo(b)) for a, b in pares):
        tipos.append("estereo_incompleta")
    if any(_estado_carga(a) != _estado_carga(b) for a, b in pares):
        tipos.append("carga")
    return tipos


# ---------------------------------------------------------------- consenso con compatibilidad
def _ganador(fuentes):
    """La estructura con más apoyo. Una fuente apoya a otra si son compatibles (misma constitución, sin
    configuraciones opuestas). Solo compiten las más específicas: una fuente que no define un centro no
    puede ganarle a una compatible que sí lo define. Devuelve (nombre, [nombres que la apoyan])."""
    nombres = list(fuentes)
    comp = {(a, b): compatibles(fuentes[a], fuentes[b]) for a in nombres for b in nombres}
    esp = {a: especificidad(fuentes[a]) for a in nombres}
    cand = [a for a in nombres if not any(comp[a, b] and esp[b] > esp[a] for b in nombres)]
    apoyo = {a: [b for b in nombres if comp[a, b]] for a in cand}
    mejor = max(cand, key=lambda a: (len(apoyo[a]), esp[a], -nombres.index(a)))
    return mejor, apoyo[mejor]


# ---------------------------------------------------------------- similitud
_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


@lru_cache(maxsize=None)
def _huella(inchi):
    m = _mol_neutra(inchi)
    return None if m is None else _GEN.GetFingerprint(m)


def tanimoto_minimo(inchis):
    fps = [fp for fp in (_huella(v) for v in inchis.values()) if fp is not None]
    if len(fps) < 2:
        return np.nan
    return float(min(min(DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1:])) for i in range(len(fps) - 1)))


# ---------------------------------------------------------------- variables por metabolito
def _frac_moda(valores):
    if not valores:
        return np.nan
    return Counter(valores).most_common(1)[0][1] / len(valores)


def variables(inchis, formula_gem, carga_gem):
    """Variables de consistencia de un metabolito (las que lista el tutor)."""
    p = _presentes(inchis)
    fila = {
        "n_fuentes": len(p),
        "n_externas": sum(1 for f in EXTERNAS if f in p),
        "n_inchis_distintos": len(set(p.values())),
        "n_estructuras_distintas": len({neutro(v) for v in p.values()}),
        "consenso_inchi_completo": _frac_moda(list(p.values())),
        "consenso_capa_principal": _frac_moda([principal(v) for v in p.values()]),
        "consenso_estructura": (len(_ganador(p)[1]) / len(p)) if p else np.nan,
        "n_capas_media": float(np.mean([n_capas(v) for v in p.values()])) if p else np.nan,
        "tanimoto_minimo": tanimoto_minimo(p),
    }
    fg = [coincide_formula_gem(v, formula_gem) for v in p.values()]
    fg = [x for x in fg if x is not None]
    fila["frac_formula_gem"] = float(np.mean(fg)) if fg else np.nan
    fp = [coincide_pesados_gem(v, formula_gem) for v in p.values()]
    fp = [x for x in fp if x is not None]
    fila["frac_pesados_gem"] = float(np.mean(fp)) if fp else np.nan
    fila["formula_gem_generica"] = formula_generica(formula_gem)
    try:
        cg = int(carga_gem)
        fila["frac_carga_gem"] = float(np.mean([carga_neta(v) == cg for v in p.values()])) if p else np.nan
    except (TypeError, ValueError):
        fila["frac_carga_gem"] = np.nan
    tipos = tipos_inconsistencia(p, normalizar=False)
    for t in TIPOS:
        fila[f"inc_{t}"] = np.nan if tipos is None else int(t in tipos)
    fila["tipos_inconsistencia"] = "sin_comparacion" if tipos is None else (", ".join(tipos) or "ninguna")
    fila["tipo_principal"] = "sin_comparacion" if tipos is None else (tipos[0] if tipos else "consistente")
    interno = tipos_inconsistencia({k: p[k] for k in DEL_MODELO if k in p}, normalizar=False)
    fila["gem_interno"] = "sin_comparacion" if interno is None else (", ".join(interno) or "ninguna")
    return fila


# ---------------------------------------------------------------- recomendación (Fase 7)
def recomendar(inchis, carga_gem):
    """Propuesta de InChI recomendado. Reglas en fase6_7/PLAN.md (pendientes de aprobar por el tutor)."""
    p = _presentes(inchis)
    votos = {}
    gem = {k: p[k] for k in DEL_MODELO if k in p}
    if gem:   # el GEM vota una vez: su InChI, o el más específico de sus dos formas si son compatibles
        g = gem.get("INCHI_MODELO") or gem["INCHI_SMILES"]
        if len(gem) == 2 and compatibles(*gem.values()):
            g = max(gem.values(), key=especificidad)
        votos["GEM"] = g
    for f in EXTERNAS:
        if f in p:
            votos[f] = p[f]
    sal = {"inchi": None, "confianza": "revision_manual", "prioridad_revision": "baja",
           "motivo": "", "gem_coincide": None, "n_votantes": len(votos), "apoyo": 0, "fuentes_a_favor": ""}
    if not votos:
        sal["confianza"], sal["motivo"] = "sin_datos", "sin InChI en ninguna fuente"
        return sal
    mejor, a_favor = _ganador(votos)
    ext_a_favor = [k for k in a_favor if k != "GEM"]
    sal["apoyo"], sal["fuentes_a_favor"] = len(a_favor), ", ".join(k.replace("INCHI_", "") for k in a_favor)
    if len(votos) < 2:
        sal["motivo"] = "solo fuentes del modelo" if "GEM" in votos else "una sola base externa"
        return sal
    if len(a_favor) * 2 <= len(votos) or not ext_a_favor:
        sal["motivo"] = "fuentes en conflicto, sin mayoría"
        sal["prioridad_revision"] = "media"
        return sal
    if "GEM" in votos:
        sal["gem_coincide"] = "GEM" in a_favor
    elegido = votos[mejor]
    try:
        cg = int(carga_gem)
    except (TypeError, ValueError):
        cg = None
    misma = [v for v in p.values() if neutro(v) == neutro(elegido)]
    con_carga = [v for v in misma if cg is not None and carga_neta(v) == cg]
    sal["inchi"] = (con_carga or [elegido])[0]
    n_ext = sum(1 for k in votos if k != "GEM")
    sal["confianza"] = "alta" if (len(a_favor) == len(votos) and n_ext >= 2) else "media"
    if "GEM" in votos and "GEM" not in a_favor and len(ext_a_favor) >= 2:
        sal["prioridad_revision"] = "alta"
        sal["motivo"] = "la estructura del GEM difiere del consenso de las bases externas"
    elif len(a_favor) == len(votos):
        identicas = len({neutro(votos[k]) for k in a_favor}) == 1
        sal["motivo"] = "todas las fuentes coinciden" if identicas else \
            "coinciden; alguna fuente no define toda la estereoquímica"
    else:
        sal["motivo"] = "mayoría a favor"
    return sal


def base_id(met_id):
    """'atp[c]' -> 'atp' (quita el compartimento entre corchetes al final)."""
    return re.sub(r"\[[^\]]+\]$", "", str(met_id))
