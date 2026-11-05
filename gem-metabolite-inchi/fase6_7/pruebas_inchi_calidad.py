"""Batería de casos de inchi_calidad.py con InChIs reales conocidos.
Correr: .venv/bin/python fase6_7/pruebas_inchi_calidad.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from inchi_calidad import (capas, carga_neta, formula_corregida, coincide_formula_gem,
                           tipos_inconsistencia, recomendar, base_id, a_estandar, variables)

ACIDO_ACETICO = "InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)"
ACETATO = ACIDO_ACETICO + "/p-1"
L_ALA = "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1"
D_ALA = "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m1/s1"
ALA_SIN_ESTEREO = "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)"
SARCOSINA = "InChI=1S/C3H7NO2/c1-4-2-3(5)6/h4H,2H2,1H3,(H,5,6)"   # misma fórmula que alanina
ETANOL = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
TMA = "InChI=1S/C4H12N/c1-5(2,3)4/h1-4H3/q+1"   # tetrametilamonio, carga fija +1

M, S, P, C, K = "INCHI_MODELO", "INCHI_SMILES", "INCHI_PUBCHEM", "INCHI_CHEBI", "INCHI_KEGG"

casos = []
def caso(nombre, obtenido, esperado):
    casos.append((nombre, obtenido == esperado, obtenido, esperado))

# Capas
c = capas(ACETATO)
caso("fórmula del acetato", c["formula"], "C2H4O2")
caso("protonación del acetato", c["p"], "-1")
caso("estereo de L-alanina", (c2 := capas(L_ALA))["t"] + "|" + c2["m"], "2-|0")
caso("InChI vacío no truena", capas(None), None)

# Carga y fórmula contra el GEM
caso("carga neta acetato = -1", carga_neta(ACETATO), -1)
caso("carga neta ácido acético = 0", carga_neta(ACIDO_ACETICO), 0)
caso("carga neta con /q (tetrametilamonio) = +1", carga_neta(TMA), 1)
caso("fórmula corregida acetato = C2H3O2", formula_corregida(ACETATO), "C2H3O2")
caso("acetato coincide con GEM C2H3O2", coincide_formula_gem(ACETATO, "C2H3O2"), True)
caso("ácido acético NO coincide con GEM C2H3O2", coincide_formula_gem(ACIDO_ACETICO, "C2H3O2"), False)
caso("fórmula del GEM con R es genérica", coincide_formula_gem(ACETATO, "C2H3O2R"), None)

from inchi_calidad import coincide_pesados_gem, formula_generica
PROTON = "InChI=1S/p+1"
caso("protón: sin fórmula", capas(PROTON)["formula"], "")
caso("protón: carga +1 (caso real que fallaba)", carga_neta(PROTON), 1)
caso("protón coincide con GEM H (fórmula)", coincide_formula_gem(PROTON, "H"), True)
caso("protón coincide con GEM H (pesados)", coincide_pesados_gem(PROTON, "H"), True)
caso("pesados: acetato vs GEM C2H3O2", coincide_pesados_gem(ACIDO_ACETICO, "C2H3O2"), True)
caso("pesados: carnitina C31 vs GEM C27 (caso real arachdcrn)",
     coincide_pesados_gem("InChI=1/C31H53NO4/c1-5-6/h1H", "C27H45NO4"), False)
caso("fórmula FULLR es genérica", formula_generica("C3H4O5PFULLRCO2"), True)
caso("fórmula normal no es genérica", formula_generica("C2H3O2"), False)

caso("carga con notación N*q (2*-1;+2) = 0 (hallazgo de la revisión)",
     carga_neta("InChI=1S/2CH3O.Ca/c2*1-2;/h2*1H3;/q2*-1;+2"), 0)
caso("carga con notación N*q sin contraion (2*-1) = -2", carga_neta("InChI=1S/2CH3O/c2*1-2/h2*1H3/q2*-1"), -2)

# Tipos de inconsistencia
caso("todas iguales = consistente", tipos_inconsistencia({M: L_ALA, P: L_ALA, C: L_ALA}), [])
caso("solo protonación", tipos_inconsistencia({M: ACETATO, P: ACIDO_ACETICO, K: ACIDO_ACETICO}), ["carga"])
caso("L vs D = estereo contradictoria", tipos_inconsistencia({M: L_ALA, P: D_ALA}), ["estereo_contradictoria"])
caso("estereo ausente = estereo incompleta", tipos_inconsistencia({M: L_ALA, P: ALA_SIN_ESTEREO}), ["estereo_incompleta"])
caso("alanina vs sarcosina = solo conectividad", tipos_inconsistencia({M: L_ALA, C: SARCOSINA}), ["conectividad"])
caso("acetato vs etanol = solo fórmula", tipos_inconsistencia({M: ACETATO, K: ETANOL}), ["formula"])
caso("dos tipos a la vez (fórmula y estereo)", tipos_inconsistencia({M: ETANOL, P: L_ALA, C: ALA_SIN_ESTEREO}),
     ["formula", "estereo_incompleta"])
caso("una sola fuente no se compara", tipos_inconsistencia({M: L_ALA}), None)

# Recomendación
r = recomendar({M: ACETATO, S: ACETATO, P: ACIDO_ACETICO, C: ACIDO_ACETICO, K: ACIDO_ACETICO}, carga_gem=-1)
caso("acetato: confianza alta", r["confianza"], "alta")
caso("acetato: se recomienda la forma del GEM", r["inchi"], ACETATO)
r = recomendar({M: D_ALA, S: D_ALA, P: L_ALA, C: L_ALA, K: L_ALA}, carga_gem=0)
caso("modelo con D y 3 bases con L: recomienda L", r["inchi"], L_ALA)
caso("modelo con D y 3 bases con L: confianza media", r["confianza"], "media")
caso("modelo con D y 3 bases con L: prioridad alta", r["prioridad_revision"], "alta")
r = recomendar({M: L_ALA, S: L_ALA}, carga_gem=0)
caso("solo fuentes del modelo: revisión manual", r["confianza"], "revision_manual")
r = recomendar({M: L_ALA, P: D_ALA}, carga_gem=0)
caso("empate 1 a 1: revisión manual", r["confianza"], "revision_manual")
caso("empate: sin InChI recomendado", r["inchi"], None)
r = recomendar({M: L_ALA, S: L_ALA, P: L_ALA, C: SARCOSINA, K: D_ALA}, carga_gem=0)
caso("MODELO+SMILES cuentan como una: sin mayoría → revisión", r["confianza"], "revision_manual")

caso("solo modelo: motivo", recomendar({M: L_ALA, S: L_ALA}, carga_gem=0)["motivo"], "solo fuentes del modelo")
caso("solo modelo: prioridad baja", recomendar({M: L_ALA, S: L_ALA}, carga_gem=0)["prioridad_revision"], "baja")
caso("conflicto sin mayoría: prioridad media", recomendar({M: L_ALA, P: D_ALA}, carga_gem=0)["prioridad_revision"], "media")
caso("sin ningún InChI: sin_datos", recomendar({}, carga_gem=0)["confianza"], "sin_datos")

caso("acetato: el GEM coincide", recomendar({M: ACETATO, P: ACIDO_ACETICO, C: ACIDO_ACETICO}, carga_gem=-1)["gem_coincide"], True)
caso("modelo D vs bases L: el GEM no coincide", recomendar({M: D_ALA, P: L_ALA, C: L_ALA}, carga_gem=0)["gem_coincide"], False)

# Similitud: la protonación no debe bajar el Tanimoto
from inchi_calidad import tanimoto_minimo, sin_carga
caso("sin_carga quita /p", sin_carga(ACETATO), ACIDO_ACETICO)
caso("sin_carga conserva la carga fija /q", sin_carga(TMA), TMA)
COLINA = "InChI=1S/C5H14NO/c1-6(2,3)4-5-7/h7H,4-5H2,1-3H3/q+1"
caso("Tanimoto de colina consigo misma = 1 (caso real que fallaba)", tanimoto_minimo({M: COLINA, P: COLINA}), 1.0)
caso("Tanimoto acetato vs ácido acético = 1", tanimoto_minimo({M: ACETATO, P: ACIDO_ACETICO}), 1.0)
caso("Tanimoto alanina vs sarcosina < 1", tanimoto_minimo({M: L_ALA, P: SARCOSINA}) < 1, True)

# Casos reales encontrados al revisar a mano (Fase 6, 3-nov-2026)
LAC_D_GEM = "InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)/p-1/t2-/m0/s1"     # ecoli lac_D: es (S) = L
LAC_D_BASES = "InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)/t2-/m1/s1"       # (R) = D
R5P_GEM = "InChI=1S/C5H11O8P/c6-3-2(1-12-14(9,10)11)13-5(8)4(3)7/h2-8H,1H2,(H2,9,10,11)/p-2/t2-,3-,4-,5+/m0/s1"
R5P_BASES = "InChI=1S/C5H11O8P/c6-3-2(1-12-14(9,10)11)13-5(8)4(3)7/h2-8H,1H2,(H2,9,10,11)/t2-,3-,4-,5?/m1/s1"
ARSENIOSO = "InChI=1S/AsH3O3/c2-1(3)4/h2-4H"
ARSENITO = "InChI=1S/AsO3/c2-1(3)4/q-3"
caso("lac_D: estereo contradictoria y carga", tipos_inconsistencia({M: LAC_D_GEM, P: LAC_D_BASES}),
     ["estereo_contradictoria", "carga"])
r = recomendar({M: LAC_D_GEM, S: LAC_D_GEM, P: LAC_D_BASES, C: LAC_D_BASES, K: LAC_D_BASES}, carga_gem=-1)
caso("lac_D: recomienda la forma (R) de las bases", r["inchi"], LAC_D_BASES)
caso("lac_D: prioridad alta", r["prioridad_revision"], "alta")
caso("r5p: estereo contradictoria (R/S opuestos, no solo el centro ?)",
     "estereo_contradictoria" in tipos_inconsistencia({M: R5P_GEM, P: R5P_BASES}), True)
caso("arsenito vs ácido arsenioso = solo carga (no fórmula)", tipos_inconsistencia({M: ARSENIOSO, C: ARSENITO, K: ARSENITO}), ["carga"])
r = recomendar({M: ARSENIOSO, S: ARSENIOSO, C: ARSENITO, K: ARSENITO}, carga_gem=0)
caso("arsenito: confianza alta", r["confianza"], "alta")
caso("arsenito: recomienda la forma con la carga del GEM (0)", r["inchi"], ARSENIOSO)

MTHF_REL = ("InChI=1/C20H25N7O6/c1-27-12(9-23-16-15(27)18(31)26-20(21)25-16)8-22-11-4-2-10(3-5-11)17(30)24-13(19(32)33)6-7-14(28)29"
            "/h2-5,12-13,22H,6-9H2,1H3,(H,24,30)(H,28,29)(H,32,33)(H4,21,23,25,26,31)/p-1/t12?,13-/s2")
MTHF_ABS = ("InChI=1S/C20H25N7O6/c1-27-12(9-23-16-15(27)18(31)26-20(21)25-16)8-22-11-4-2-10(3-5-11)17(30)24-13(19(32)33)6-7-14(28)29"
            "/h2-5,12-13,22H,6-9H2,1H3,(H,24,30)(H,28,29)(H,32,33)(H4,21,23,25,26,31)/p-1/t12?,13-/m1/s1")
caso("/s2 relativa: no se inventa la absoluta", "/m" in a_estandar(MTHF_REL), False)
caso("5mthf /s2 vs SMILES absoluto = incompleta, no contradictoria (caso real)",
     tipos_inconsistencia({M: MTHF_REL, S: MTHF_ABS}), ["estereo_incompleta"])

F6P_ANILLO = "InChI=1S/C6H13O9P/c7-2-6(10)5(9)4(8)3(15-6)1-14-16(11,12)13/h3-5,7-10H,1-2H2,(H2,11,12,13)/p-2/t3-,4-,5+,6?/m1/s1"
F6P_ABIERTA = "InChI=1S/C6H13O9P/c7-1-3(8)5(10)6(11)4(9)2-15-16(12,13)14/h4-7,9-11H,1-2H2,(H2,12,13,14)/t4-,5-,6-/m1/s1"
caso("f6p anillo vs abierta = anillo_cadena, no conectividad (caso real)",
     tipos_inconsistencia({M: F6P_ANILLO, P: F6P_ABIERTA})[0], "anillo_cadena")
caso("alanina vs sarcosina sigue siendo conectividad", tipos_inconsistencia({M: L_ALA, C: SARCOSINA}), ["conectividad"])

# El voto con estereoquímica incompleta
r = recomendar({M: L_ALA, S: L_ALA, P: ALA_SIN_ESTEREO, C: ALA_SIN_ESTEREO, K: ALA_SIN_ESTEREO}, carga_gem=0)
caso("GEM más específico que las bases: confianza alta", r["confianza"], "alta")
caso("GEM más específico: recomienda el del GEM", r["inchi"], L_ALA)
caso("GEM más específico: el GEM coincide", r["gem_coincide"], True)
caso("GEM más específico: motivo", r["motivo"], "coinciden; alguna fuente no define toda la estereoquímica")
r = recomendar({M: D_ALA, P: L_ALA, C: L_ALA, K: ALA_SIN_ESTEREO}, carga_gem=0)
caso("trampa: la fuente sin estereo no gana aunque sea compatible con todos", r["inchi"], L_ALA)
caso("trampa: el GEM no coincide", r["gem_coincide"], False)
caso("trampa: prioridad alta", r["prioridad_revision"], "alta")
r = recomendar({M: ALA_SIN_ESTEREO, S: ALA_SIN_ESTEREO, P: L_ALA, C: L_ALA}, carga_gem=0)
caso("GEM sin estereo y bases con L: recomienda L", r["inchi"], L_ALA)
caso("GEM sin estereo: el GEM coincide (es compatible)", r["gem_coincide"], True)
v = variables({M: L_ALA, P: ALA_SIN_ESTEREO, C: ALA_SIN_ESTEREO}, "C3H7NO2", 0)
caso("consenso de estructura con estereo incompleta = 1", v["consenso_estructura"], 1.0)
caso("GEM interno: MODELO vs SMILES L/D", variables({M: L_ALA, S: D_ALA}, "C3H7NO2", 0)["gem_interno"], "estereo_contradictoria")

# Utilidades
caso("base_id quita compartimento", base_id("atp[c]"), "atp")
caso("base_id con guion bajo de compartimento", base_id("atp_c"), "atp_c")
caso("InChI no estándar pasa a estándar", a_estandar("InChI=1/C2H4O2/c1-2(3)4/h1H3,(H,3,4)/f/h3H"), ACIDO_ACETICO)

ok = sum(1 for _, b, _, _ in casos if b)
for n, b, o, e in casos:
    print(("  ok    " if b else "  FALLA ") + n + ("" if b else f"  → salió {o!r}, esperaba {e!r}"))
print(f"\nResultado: {ok} bien, {len(casos) - ok} mal")
sys.exit(0 if ok == len(casos) else 1)
