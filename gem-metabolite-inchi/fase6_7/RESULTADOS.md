# Fases 6 y 7: resultados (corrida del 4-nov-2026)

Cifras de `DataMining_GEM.ipynb` con su configuración: ecoli_core_model, iDopaNeuroC y Recon3D_301.
Corriendo `fase6_7_analisis.py` solo, con smallModel incluido, los totales suben en 22 metabolitos y
los grupos no cambian. Excel completo: `resultados/fase6_7/tabla_final_fase7.xlsx`; figuras `fig1…fig5`.

## Unidad de análisis
Las 9,715 filas de la Fase 3 son **4,858 metabolitos** al fusionar compartimentos (`atp[c]` = `atp[m]`;
se comprobó que ningún metabolito cambia de InChI entre compartimentos). Recon3D pasa de 8,399 a 4,140.

| | Metabolitos | Sin InChI | 1 fuente | Comparables (≥2) | InChI idéntico | Con inconsistencias |
|---|---|---|---|---|---|---|
| ecoli_core_model | 54 | 4 | 0 | 50 | 9 | 41 |
| iDopaNeuroC | 664 | 52 | 1 | 611 | 166 | 445 |
| Recon3D_301 | 4,140 | 1,203 | 12 | 2,925 | 1,229 | 1,696 |
| **Total** | **4,858** | **1,259** | **13** | **3,586** | **1,404** | **2,182** |

## Fase 6.1: qué tipo de inconsistencia aparece (de 3,586 comparables; uno puede tener varios)
| Tipo | Metabolitos | Lectura química |
|---|---|---|
| Carga / protonación | 1,294 | Misma molécula a otro pH: el GEM usa la forma a pH fisiológico y las bases la neutra. No es error. |
| Estereoquímica incompleta | 1,222 | Una fuente no define un centro quiral que otra sí. Sobre todo el InChI del GEM (18.5 %, InChIs `/s2` de Recon3D) y PubChem (15.9 %). |
| Estereoquímica contradictoria | 137 | Configuración opuesta (R vs S). Es el error más serio: puede ser el enantiómero equivocado. |
| Fórmula | 61 | Distinto número de átomos pesados: otra molécula, polímeros/poliglutamatos representados distinto. |
| Anillo / cadena | 32 | El mismo azúcar en forma abierta en una fuente y cíclica en otra (p. ej. fructosa-6-P). |
| Hidrógenos | 8 | Tautómeros. |
| Conectividad | 3 | Isómeros constitucionales distintos. |

**La diferencia de carga entre el GEM y las bases no es un problema de anotación:** al neutralizar las
estructuras, la estructura coincide. Eso explica el 19–39 % de coincidencia en carga de la Fase 3.

## Fase 6.2: clustering (K-means, k = 8; 2,911 metabolitos únicos con ≥2 fuentes)
Variables: nº de bases externas, consenso del InChI completo, consenso de estructura (compatible en
estereo), consenso de la capa principal, Tanimoto mínimo (estructura neutralizada), fracción con la carga
del GEM y fracción con los átomos pesados de la fórmula del GEM. Se eligió k = 8: silhouette 0.671 y
Davies-Bouldin 0.583, los mejores entre k ≥ 3; K-means superó a Ward en todos los k.

| Grupo | n | Qué lo caracteriza | Qué dice de la calidad del GEM |
|---|---|---|---|
| Alta consistencia | 381 | ~2.8 bases externas, todas coinciden, incluso en carga. Moléculas pequeñas y sin grupos ionizables (agua, CO2, etanol). | Bien anotados. |
| Misma estructura, difiere solo la protonación | 603 | ~2.8 bases, estructura idéntica, solo 42 % con la carga del GEM. NAD, CoA, PEP. | Bien anotados; la diferencia es de pH. |
| Estereoquímica contradictoria | 110 | 97 % con R/S opuesto entre fuentes. | Revisar: posible enantiómero equivocado. |
| Discrepancias de constitución | 70 | Fórmula o conectividad distintas; Tanimoto 0.35. Azúcares abierta/cíclica, poliglutamatos, algunos iones. | Mezcla de representación y errores reales. |
| Carga/fórmula del GEM ionizadas, estructura neutra | 292 | Ninguna fuente tiene la carga del GEM, ni su propio InChI. Casi todo Recon3D. | Inconsistencia **interna** del GEM: los campos de fórmula y carga no corresponden a su estructura. |
| La fórmula del GEM no corresponde a su estructura | 46 | Los átomos pesados no coinciden (ej. `arachdcrn`: GEM C27, InChI C31). | Error de anotación en el GEM. |
| Estereoquímica incompleta, poca evidencia externa | 740 | Casi sin bases externas; el InChI del GEM trae estereo relativa (`/s2`). | El GEM no fija la configuración absoluta. |
| Información insuficiente (solo el GEM) | 669 | ~0.3 bases externas; el GEM es consistente consigo mismo. | No se puede validar. |

Fuera del clustering, por regla: 1,259 sin ningún InChI, 13 con una sola fuente y 30 con fórmula genérica
en el GEM (clases de lípidos con R).

## Fase 7: tabla final
| Estado de la anotación | Metabolitos |
|---|---|
| Correctamente anotado | 1,924 |
| Inconsistente: carga/fórmula del GEM no coinciden con la estructura | 238 |
| Inconsistente: el GEM difiere del consenso externo (**prioridad alta**) | 23 |
| Inconsistente: fuentes en conflicto (revisión manual) | 89 |
| Sin evidencia externa suficiente | 2,584 |

Confianza de la recomendación: alta 1,512 · media 673 · revisión manual 1,414 · sin datos 1,259.

Casos comprobados a mano entre los prioritarios: `lac_D` de ecoli (el campo InChI del GEM es L-lactato y
su propio SMILES es D-lactato), `r5p` de ecoli (R/S opuestos en 3 centros frente a las 3 bases) y `galam`
de Recon3D (las bases dan galactosamina-6-fosfato en lugar de galactosamina).

## Limitaciones a declarar
- Iones de un átomo (K⁺, Cl⁻): la huella Morgan no sirve para medir su similitud; pocos casos.
- La estereo relativa (`/s2`) se trata como "absoluta no definida".
- El análisis usa la versión de RDKit de la Mac (2025.09.2); Colab trae otra y pequeñas diferencias son posibles.
- HMDB sigue sin consultarse (decisión pendiente de la Fase 3).
