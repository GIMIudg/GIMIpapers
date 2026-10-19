"""
PASO 2: Validación del formulario con LLM-as-judge
Evalúa cada pregunta del schema en 4 dimensiones:
  - Claridad
  - Relevancia clínica
  - Completitud
  - Granularidad clínica (nueva)

Características:
  - Guardado incremental tras cada pregunta (tolerante a interrupciones)
  - Reintentos automáticos en caso de error de API
  - Compatible con schema_evaluation.json (306 preguntas)

Produce: form_validation_report.json  +  form_validation_report.csv
"""

import json
import time
import os
import pandas as pd
import anthropic

# ── Configuración ──────────────────────────────────────────────────────────────
API_KEY     = "tu_api_key_aquí"
FOLDER_PAHT = r"ruta_salida_aquí"
SCHEMA_PATH = "ruta_al_esquema_pregunta-contexto"
MODEL_NAME  = "claude-haiku-4-5-20251001"
MAX_RETRIES = 3
DELAY_SECONDS = 0
N_corridas = 9
client = anthropic.Anthropic(api_key=API_KEY)

for corrida in range(1,N_corridas+1):
    tag = f"c1_r{corrida:02d}"

    OUTPUT_JSON = os.path.join(FOLDER_PAHT, f"form_validation_report_v2_enhanced_{tag}.json")
    OUTPUT_CSV  = os.path.join(FOLDER_PAHT, f"form_validation_report_v2_enhanced_{tag}.csv")
    summary     = os.path.join(FOLDER_PAHT, f"section_summary_v2_enhanced_{tag}.json")
    
    print(f"\n{'='*60}")
    print(f"CORRIDA {corrida}/{N_corridas} — {tag}")
    print(f"{'='*60}")




    # ── Prompt de evaluación ───────────────────────────────────────────────────────
    EVAL_SYSTEM = """Eres un experto en diseño de instrumentos clínicos y en generación de datos sintéticos 
    para investigación en salud neurológica. Tu tarea es evaluar preguntas de una historia clínica estructurada.

    Para cada pregunta responde ÚNICAMENTE con un objeto JSON con esta estructura exacta, sin texto adicional:
    {
    "clarity_score": <1-5>,
    "clarity_comment": "<comentario breve>",
    "clinical_relevance_score": <1-5>,
    "clinical_relevance_comment": "<comentario breve>",
    "completeness_score": <1-5>,
    "completeness_comment": "<comentario breve>",
    "granularity_score": <1-5>,
    "granularity_comment": "<comentario breve>",
    "synthetic_generation_issues": "<problemas potenciales al generar datos sintéticos, o 'ninguno'>",
    "suggested_improvement": "<mejora concreta al enunciado o las opciones, o 'ninguna'>"
    }

    Escala general: 1=muy deficiente, 3=aceptable, 5=excelente

    CLARIDAD (1-5): ¿El enunciado es preciso y no ambiguo?
    1 = confuso o con múltiples interpretaciones
    3 = comprensible pero mejorable
    5 = claro, preciso y sin ambigüedad

    RELEVANCIA CLÍNICA (1-5): ¿Es pertinente para una historia clínica neurológica?
    1 = irrelevante para neurología
    3 = moderadamente relevante
    5 = directamente relevante para diagnóstico o seguimiento neurológico

    COMPLETITUD (1-5): ¿El universo de opciones cubre todos los casos importantes?
    1 = faltan opciones críticas
    3 = opciones aceptables pero con huecos notables
    5 = opciones exhaustivas, ningún caso clínico queda sin cubrir

    GRANULARIDAD CLÍNICA (1-5): ¿Las opciones son suficientemente específicas para distinguir casos clínicos?
    1 = opciones demasiado genéricas (ej: solo sí/no donde se necesita detalle)
    3 = granularidad aceptable pero subdivisible para mayor utilidad
    5 = nivel de detalle óptimo para análisis clínico y generación de fenotipos

    Nota: completitud alta no implica granularidad alta.
    Baja granularidad: tratamiento_diabetes → ['medicamento', 'dieta', 'otro']
    Alta granularidad: tratamiento_diabetes → ['metformina', 'insulina_basal', 'GLP-1', 'SGLT2', 'dieta']"""

    # ── Carga del schema ───────────────────────────────────────────────────────────
    with open(SCHEMA_PATH, encoding='utf-8') as f:
        schema = json.load(f)

    questions = schema['questions']
    total = len(questions)

    # ── Carga de progreso previo (si existe) ───────────────────────────────────────
    completed_names = set()
    results = []

    if os.path.exists(OUTPUT_JSON):
        with open(OUTPUT_JSON, encoding='utf-8') as f:
            results = json.load(f)
        completed_names = {r['name'] for r in results}
        print(f"Retomando desde pregunta anterior — {len(completed_names)}/{total} ya completadas")
    else:
        print(f"Iniciando evaluación — {total} preguntas")

    # ── Contexto clínico por sección iterativa ────────────────────────────────────
    REPEAT_CONTEXT = {
        'Registro de familiares con enfermedad':
            'Bloque iterativo para registrar antecedentes heredofamiliares. '
            'Se repite por cada familiar con enfermedad relevante. '
            'El objetivo es capturar el perfil genético-familiar del paciente.',

        'Registro de sustancias':
            'Bloque iterativo para registrar historial de consumo de sustancias. '
            'Se repite por cada sustancia consumida (tabaco, alcohol, drogas, etc.). '
            'Incluye patrón, frecuencia, vía de administración y temporalidad.',

        'Registro de animales':
            'Bloque iterativo para registrar convivencia con animales. '
            'Se repite por cada tipo de animal. Relevante para exposiciones '
            'alérgicas, zoonosis y factores ambientales.',

        'Registro de exposiciones':
            'Bloque iterativo para exposiciones ocupacionales o ambientales. '
            'Se repite por cada agente de exposición identificado. '
            'Relevante para factores de riesgo neurológico y sistémico.',

        '6.2.2 Registro de métodos anticonceptivos':
            'Bloque iterativo para historial de métodos anticonceptivos. '
            'Se repite por cada método usado. Relevante para antecedentes '
            'ginecológicos y farmacológicos de la paciente.',

        '8.1.2 Registro de cirugías':
            'Bloque iterativo para antecedentes quirúrgicos. '
            'Se repite por cada procedimiento previo. Incluye tipo de cirugía, '
            'periodo, complicaciones y comorbilidades perioperatorias.',

        '8.2.2 Registro de alergias':
            'Bloque iterativo para alergias a medicamentos o sustancias. '
            'Se repite por cada alérgeno identificado. Incluye tipo de reacción, '
            'gravedad y confirmación diagnóstica.',

        '8.3.2 Registro de transfusiones':
            'Bloque iterativo para antecedentes transfusionales. '
            'Se repite por cada evento de transfusión. Incluye producto, '
            'indicación, complicaciones y periodo.',

        '9.12.2 Detalle por tipo de cáncer (repita por cada tipo seleccionado)':
            'Bloque iterativo para caracterizar cada tipo de cáncer diagnosticado. '
            'Se repite por cada neoplasia. Incluye estadio, periodo de diagnóstico, '
            'presencia de metástasis y estado actual.',

        '9.15.1.a Otras enfermedades (agrega una por fila)':
            'Bloque iterativo para registrar enfermedades crónicas no contempladas '
            'en las secciones anteriores. Se repite por cada enfermedad adicional. '
            'Permite capturar comorbilidades atípicas o poco frecuentes.',

        '11.2 Registro de medicamentos':
            'Bloque iterativo para el tratamiento farmacológico actual o reciente. '
            'Se repite por cada medicamento. Incluye identificación, vía, dosis, '
            'frecuencia, adherencia, efectos adversos e indicación.',

        '13.2 Registrar estudio complementario':
            'Bloque iterativo para estudios diagnósticos mencionados en la consulta. '
            'Se repite por cada estudio (laboratorio, imagen, bioseñales, valoraciones). '
            'Permite registrar resultados y su mapeo a catálogos estándar (LOINC).',
    }

    # ── Función de construcción del prompt ────────────────────────────────────────
    def build_question_prompt(q):
        choices_text = "\n".join([f"  - {c['value']}: {c['label']}" for c in q['choices']]) \
                    if q['choices'] else "  (sin opciones — catálogo externo reducido)"
        relevant_text  = f"\nCondición de activación: {q['relevant']}" if q['relevant'] else ""

        # Contexto iterativo: solo el del repeat al que pertenece esta pregunta
        if q.get('in_repeat') and q.get('repeat_name'):
            repeat_ctx = REPEAT_CONTEXT.get(q['repeat_name'], '')
            in_repeat_text = (
                f"\nPertenece a sección iterativa: {q['repeat_name']}"
                + (f"\nContexto del bloque: {repeat_ctx}" if repeat_ctx else "")
            )
        else:
            in_repeat_text = ""

        return f"""Evalúa la siguiente pregunta de historia clínica neurológica:

    Sección: {q['section']}
    Variable: {q['name']}
    Enunciado: {q['label']}
    Tipo: {q['type']}{relevant_text}{in_repeat_text}
    Opciones de respuesta:
    {choices_text}"""

    # ── Evaluación ─────────────────────────────────────────────────────────────────
    pending = [q for q in questions if q['name'] not in completed_names]
    print(f"Preguntas pendientes: {len(pending)}\n")

    for i, q in enumerate(pending):
        evaluation = None

        for intento in range(MAX_RETRIES):
            try:
                response = client.messages.create(
                    model=MODEL_NAME,
                    max_tokens=1024,
                    system=EVAL_SYSTEM,
                    messages=[{
                        "role": "user",
                        "content": build_question_prompt(q)
                    }]
                )
                raw_text = response.content[0].text.strip()

                # Limpiar markdown si el modelo lo agrega
                if raw_text.startswith("```"):
                    raw_text = raw_text.split("```")[1]
                    if raw_text.startswith("json"):
                        raw_text = raw_text[4:]

                evaluation = json.loads(raw_text)
                break

            except json.JSONDecodeError:
                evaluation = {
                    "clarity_score": None, "clarity_comment": "Error de parseo",
                    "clinical_relevance_score": None, "clinical_relevance_comment": "Error de parseo",
                    "completeness_score": None, "completeness_comment": "Error de parseo",
                    "granularity_score": None, "granularity_comment": "Error de parseo",
                    "synthetic_generation_issues": raw_text[:200],
                    "suggested_improvement": "revisar manualmente"
                }
                break

            except Exception as e:
                if intento < MAX_RETRIES - 1:
                    wait = 30 * (intento + 1)
                    print(f"  Error en intento {intento+1}: {str(e)[:60]} — esperando {wait}s...")
                    time.sleep(wait)
                else:
                    evaluation = {
                        "clarity_score": None, "clarity_comment": str(e)[:100],
                        "clinical_relevance_score": None, "clinical_relevance_comment": str(e)[:100],
                        "completeness_score": None, "completeness_comment": str(e)[:100],
                        "granularity_score": None, "granularity_comment": str(e)[:100],
                        "synthetic_generation_issues": str(e)[:200],
                        "suggested_improvement": "revisar manualmente"
                    }

        # Construir resultado
        result = {
            'section':     q['section'],
            'name':        q['name'],
            'label':       q['label'],
            'n_choices':   len(q['choices']),
            'in_repeat':   q.get('in_repeat', False),
            'repeat_name': q.get('repeat_name'),
            **evaluation
        }
        results.append(result)

        # Score promedio para mostrar progreso
        scores = [evaluation[k] for k in
                ['clarity_score','clinical_relevance_score','completeness_score','granularity_score']
                if evaluation.get(k) is not None]
        avg = round(sum(scores)/len(scores), 1) if scores else None

        global_idx = len(completed_names) + i + 1
        print(f"[{global_idx}/{total}] {q['name']} — avg: {avg} "
            f"(C:{evaluation.get('clarity_score')} "
            f"R:{evaluation.get('clinical_relevance_score')} "
            f"Co:{evaluation.get('completeness_score')} "
            f"G:{evaluation.get('granularity_score')})")

        # ── Guardado incremental ───────────────────────────────────────────────
        with open(OUTPUT_JSON, 'w', encoding='utf-8') as f:
            json.dump(results, f, ensure_ascii=False, indent=2)

        time.sleep(DELAY_SECONDS)

    # ── CSV final ──────────────────────────────────────────────────────────────────
    df = pd.DataFrame(results)
    df.to_csv(OUTPUT_CSV, index=False, encoding='utf-8-sig')

    # ── Resumen ────────────────────────────────────────────────────────────────────
    df_ok = df[df['clarity_score'].notna()].copy()
    print(f"\n{'='*60}")
    print(f"Reporte guardado: {OUTPUT_JSON} y {OUTPUT_CSV}")
    print(f"Evaluadas: {len(df_ok)}/{total}")

    print(f"\nScores globales:")
    for col, label in [
        ('clarity_score',            'Claridad    '),
        ('clinical_relevance_score', 'Relevancia  '),
        ('completeness_score',       'Completitud '),
        ('granularity_score',        'Granularidad'),
    ]:
        if col in df_ok.columns:
            print(f"  {label}: {df_ok[col].mean():.2f} ± {df_ok[col].std():.2f}")

    print(f"\nResumen por sección:")
    for sec in df['section'].unique():
        sec_df = df_ok[df_ok['section'] == sec]
        if len(sec_df) == 0:
            continue
        c  = sec_df['clarity_score'].mean()
        r  = sec_df['clinical_relevance_score'].mean()
        co = sec_df['completeness_score'].mean()
        g  = sec_df['granularity_score'].mean()
        print(f"\n  {sec} ({len(sec_df)} preguntas)")
        print(f"    Claridad:     {c:.2f}/5")
        print(f"    Relevancia:   {r:.2f}/5")
        print(f"    Completitud:  {co:.2f}/5")
        print(f"    Granularidad: {g:.2f}/5")

    # ── Resumen por sección en JSON ────────────────────────────────────────────────
    section_summary = {}
    for sec in df['section'].unique():
        sec_df = df_ok[df_ok['section'] == sec]
        if len(sec_df) == 0:
            continue
        section_summary[sec] = {
            'n_questions':    int(len(sec_df)),
            'n_in_repeat':    int(sec_df['in_repeat'].sum()) if 'in_repeat' in sec_df.columns else 0,
            'clarity': {
                'mean': round(float(sec_df['clarity_score'].mean()), 2),
                'std':  round(float(sec_df['clarity_score'].std()), 2),
            },
            'relevance': {
                'mean': round(float(sec_df['clinical_relevance_score'].mean()), 2),
                'std':  round(float(sec_df['clinical_relevance_score'].std()), 2),
            },
            'completeness': {
                'mean': round(float(sec_df['completeness_score'].mean()), 2),
                'std':  round(float(sec_df['completeness_score'].std()), 2),
            },
            'granularity': {
                'mean': round(float(sec_df['granularity_score'].mean()), 2),
                'std':  round(float(sec_df['granularity_score'].std()), 2),
            },
            'avg_overall': round(float(sec_df[[
                'clarity_score', 'clinical_relevance_score',
                'completeness_score', 'granularity_score'
            ]].mean().mean()), 2),
            'n_needs_improvement': int(
                (sec_df[['completeness_score', 'granularity_score']] <= 2).any(axis=1).sum()
            ),
        }

    with open(summary, 'w', encoding='utf-8') as f:
        json.dump(section_summary, f, ensure_ascii=False, indent=2)

    print(f"\nResumen por sección guardado: section_summary.json")
