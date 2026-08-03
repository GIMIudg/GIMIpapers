# ============================================================
# PLOS ONE — ESPECIFICACIONES DE EXPORTACIÓN
# ============================================================
# Formato:     TIFF (compresión LZW)
# Resolución:  300 dpi
# Ancho:       789–2250 px @ 300 dpi  →  2.63–7.5 in
# Alto máx.:   2625 px @ 300 dpi      →  8.75 in
# Fuente:      Times New Roman (ya usada), 8–12 pt
# Color:       RGB
# Fondo:       blanco
# ============================================================

# === Librerías ===
if(!require(ggplot2))  install.packages("ggplot2",  dependencies=TRUE)
if(!require(scales))   install.packages("scales",   dependencies=TRUE)
if(!require(dplyr))    install.packages("dplyr",    dependencies=TRUE)

library(ggplot2)
library(scales)
library(dplyr)

# ── Directorio de salida: misma carpeta que el script ──────
OUTPUT_DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
# Si se ejecuta fuera de RStudio, usar directorio de trabajo:
# OUTPUT_DIR <- getwd()

# ── Helper: guardar como TIFF LZW a 300 dpi ────────────────
# w, h en pulgadas. Límites PLOS: w ≤ 7.5, h ≤ 8.75
save_plos <- function(plot_obj, filename, w, h) {
  path <- file.path(OUTPUT_DIR, filename)
  tiff(
    filename    = path,
    width       = w,
    height      = h,
    units       = "in",
    res         = 300,
    compression = "lzw",
    bg          = "white",
    type        = "cairo"      # RGB de alta calidad; quitar si no está disponible
  )
  print(plot_obj)
  dev.off()
  size_mb <- file.info(path)$size / 1048576
  message(sprintf("  ✅  %s  (%.2f MB)", filename, size_mb))
}

# ── Detección automática de Arial (igual que en el script Python) ──
if(!require(systemfonts)) install.packages("systemfonts", dependencies=TRUE)
library(systemfonts)

arial_available <- nrow(system_fonts()[system_fonts()$family == "Arial", ]) > 0
FONT_FAMILY <- ifelse(arial_available, "Arial", "DejaVu Sans")
message(sprintf("  Fuente: %s%s", FONT_FAMILY,
                ifelse(arial_available, "  ✔", "  ⚠ Arial no encontrada — usando DejaVu Sans")))

# ── Tema base PLOS One ──────────────────────────────────────
# Fuente Times, tamaños 8–12 pt, sin elementos innecesarios
theme_plos <- function() {
  theme_minimal(base_size = 10, base_family = FONT_FAMILY) +
    theme(
      plot.title         = element_text(size = 12, face = "bold",
                                        family = FONT_FAMILY, margin = margin(b = 6)),
      axis.title         = element_text(size = 11, family = FONT_FAMILY),
      axis.text          = element_text(size = 10, family = FONT_FAMILY),
      legend.title       = element_text(size = 10, family = FONT_FAMILY),
      legend.text        = element_text(size = 10, family = FONT_FAMILY),
      panel.grid.minor   = element_blank(),
      plot.background    = element_rect(fill = "white", color = NA),
      panel.background   = element_rect(fill = "white", color = NA)
    )
}

# === Leer CSV ===
df <- read.csv("/Users/eduardoruiz/Documents/GitHub/Precision-Oncology-for-Breast-Cancer-Diagnosis/src/GEMs construction MATLAB/GEMS_Exploratory_Analysis/Resumen_Modelos_ThermoKernel.csv")

# Colores
color_reacciones <- "#1f77b4"
color_similitud  <- "#2ca02c"

# ============================================================
# FIG 1 — Número de reacciones por modelo
# ============================================================
fig1 <- ggplot(df, aes(x = reorder(Model, NumReactions), y = NumReactions)) +
  geom_col(fill = color_reacciones, width = 0.7) +
  labs(
    title = "Number of reactions per model (ascending)",
    x     = "Models",
    y     = "Number of reactions"
  ) +
  theme_plos() +
  theme(
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    panel.grid.major.x = element_blank()
  )

save_plos(fig1, "Fig1_reactions_per_model.tif", w = 5.2, h = 3.5)

# ============================================================
# FIG 2 — Similaridad promedio por modelo
# ============================================================
fig2 <- ggplot(df, aes(x = reorder(Model, AverageSimilarity), y = AverageSimilarity)) +
  geom_col(fill = color_similitud, width = 0.7) +
  labs(
    title = "Average similarity per model (ascending)",
    x     = "Models",
    y     = "Average similarity"
  ) +
  theme_plos() +
  theme(
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    panel.grid.major.x = element_blank()
  )

save_plos(fig2, "Fig2_similarity_per_model.tif", w = 5.2, h = 3.5)

# ============================================================
# FIG 3 — Histograma: Número de reacciones
# ============================================================
fig3 <- ggplot(df, aes(x = NumReactions)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 15, fill = "lightgreen", color = "black", linewidth = 0.3) +
  geom_density(color = "darkgreen", linewidth = 0.8) +
  geom_vline(xintercept = mean(df$NumReactions),
             color = "red", linetype = "dashed", linewidth = 0.8) +
  labs(
    title = "Distribution of number of reactions",
    x     = "Number of reactions",
    y     = "Density"
  ) +
  theme_plos()

save_plos(fig3, "Fig3_histogram_reactions.tif", w = 5.2, h = 3.5)

# ============================================================
# FIG 4 — Histograma: Similaridad promedio
# ============================================================
fig4 <- ggplot(df, aes(x = AverageSimilarity)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 15, fill = "lightcoral", color = "black", linewidth = 0.3) +
  geom_density(color = "darkred", linewidth = 0.8) +
  geom_vline(xintercept = mean(df$AverageSimilarity),
             color = "blue", linetype = "dashed", linewidth = 0.8) +
  labs(
    title = "Distribution of average similarity",
    x     = "Average similarity",
    y     = "Density"
  ) +
  theme_plos()

save_plos(fig4, "Fig4_histogram_similarity.tif", w = 5.2, h = 3.5)

# ============================================================
# FIG 5 — Scatter: NumReactions vs AverageSimilarity
# ============================================================
fig5 <- ggplot(df, aes(x = NumReactions, y = AverageSimilarity, label = Model)) +
  geom_point(aes(color = AverageSimilarity, size = NumReactions), alpha = 0.85) +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  scale_size_continuous(range = c(1.5, 5)) +   # ✅ tamaño de puntos legible a 300 dpi
  labs(
    title = "Reactions vs. Similarity (1226 Metabolic Models)",
    x     = "Number of reactions",
    y     = "Average similarity"
  ) +
  theme_plos() +
  theme(plot.margin = margin(t = 4, r = 6, b = 6, l = 4))

save_plos(fig5, "Fig5_scatter_reactions_similarity.tif", w = 5.2, h = 4.0)

# ============================================================
# RESUMEN
# ============================================================
message("\n", strrep("=", 50))
message("  TODAS LAS FIGURAS GENERADAS")
message(strrep("=", 50))
message("  DPI      : 300")
message("  Formato  : TIFF (LZW)")
message("  Fuente   : Times New Roman, 10–12 pt")
message("  Color    : RGB")
message("  Archivos : Fig1–Fig5.tif")
message(strrep("=", 50))