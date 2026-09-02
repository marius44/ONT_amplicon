# ============================================================
# Resultados para diversidad de 16S RNA obtenida desde EMU
# LMHS - Marzo 2026
# Script corregido: limpieza de corchetes en tax_table y taxa_names
# ============================================================

library(tidyverse)
library(phyloseq)
library(ggplot2)
library(ggsignif)
library(FSA)
library(ggprism)
library(vegan)
library(RColorBrewer)
library(pheatmap)
library(ggsci)
library(officer)
library(rvg)
library(flextable)

setwd("/home/user/Escritorio/elena/emu_results")
getwd()

# ============================================================
# CARGAR BASE DE DATOS DE TAXONOMÍA EMU
# ============================================================

tax_file <- read.delim("/data/databases/emu/taxonomy.tsv", check.names = FALSE)

# Seleccionar columnas requeridas por phyloseq
tax_file <- tax_file %>%
  dplyr::select(superkingdom, phylum, class, order, family, genus, species) %>%
  as.matrix()

# Convertir todo a caracteres
tax_file <- apply(tax_file, 2, as.character)

# ============================================================
# CARGAR RESULTADOS DE EMU (ABUNDANCIAS)
# ============================================================

files <- list.files(pattern = "rel-abundance.tsv", full.names = TRUE)
files

# Combinar en una sola tabla
emu <- map_df(files, function(f) {
  sample_name <- basename(f) %>%
    sub("\\.tsv$", "", .) %>%
    sub("_rel-abundance$", "", .) %>%
    sub("\\.fastq\\.filtered$", "", .)
  
  read_tsv(f, show_col_types = FALSE) %>%
    select(species, abundance) %>%
    mutate(sample = sample_name)
}) %>%
  filter(!is.na(species) & species != "") %>%
  group_by(species, sample) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  pivot_wider(
    names_from  = sample,
    values_from = abundance,
    values_fill = 0
  ) %>%
  column_to_rownames("species")

# Preparar matriz OTU
emu_otu_mat           <- as.matrix(emu)
emu_otu_mat           <- apply(emu_otu_mat, 2, as.numeric)
rownames(emu_otu_mat) <- rownames(emu)

OTU <- otu_table(emu_otu_mat, taxa_are_rows = TRUE)

# ============================================================
# PREPARAR TAXONOMÍA
# ============================================================

tax_df          <- as.data.frame(tax_file)
tax_df$species  <- as.character(tax_df$species)
tax_df          <- tax_df[!is.na(tax_df$species) & tax_df$species != "", ]
rownames(tax_df) <- tax_df$species
tax_df$species  <- NULL

tax_mat <- as.matrix(tax_df)
TAX     <- tax_table(tax_mat)

# ============================================================
# INTERSECCIÓN Y OBJETO PHYLOSEQ
# ============================================================

common_taxa <- intersect(rownames(OTU), rownames(TAX))
length(common_taxa)

OTU2 <- prune_taxa(common_taxa, OTU)
TAX2 <- prune_taxa(common_taxa, TAX)

physeq <- phyloseq(OTU2, TAX2)
physeq

# ============================================================
# LIMPIEZA DE CORCHETES (CORRECCIÓN PRINCIPAL)
# Aplica AMBOS: taxa_names y el contenido de tax_table
# ============================================================

# 1. Limpiar los taxa_names (rownames del OTU)
taxa_names(physeq) <- taxa_names(physeq) %>%
  gsub("\\[|\\]", "", .) %>%
  gsub(" ", "_", .)

# 2. Limpiar el CONTENIDO de tax_table (genus, species, etc.)
#    Esto es lo que causaba los sufijos _2, _17 en el heatmap
clean_tax         <- apply(tax_table(physeq), 2, function(x) gsub("\\[|\\]", "", x))
tax_table(physeq) <- tax_table(clean_tax)

# Verificar integridad
table(tax_table(physeq)[, "superkingdom"])

# Filtrar taxa ausentes
physeq_filt <- filter_taxa(physeq, function(x) sum(x) > 0, TRUE)
physeq_filt

# ============================================================
# METADATA MÍNIMA (extraída desde nombres de muestras)
# ============================================================

samples <- sample_names(physeq)
meta    <- data.frame(Sample = samples)
meta$Group        <- sub("^([A-Z]+).*", "\\1", meta$Sample)
rownames(meta)    <- meta$Sample
sample_data(physeq) <- sample_data(meta)

# ============================================================
# 1. NORMALIZAR A ABUNDANCIA RELATIVA
# ============================================================

physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))

# Aplicar justo después de crear physeq_rel, antes de los tax_glom
# Rellenar NA y cadenas vacías en tax_table con "Unassigned"
tax_filled <- apply(tax_table(physeq_rel), 2, function(x) {
  x[is.na(x) | trimws(x) == ""] <- "Unassigned"
  return(x)
})
tax_table(physeq_rel) <- tax_table(tax_filled)

# ============================================================
# 2. BARPLOT APILADO — PHYLUM
# ============================================================

physeq_phylum <- tax_glom(physeq_rel, taxrank = "phylum")
df_phylum     <- psmelt(physeq_phylum)

top_phyla <- df_phylum %>%
  group_by(phylum) %>%
  summarise(MeanAbund = mean(Abundance)) %>%
  arrange(desc(MeanAbund)) %>%
  slice(1:10) %>%
  pull(phylum)

df_phylum$phylum <- ifelse(df_phylum$phylum %in% top_phyla, df_phylum$phylum, "Other")

p_bar_phylum <- ggplot(df_phylum, aes(x = Sample, y = Abundance, fill = phylum)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y = "Relative abundance", x = "Sample", fill = "Phylum")+
  scale_fill_igv()

p_bar_phylum

# ============================================================
# 3. BARPLOT APILADO — FAMILIA
# ============================================================

physeq_family <- tax_glom(physeq_rel, taxrank = "family")
df_family     <- psmelt(physeq_family)

top_family <- df_family %>%
  group_by(family) %>%
  summarise(MeanAbund = mean(Abundance)) %>%
  arrange(desc(MeanAbund)) %>%
  slice(1:15) %>%
  pull(family)

df_family$family <- ifelse(df_family$family %in% top_family, df_family$family, "Other")

p_bar_family <- ggplot(df_family, aes(x = Sample, y = Abundance, fill = family)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y = "Relative abundance", x = "Sample", fill = "Family")+
  scale_fill_igv()

p_bar_family

# ============================================================
# 4. BARPLOT APILADO — GÉNERO
# ============================================================

physeq_genus <- tax_glom(physeq_rel, taxrank = "genus")
df_genus     <- psmelt(physeq_genus)

top_genus <- df_genus %>%
  group_by(genus) %>%
  summarise(MeanAbund = mean(Abundance)) %>%
  arrange(desc(MeanAbund)) %>%
  slice(1:15) %>%
  pull(genus)

df_genus$genus <- ifelse(df_genus$genus %in% top_genus, df_genus$genus, "Other")

p_bar_genus <- ggplot(df_genus, aes(x = Sample, y = Abundance, fill = genus)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y = "Relative abundance", x = "Sample", fill = "Genus")+
  scale_fill_igv()

p_bar_genus

# ============================================================
# 5. DIVERSIDAD ALFA (SHANNON / SIMPSON)
# ============================================================

alpha_div         <- estimate_richness(physeq, measures = c("Shannon", "Simpson"))
alpha_div$Sample  <- rownames(alpha_div)
View(alpha_div)
alpha_div
write.table(alpha_div, file = "alpha.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

p_alpha <- ggplot(alpha_div, aes(x = Sample, y = Shannon)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = round(Shannon, 2)), vjust = -0.5, size = 3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y = "Shannon diversity", x = "Sample")

p_alpha

# ============================================================
# 6. BETA DIVERSIDAD (PCoA — BRAY-CURTIS)
# ============================================================

dist_bc <- phyloseq::distance(physeq_rel, method = "bray")
ord     <- ordinate(physeq_rel, method = "PCoA", distance = dist_bc)

ord_df         <- as.data.frame(ord$vectors)
ord_df$Sample  <- rownames(ord_df)
ord_df         <- left_join(ord_df, meta, by = "Sample")

ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(
    aes(label = Sample),    # ← label SOLO aquí, no en el aes() global
    vjust = -0.8, 
    size = 3,
    show.legend = FALSE     # ← evita que geom_text agregue entradas a la leyenda
  ) +
  theme_classic() +
  labs(
    title = "PCoA (Bray-Curtis)",
    x     = paste0("PCoA1 (", round(ord$values$Relative_eig[1] * 100, 1), "%)"),
    y     = paste0("PCoA2 (", round(ord$values$Relative_eig[2] * 100, 1), "%)")
  )

# ============================================================
# 7. HEATMAP (TOP 20 ESPECIES) — CORREGIDO
# ============================================================

abund_mat <- as(otu_table(physeq_rel), "matrix")
top_taxa  <- names(sort(rowSums(abund_mat), decreasing = TRUE))[1:20]
abund_top <- abund_mat[top_taxa, ]

# Los rownames YA son "Genus_species" — solo reemplazar _ por espacio
labels <- gsub("_", " ", top_taxa)

# make.unique como red de seguridad (no debería numerar si los taxa son distintos)
labels <- make.unique(labels)
rownames(abund_top) <- labels

pheatmap(
  abund_top,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale        = "row",
  fontsize_row = 8,
  main         = "Top 20 species"
)

###
# NO CORRER A PARTIR DE ESTA LINEA
# NO CORRER A PARTIR DE ESTA LINEA
# NO CORRER A PARTIR DE ESTA LINEA
# NO CORRER A PARTIR DE ESTA LINEA
# NO CORRER A PARTIR DE ESTA LINEA
# NO CORRER A PARTIR DE ESTA LINEA
# NO CORRER A PARTIR DE ESTA LINEA
# NO CORRER A PARTIR DE ESTA LINEA

###############################################################3

tax <- as.data.frame(tax_table(physeq_rel))

# Ver todas las filas de Blautia y Roseburia
tax[grepl("Blautia|Roseburia", rownames(tax)), c("genus", "species")]
##

# Ver nombres de columnas reales
colnames(as.data.frame(tax_table(physeq_rel)))
# Ver las primeras filas completas
head(as.data.frame(tax_table(physeq_rel)))


#
# Construir el plot por capas y ver en cuál aparece la "a"

# Solo puntos
ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  theme_classic()

# Agregar el texto
ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(aes(label = Sample), vjust = -0.8, size = 3) +
  theme_classic()
# Ver exactamente qué se está graficando como label
ord_df$Sample

# Buscar si algún Sample contiene "a" minúscula
ord_df[grepl("a", ord_df$Sample, ignore.case = FALSE), "Sample"]

# ============================================================
# EXPORTAR FIGURAS A POWERPOINT (OFFICER)
# ============================================================


# ============================================================
# 0. CREAR POWERPOINT
# ============================================================

ppt <- read_pptx()

# Función helper para añadir slide con ggplot (vectorial)
add_plot_slide <- function(ppt, title, plot_obj) {
  ppt %>%
    add_slide(layout = "Title and Content", master = "Office Theme") %>%
    ph_with(value = title, location = ph_location_type(type = "title")) %>%
    ph_with(dml(ggobj = plot_obj), location = ph_location_type(type = "body"))
}

# ============================================================
# 1. ABUNDANCIA RELATIVA
# ============================================================

ppt <- add_plot_slide(ppt, "Relative abundance — Phylum",  p_bar_phylum)
ppt <- add_plot_slide(ppt, "Relative abundance — Family",  p_bar_family)
ppt <- add_plot_slide(ppt, "Relative abundance — Genus",   p_bar_genus)

# ============================================================
# 2. ALPHA DIVERSITY (GRÁFICA)
# ============================================================

ppt <- add_plot_slide(ppt, "Alpha diversity (Shannon)", p_alpha)

# ============================================================
# 3. ALPHA DIVERSITY (TABLA REAL, SIN CAMBIOS)
# ============================================================

ft_alpha <- flextable(alpha_div)
ft_alpha <- autofit(ft_alpha)

ppt <- ppt %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with("Alpha diversity table", location = ph_location_type(type = "title")) %>%
  ph_with(ft_alpha, location = ph_location_type(type = "body"))

# ============================================================
# 4. PCoA (USAR EXACTAMENTE TU MISMO PLOT)
# ============================================================

p_pcoa <- ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(
    aes(label = Sample),
    vjust = -0.8,
    size = 3,
    show.legend = FALSE
  ) +
  theme_classic() +
  labs(
    title = "PCoA (Bray-Curtis)",
    x = paste0("PCoA1 (", round(ord$values$Relative_eig[1] * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(ord$values$Relative_eig[2] * 100, 1), "%)")
  )

ppt <- add_plot_slide(ppt, "PCoA (Bray-Curtis)", p_pcoa)

# ============================================================
# 5. HEATMAP (SIN CAMBIAR TU OBJETO)
# pheatmap NO es ggplot → se exporta como imagen
# ============================================================

# Guardar heatmap EXACTAMENTE como lo generaste
png("heatmap.png", width = 3000, height = 2500, res = 300)

pheatmap(
  abund_top,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale        = "row",
  fontsize_row = 8,
  main         = "Top 20 species"
)

dev.off()

# Insertar imagen
ppt <- ppt %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with("Heatmap — Top 20 species", location = ph_location_type(type = "title")) %>%
  ph_with(external_img("heatmap.png", width = 8, height = 6),
          location = ph_location_type(type = "body"))

# ============================================================
# 6. GUARDAR
# ============================================================

print(ppt, target = "Microbiome_EMU_results.pptx")

