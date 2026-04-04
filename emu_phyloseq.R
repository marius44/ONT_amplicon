# Emu2phyloseq
# Resultados para diversidad de 16s RNAobtenida desde EMU
# LMHS
# Marzo 2026

library(tidyverse)
library(phyloseq)
library(ggplot2)
library(ggsignif)
library(FSA)        # Para dunnTest
library(ggprism)
library(vegan)
library(RColorBrewer)
library(pheatmap)
#install.packages("ggsignif")
#install.packages("FSA")
#install.packages("ggprism")

setwd("/home/user/Escritorio/elena/emu_results")
getwd()
files <- list.files(pattern="rel-abundance.tsv", full.names=TRUE)
View(files)
files

# =========================================
# PHYLOSEQ PIPELINE PARA EMU 16S ONT SIN METADATA

# Cargar base de datos de taxonomía Emu taxonomy.tsv
tax_file <- read.delim("/data/databases/emu/taxonomy.tsv", check.names = FALSE)
tax_file
#View(tax_file)
#str(tax_file)

# Asegúrate de que las columnas sigan el formato requerido por phyloseq
# superkingdom, phylum, class, order, family, genus, species
tax_file <- tax_file %>%
  dplyr::select(superkingdom, phylum, class, order, family, genus, species) %>%
  as.matrix()

# Convertir todo a caracteres
tax_file <- apply(tax_file, 2, as.character)
tax_file
#View(tax_file)



# Cargar resultados de Emu (abundancias)
# Archivos rel-abundance.tsv en la carpeta "emu_results"
files <- list.files(pattern = "rel-abundance.tsv", full.names = TRUE)
files

# Combinar en una sola tabla quitando la terminación "_rel-abundance.tsv"

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
    names_from = sample,
    values_from = abundance,
    values_fill = 0
  ) %>%
  column_to_rownames("species")

str(emu)
View(emu)


# Preparar OTUs en matriz
emu_otu_mat <- as.matrix(emu)
emu_otu_mat <- apply(emu_otu_mat, 2, as.numeric)
emu_otu_mat

# IMPORTANTE: recuperar nombres de especies
rownames(emu_otu_mat) <- rownames(emu)

# Integrarlo
OTU <- otu_table(emu_otu_mat, taxa_are_rows = TRUE)


#Preparar taxonomía correctamente
tax_df <- as.data.frame(tax_file)
tax_df

# Asegurar que species está presente
tax_df$species <- as.character(tax_df$species)
tax_df

# Usar species como rownames
tax_df <- tax_df[!is.na(tax_df$species) & tax_df$species != "", ]
rownames(tax_df) <- tax_df$species

# Quitar columna redundante
tax_df$species <- NULL

# Convertir a matriz
tax_mat <- as.matrix(tax_df)
TAX <- tax_table(tax_mat)

# INTERSECCIÓN REAL (clave)
common_taxa <- intersect(rownames(OTU), rownames(TAX))
length(common_taxa)  # debería ser alto (~200+)

OTU2 <- prune_taxa(common_taxa, OTU)
OTU2
TAX2 <- prune_taxa(common_taxa, TAX)
TAX2

# Crear objeto phyloseq 
physeq <- phyloseq(OTU2, TAX2)
physeq

# Verificar que no este roto
table(tax_table(physeq)[, "superkingdom"])

# Filtrar taxa ausentes (opcional)
physeq_filt <- filter_taxa(physeq, function(x) sum(x) > 0, TRUE)
physeq_filt

# Diversidad alfa (Shannon y Simpson)
abund_mat <- as.matrix(otu_table(physeq_filt))
                           
# =========================================
# 1. NORMALIZAR (ya tienes abundancias relativas, pero aseguramos)
# =========================================
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))


# =========================================
# 2. BARPLOT APILADO (PHYLUM)
# =========================================
physeq_phylum <- tax_glom(physeq_rel, taxrank = "phylum")

df_phylum <- psmelt(physeq_phylum)

# Top phyla
top_phyla <- df_phylum %>%
  group_by(phylum) %>%
  summarise(MeanAbund = mean(Abundance)) %>%
  arrange(desc(MeanAbund)) %>%
  slice(1:10) %>%
  pull(phylum)

df_phylum$phylum <- ifelse(df_phylum$phylum %in% top_phyla,
                           df_phylum$phylum, "Other")

p_bar_phylum <- ggplot(df_phylum,
                       aes(x = Sample, y = Abundance, fill = phylum)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y = "Relative abundance", x = "Sample", fill = "Phylum")

p_bar_phylum
# =========================================
# 3. BARPLOT APILADO (familia)
# =========================================
physeq_family <- tax_glom(physeq_rel, taxrank = "family")
df_family <- psmelt(physeq_family)

top_family <- df_family %>%
  group_by(family) %>%
  summarise(MeanAbund = mean(Abundance)) %>%
  arrange(desc(MeanAbund)) %>%
  slice(1:15) %>%
  pull(family)

df_family$family <- ifelse(df_family$family %in% top_family,
                         df_family$family, "Other")

p_bar_family <- ggplot(df_family,
                      aes(x = Sample, y = Abundance, fill = family)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y = "Relative abundance", x = "Sample", fill = "family")

p_bar_family




# =========================================
# 3. BARPLOT APILADO (GÉNERO)
# =========================================
physeq_genus <- tax_glom(physeq_rel, taxrank = "genus")
df_genus <- psmelt(physeq_genus)

top_genus <- df_genus %>%
  group_by(genus) %>%
  summarise(MeanAbund = mean(Abundance)) %>%
  arrange(desc(MeanAbund)) %>%
  slice(1:15) %>%
  pull(genus)

df_genus$genus <- ifelse(df_genus$genus %in% top_genus,
                         df_genus$genus, "Other")

p_bar_genus <- ggplot(df_genus,
                      aes(x = Sample, y = Abundance, fill = genus)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y = "Relative abundance", x = "Sample", fill = "Genus")

p_bar_genus


# =========================================
# 4. ALPHA DIVERSIDAD (SHANNON)
# =========================================
alpha_div <- estimate_richness(physeq, measures = c("Shannon", "Simpson"))

alpha_div$Sample <- rownames(alpha_div)

p_alpha <- ggplot(alpha_div, aes(x = Sample, y = Shannon)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = round(Shannon, 2)), vjust = -0.5, size = 3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y = "Shannon diversity", x = "Sample")

p_alpha


# =========================================
# 5. BETA DIVERSIDAD (PCoA - BRAY)
# =========================================
# -----------------------------------------
# 1. Extraer nombres de muestras
# -----------------------------------------
samples <- sample_names(physeq)

# Crear metadata mínima
meta <- data.frame(Sample = samples)

# OPCIONAL: extraer grupo desde nombre (ej. EN, EP, ED)
meta$Group <- sub("^([A-Z]+).*", "\\1", meta$Sample)

rownames(meta) <- meta$Sample

# Añadir metadata al phyloseq
sample_data(physeq) <- sample_data(meta)

# -----------------------------------------
# 2. Transformar a abundancia relativa
# -----------------------------------------
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))

# -----------------------------------------
# 3. Distancia Bray-Curtis
# -----------------------------------------
dist_bc <- phyloseq::distance(physeq_rel, method = "bray")

# -----------------------------------------
# 4. PCoA
# -----------------------------------------
ord <- ordinate(physeq_rel, method = "PCoA", distance = dist_bc)

# -----------------------------------------
# 5. Extraer coordenadas
# -----------------------------------------
ord_df <- as.data.frame(ord$vectors)
ord_df$Sample <- rownames(ord_df)

# Unir metadata
ord_df <- left_join(ord_df, meta, by = "Sample")

# -----------------------------------------
# 6. Plot PRO
# -----------------------------------------
p_pcoa <- ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  
  # etiquetas con nombre de muestra
  geom_text(aes(label = Sample), vjust = -0.8, size = 3) +
  
  theme_classic() +
  labs(
    title = "PCoA (Bray-Curtis)",
    x = paste0("PCoA1 (", round(ord$values$Relative_eig[1] * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(ord$values$Relative_eig[2] * 100, 1), "%)")
  )

p_pcoa


# =========================================
# 6. HEATMAP (TOP 20 ESPECIES)
# =========================================
abund_mat <- as(otu_table(physeq_rel), "matrix")

top_taxa <- names(sort(rowSums(abund_mat), decreasing = TRUE))[1:20]
abund_top <- abund_mat[top_taxa, ]

# Etiquetas taxonómicas bonitas
tax <- as.data.frame(tax_table(physeq_rel))

labels <- paste(
  tax[top_taxa, "genus"],
  tax[top_taxa, "species"]
)

rownames(abund_top) <- labels

pheatmap(abund_top,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         fontsize_row = 8,
         main = "Top 20 species")








