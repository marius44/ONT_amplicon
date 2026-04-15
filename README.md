# Análisis de secuencias 16s rRNA originadas desde ONT (Oxford Nanopore Technologies)
## En este repositorio encontrara 3 archivos:

* Demultiplexing.txt
* Analisis_emu_(pasarlo_a_R).md
* emu_phyloseq.R

El primero, Demultiplexing.txt, solo lo usa el operador del secuenciador. Si esta analizando secuencias que ya le fueron demultiplexadas, **ignorelo**.

El segundo, Analisis_emu_(pasarlo_a_R).md, es un pipeline corto que utiliza al programa Emu. Con este comienza propiamente el análisis de sus secuencias.

El tercero, emu_phyloseq.R, es un script de R. Este le generará gráficas de abundancia relativa (por default estan a nivel phylum, familia y género pero puede agregar niveles a su gusto), un PCoA (Bray-Curtis), un heatmap con los 20 taxones más abundantes y generará uan gráfica de barras y una tabla en *.txt con los índices de diversidad calculados (Shannon y Simpson). 


