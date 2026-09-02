# Análisis de amplicones 16s rRNA obtenido por ONT 

## Filtrado de secuencias del tamaño del 16S rRNA (1500 bp). Filtrado con filtlong
### Esto debe hacerse por cada muestra:
filtlong --min_length 1300 --max_length 1700 MUESTRA1.fastq > MUESTRA1.filtered.fastq

## Se activa el entorno de emu
conda activate emu

## Después se corre emu con las secuencias filtradas. Se hace uno por uno, con los archivos que terminan en .filtered.fastq  
emu abundance MUESTRA1.filtered.fastq --db /data/databases/emu/ --keep-counts --threads 39 --type map-ont --output-dir emu_results

## Esto generará una serie de archivos que terminan en .filtered_rel-abundance.tsv. Abra uno para revisar el contenido
cd emu_results
less Muestra.filtered_rel-abundance.tsv
## Para salir de esa pantalla presione la tecla q

## Se desactiva el entorno de emu
conda deactivate

## Ahora se puede pasar al script de R _emu_phyloseq.R_
