# Análisis de amplicones 16s rRNA obtenido por ONT 
## Filtrado de secuencias del tamaño del 16s (1500 bp) 
filtlong --min_length 1300 --max_length 1700 *.fastq > *.filtered.fastq

## Se activa el entorno de emu
conda activate emu

## Después se corre. Se sugiere que sea 
emu abundance *.filtered.fastq --db /data/databases/emu/ --threads 39 --type map-ont --output-dir emu_results

## Se desactiva el entorno de emu
conda deactivate

