# El analisis es en dos fases
## La primera es simple, 
filtlong --min_length 1300 --max_length 1700 01.fastq > 01.filtered.fastq

## Se activa el entorno de emu
conda activate emu

## Después se corre. Se sugiere que sea 
emu abundance *.filtered.fastq --db /data/databases/emu/ --threads 39 --type map-ont --output-dir emu_results
