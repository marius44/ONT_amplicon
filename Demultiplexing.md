# Demultiplex de lecturas
## Anotación importante:
Hacer esto solo si es el operador del secuenciador o bien si el operador no lo hizo
Para trabajar con las secuencias de amplicon de ONT, hay que dedemultiplexarlas:

## Usar sudo porque la carpeta esta bloqueado por la MinKnow
Eg ruta
```
~/Escritorio/minknow/data/ajolotes_murcielagos/ajolotes_murcielagos/20250328_1209_MN47942_FBA39123_f311223a/fastq_pass
```

```bash
sudo cat barcode17/*.fastq.gz > P53C.fastq.gz
```
