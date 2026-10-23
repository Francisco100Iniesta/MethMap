# Estado de validación

Validación local realizada el **22 de octubre de 2026** con datos sintéticos.

## Ejecutado

- Construcción de la demo desde bedMethyl: seis muestras, dos cromosomas sintéticos y un cromosoma sin observaciones.
- Ampliación con una séptima muestra y una posición nueva, conservando la referencia anterior.
- Cuantiles conocidos: `0,20,40,60,80,100` producen p5=5, p50=50, p95=95 e IQR=50.
- Pruebas de cobertura válida, porcentajes 0–100, gzip, coordenadas, datos ausentes y valores cero.
- Rechazo de posiciones duplicadas, muestras repetidas, cambios de FASTA y carpetas de salida ya existentes.
- Construcción y ampliación desde una copia del código en una ruta con espacios.
- Generación y revisión visual del esquema del flujo y de la figura de la demo sintética.

Actualización de documentación del **23 de octubre de 2026**: el README conserva el esquema inicial y sustituye la figura sintética por dos capturas aportadas por el autor, una global y otra regional. Estas capturas ilustran el uso de las pistas; no son resultados de las pruebas sintéticas descritas aquí.

Resultado de la suite: **12 pruebas descubiertas, 11 superadas y 1 omitida** —la integración BAM→BigWig—. También se comprobaron la sintaxis del YAML, Python, R y los lanzadores Bash, la previsualización `--dry-run` y la resolución del directorio de trabajo del lanzador SLURM. No se envió un trabajo a un clúster.

Entorno local utilizado: Python 3.13.13, R 4.6.0 y data.table 1.18.6.1. Este entorno preexistente permite verificar la ruta bedMethyl→RDS→estadísticas→bedGraph. No equivale a haber resuelto e instalado el YAML completo.

## Pendiente de ejecutar con todas las dependencias

La prueba `BamIntegrationTests.test_modified_bam_to_bigwig` construye dos BAM sintéticos con etiquetas MM/ML, ejecuta la pipeline y comprueba los cuantiles y la generación de 23 archivos BigWig. Se omite cuando faltan Modkit, Samtools o `bedGraphToBigWig`.

La instalación de un entorno temporal con Mamba se intentó, pero los repositorios Conda no respondieron —fallos de resolución y de conexión—. Por tanto, **no se ha confirmado la resolución integral de `environment.yml` ni la ejecución BAM→BigWig en esta máquina**. El nombre `ont-modkit`, la versión 0.6.4 de su receta y las opciones de la interfaz se contrastaron con documentación primaria.

Para completar la comprobación en un equipo con acceso a Conda:

```bash
mamba env create -f environment.yml
conda activate methmap
python methmap.py doctor
python -m unittest discover -s tests -v
```

Con el entorno completo, la prueba BAM→BigWig debe ejecutarse y no aparecer como omitida. Para congelar las versiones resueltas de una ejecución concreta:

```bash
conda list --explicit > environment-linux-64.lock.txt
```

No se ha vuelto a ejecutar la cohorte humana ni se ha medido tiempo, memoria o concordancia sobre el genoma completo. Las pruebas pequeñas verifican funciones y casos límite; no constituyen validación biológica o de rendimiento a escala hg38.
