# Publicar este paquete

El contenido de esta carpeta está preparado para colocarse en la **raíz del repositorio MethMap**. La preparación local no ha creado commits ni publicado cambios en GitHub.

Incluye código, README, un esquema del método, capturas de visualización global y regional aportadas por el autor, plantillas, datos sintéticos, pruebas y el YAML de dependencias. La copia de archivo de los scripts originales queda fuera de la carpeta pública. Los BAM, referencias hg38, archivos de cohortes reales, resultados de ejecución y registros están excluidos del paquete de distribución.

## Repositorio vacío

Desde una carpeta independiente que contenga estos archivos en su raíz:

```bash
git init -b main
git add README.md CITATION.cff environment.yml .gitignore methmap.py \
  full.sh adding.sh scripts config docs examples tests
git diff --cached --stat
git commit -m "Package MethMap tissue-specific methylation reference pipeline"
git remote add origin https://github.com/Francisco100Iniesta/MethMap.git
git push -u origin main
```

## Repositorio con contenido previo

Clona el repositorio, crea una rama y copia los archivos de este paquete sobre la copia de trabajo. Revisa el diff antes de confirmar y subir. No es necesario forzar el historial remoto.

```bash
git clone https://github.com/Francisco100Iniesta/MethMap.git MethMap-publicacion
cd MethMap-publicacion
git switch -c codex/package-methmap
# Copia aquí el contenido del paquete preparado.
git status --short
git diff
```

Antes de etiquetar una versión, completa la prueba BAM→BigWig indicada en [VALIDATION.md](VALIDATION.md). El archivo `CITATION.cff` identifica al autor y el repositorio; no se ha elegido automáticamente una licencia de distribución para el código propio.
