# 📝 Notes

LaTeX notes explaining the project step by step.

## 🛠️ How to Build

In terminal execute

```bash
pdflatex 00_main.tex
bibtex 00_main
pdflatex 00_main.tex
pdflatex 00_main.tex
```

## 📂 Files & Directories

Important files:

- `00_main.tex`: main file containing the LaTeX entrypoint to build.
- `ref.bib`: BibTeX/BibLaTeX format file containing bibliographic reference data.

__Any other `.tex` files are imported by `00_main.tex`__

Directories:

- `reference/`: reference LaTex files.
- `figures/`: figures file, usually created using [TikZ](https://tikz.net/)
