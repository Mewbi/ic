# 📝 CNMAC

LaTeX Template for submitting work/abstracts in category 1 of CNMAC.

## 📄 Shorter Abstract

A otimização de geometrias moleculares de reações químicas é um dos tópicos de interesse na química computacional. O presente trabalho tem como objetivo propor um método iterativo de otimização, baseado nos métodos de Newton e Secante para realizar a convergência de funções de Superfície de Energia Potencial (SEP). O método CBPD (Convergence Bases in Partial DerivativeS) proposto apresenta um algoritmo com menos custo computacional envolvido quando comparado com o Método de Newton, comumente utilizado em cenários de otimização. Para medir a performance do Método CBPD foi utilizado a função SEP da reação F + H2O -> FH + HO e seus resultados foram comparados com o Método de Newton para os mesmos cenários de convergência.

__Obs. It's wrote in Portuguese because was submited.__

## 🛠️ How to Build

In terminal execute

```bash
pdflatex resumo.tex
biber resumo
pdflatex resumo.tex
pdflatex resumo.tex
```

## 📂 Files

The two files able to edit are:

- `resumo.tex`: main file containing the LaTeX source code of the work/abstract.
- `ref.bib`: BibTeX/BibLaTeX format file containing bibliographic reference data.

The file `pssbmac.cls` contains the definitions of the default document class to be generated. This file should be used as is and should not be altered.

Other auxiliary files include:

- `resumo.pdf`: Generated PDF.
- `image.jpg`: Image used in the abstract.

## 🎓 Créditos

Developed for [SBMAC](https://www.sbmac.org.br/) - Society for Applied and Computational Mathematics.
