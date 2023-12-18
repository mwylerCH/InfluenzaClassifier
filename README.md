# InfluenzaClassifier
Avian Influenza Classifier based on EU reference lab classification method.

Usage
```
perl InfluenzaClassifier/MAIN.pl SingleGenotype.fa
```
Read in a single fasta containing 8 Influenza segments. Header has to contain "|1|", "|2|", "|3|", etc.
Output is stdout

References were trimmed to be 95% of diversity (reduced sequence redundancy).

Dependencies:
- [mafft](https://mafft.cbrc.jp/alignment/software/linux.html)
- [fasttree](http://www.microbesonline.org/fasttree/)
- [gotree](https://github.com/evolbioinfo/gotree)
- R Packages [data.table](https://cran.r-project.org/web/packages/data.table/index.html) and [tydiverse](https://cran.r-project.org/web/packages/tidyverse/index.html)
