# InfluenzaClassifier
Avian Influenza Classifier based on EU reference lab classification method.

Usage
```
perl InfluenzaClassifier/MAIN.pl SingleGenotype.fa
```
Read in a single fasta containing 8 Influenza segments. Header has to contain "|1|", "|2|", "|3|", etc.
Output is stdout

References were trimmed to be 95% of diversity (reduced sequence redundancy).
