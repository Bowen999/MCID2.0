# January

### Jan 4
* create a function to generate compounds' 2D structure, and properties, and return an html

* demo compounds info page is available: http://111.229.151.228:5000/mcid0120791

* manually checked 0 rxn plant database （correct）

### Jan 9

* hyperlink on the result page is clickable 

### Jan 10

**SQLite logic:**

* **database**: mcid2.db

* **tables**: kegg_0_rxn, kegg_1_rxn, kegg_2_rxn, plant_0_rxn, plant_1_rxn, plant_2_rxn

* **columns**: **mcid, name, exact_mass**, formula, InChIKey, SMILES, possible_reaction, substrate_id

**Merged plant database**



# Feb 

### Feb 7

* SMARTS and SMIRKS rules are figured out, 23/76 rxn rules are corrected

* change some UI (menu and footer)

### Feb 12

The compound info and structure have been download from KEGG
 List of compounds in KEGG: https://rest.kegg.jp/list/compound

If copyright is needed, we can say we use this list to exact SMILES from PubChem*

### Feb 15

* Change the layout of the FAQ

* Fix some UI bugs


### Feb 16

Find a problem, in some cases, need to identify this compound has a specific group, then apply rxn rule to this group 

function **predict_products_group** was written to handle this condition.

### Feb 20

4 New rxn rules



# Mar

### Mar 4
Creating a new function drop_group and down 4 rxn 

### Mar 5
down 5 rxn (all advanced loss group rules have down)

### Mar 6

**addition of glucose**
[Glycosidic bonds](https://en.wikipedia.org/wiki/Glycosidic_bond) (S-, N-, C-, and O-glycosidic bonds) can be added to a glucose group, but all compounds including S- can do this?

*example: Adenine >> Adenosine*

**addition of ATP/ADP/AMP**

[Adenylylation](https://en.wikipedia.org/wiki/Adenylylation#Mass_spectrometry)
all amino acids's -OH can be added to a AXP?  (confirmed: tyrosine and threonine, and sometimes serine).

### Mar 11
down R50 addition of glucose
#### Some ideas about RT pred and quantification
**Predicted RT** are not very accurate, which is due to the poor generalization ability of DL (Deep Learning). Taking [BioTransformer](https://doi.org/10.1093/bioinformatics/btae084) as an example, even though research results show its predictions are very accurate, but both its training and testing data employ the same experimental methods. However, different experimental methods can lead to significant changes in RT.

Using Cytosine, Cortisol ([from Evagelia C. Laiakis et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3286010/)) and Capsaicin (from the CAT Bridge dataset) as examples, their prediction results are far from the real outcomes (with Capsaicin being **702.6 vs 976.6381226**), and the results on CIL MS will be even worse.

Therefore, improving their prediction generalization capability is crucial. Theoretically, this requires users to **retrain** the model using data produced by their own methods, yet, the vast majority of researchers **lack the data necessary for training and validation**.

A potentially viable method is to use all MS1 obtained **single-matches** (if MS2 is present, MS2 should also matche with this putative identity) as the training set, to correct the level of prediction drift.



The same **"internal reference"** method could be considered for application in **quantitation**. Untargeted metabolomics can only perform relative quantification (not absolute quantification), which makes it difficult to use as a clinical indicator and for integration across studies. Therefore, if it is possible to identify a **very stable metabolite** (one that does not undergo concentration changes in the vast majority of cases, including illness) from a vast amount of data to serve as an internal standard, this could transform relative quantification into absolute quantification. This would greatly unlock the potential of metabolomics in clinical settings.



### Mar 14
Nucleosides:

[nucleoside + Phosphorylation —> nucleotide](https://deepblue.lib.umich.edu/bitstream/handle/2027.42/152758/cpnc1301.pdf?sequence=1#:~:text=Nucleosides%20can%20be%20phosphorylated%20by,a%20halogen%20or%20sulfonate%20ester)

### Mar 19
* down R52 addition of uridine, R54 addition of adenosine

### Mar 20
* Down deoxyribonucleoside
* addition of D-ribose rule is wrong, delete it
* addition of glutamate rule: [Phosphate group](https://www.kegg.jp/entry/R00239), [-NH2](https://www.kegg.jp/entry/R00494), [indole](https://www.kegg.jp/entry/R00674)

### Mar 21
* down R74 addition of glutamate (indole is ignored) 
### Mar 22
replace **R-(0)(=O)-NH2** by dropping a H in NH
https://www.kegg.jp/entry/R01149
https://www.kegg.jp/entry/R04011

https://www.kegg.jp/entry/R00577
https://www.kegg.jp/entry/R02959
**+, - all works**

### Mar 27

**Addition of CO2** 
```
[pattern 1](https://www.kegg.jp/entry/R00005)
[*:1][CH2]C(=O)O>>[*:1]C

[pattern 2](https://www.kegg.jp/entry/R00449) 
[*:1]C(~[*:2])C(=O)O>>[*:1]C[*:2]

[pattern 3](https://www.kegg.jp/entry/R11917): is very similar to the pattern 1 but has one more =O, which can be handle by addition of  O; [*:1]C(~[*:2])C(=O)O>>[*:1]C(=O)[*:2]

[pattern 4](https://www.kegg.jp/entry/R03508): using predict_products_group to do this
subgroup = 'cccccc-C(=O)O'  # Target OH group on a ring
smirks = '[*:1]-C(=O)O>>[*:1]'
 [*:1]-C(=O)O>>[*:1]-O
```
**Loss of CO2**

https://www.kegg.jp/entry/R00013

https://www.kegg.jp/entry/R04885

https://www.kegg.jp/entry/R04732



## Apr

### Apr 06 

* **Down loss of CO2**
* **loss or add H2 rules**:

N+ in 6 carbons ring (3 double bonds) are not recognized with SMARTS: https://www.kegg.jp/entry/R09508 or https://www.kegg.jp/entry/R00700

two double bounds >> one double bond:

https://www.kegg.jp/entry/R02965

https://www.kegg.jp/entry/R03025

https://www.kegg.jp/entry/R09095



### Apr 8

* loss of NH3:
```
pattern 1:
https://www.kegg.jp/entry/R09253
pattern 2:
https://www.kegg.jp/entry/R11328
pattern 3:
https://www.kegg.jp/entry/R03030
pattern 4:
https://www.kegg.jp/entry/R05758, https://www.kegg.jp/entry/R01544, https://www.kegg.jp/entry/R00996
pattern 5:
https://www.kegg.jp/entry/R00223
pattern 6 (-NH2 >> -OH):
https://www.kegg.jp/entry/R05170
```

* Addition of NH3 is 
```
pattern: [*:1]C=O>>[*:1](O)CN
pattern: [*:1]=[*:2]C([*:3])>>[*:1][*:2](-N)C(=[*:3])O
pattern: [*:1](=O)C[*:2]>>[*:1](-N)C(-O)[*:2]
```



loss O2/addition O2:
https://www.kegg.jp/entry/R04857
smirks = '[*:1]-O>>[*:1]-[NH2]'
[*:1]-[NH2]>>[*:1]-O



*addition of H20:*

=O <-> -NH2
https://www.kegg.jp/entry/R00729
https://www.kegg.jp/entry/R00688
pattern 2:
https://www.kegg.jp/entry/R00132
pattern 3:
https://www.kegg.jp/entry/R01234
pattern 4 (unavailable):
https://www.kegg.jp/entry/R01335



*loss O:*
-OH >> empty
https://www.kegg.jp/entry/R00713
https://www.kegg.jp/entry/R01296



### Apr 10

* **down R01 (dehydrogenation) & R02 (hydrogenation):**

```
Hydrogenation:
[*:1][O-]>>[*:1][O]
[*:1][n]>>[*:1][N+]
[*:1]-C=C-[*:2]>>[*:1]-C-C-[*:2]
[*:1]=[O]>>[*:1]O
[*:1]-[c]=[n]>>[*:1]-c-n
[C:1]#[C]>>[C:1]=C
[C:1]=[CH2]>>[C:1]-C

Hehydrogenation:
[*:1][O]>>[*:1][O-]
[*:1][n+]>>[*:1][N]
[*:1]-C-C-[*:2]>>[*:1]-C=C-[*:2]
[*:1]=[O]>>[*:1]O
[*:1]-c-n>>[*:1]-[c]=[n]
[C:1]=C>>[C:1]#[C]
[C:1]-C>>[C:1]=[CH2]

Reference: https://en.wikipedia.org/wiki/Hydrogenation
```
*can not do this: https://www.nad.com/nad-vs-nadcan not do this: https://www.nad.com/nad-vs-nad`h*



* **down R11 and R12 (loss/addition of water) (without considering dropping a water and connect)**

```
loss of water
smirks = '[P:1]-[OH0]>>[P:1][OH]'
https://www.kegg.jp/entry/R01234

addition of Water
smirks = '[*:1]=[*:2](-N)>>[*:1][*:2](=O)
https://www.kegg.jp/entry/R03887
```



### Apr 11

**R15&R16 loss of C2H4 addition of C2H4**
```
-C2H, -OH, -Br, Cl >> *-C2H4
substrate = 'CCCC(O)CCCC'
smirks = '[C:1]-[OH]>>[C:1]OCC'
mass_difference = '+'

substrate = 'CCCCOCC'
smirks = '[*:1]-[CH2]-[CH3]>>[*:1]'
mass_difference = '-'
```


**R05&R06 loss of NH/addition of NH**
```
substrate = 'CCCC(N=N)CCC'
smirks = '[*:1]([NH2])-[*:2]>>[*:1]-[*:2]'
smirks = '[N:1]=[N:2]>>[N:1]'
mass_difference = '-'

substrate = 'CCCC(N)CCC'
smirks = '[*:1]-[NH2]>>[*:1]-N=N'
mass_difference = '+'
```
**Add a PubChem Lite test demo**



### Apr 16

* **Down R07 loss of oxygen/R08 addition of oxygen**
* **Down R13 loss of CO/R14 addition of CO**
```
-COOH <-> -OH
```
* **Correct R44 addition of D-ribose**
```
substrate = 'C1=CN=CN=C1'
smirks = '[n:1]1cnccc1>>[n+:1]1(c2c(O)c(O)c(CO)o2)cnccc1'
smirks = '[n:1]1ccccc1>>[n+:1]1(c2c(O)c(O)c(CO)o2)ccccc' #R02137     

substrate = 'CC(CCNC1=NC=NC2=C1NC=N2)CO'
smirks = '[nH1:1]1c2[*:2][*:3][*:4][*:5]c2nc1>>[n:1]1(c3c(O)c(O)c(CO)o3)c2[*:2][*:3][*:4][*:5]c2nc1' #R08074 

substrate = 'C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N'
smirks = '[*:1]-P(=O)([OH1])[OH1]>>[*:1]-P(=O)(O)c1c(O)c(O)c(CO)o1' # R01051 
mass_difference = '+'
```
* **Correct R68/70/72 addition of AMP, ATP, ADP**

# May

### May 15

* complete unrelated compounds list
* 3 new rxn need to be added: Deamination, Transamination, loss/add CDP, and one possible: oxidative cyclization: https://www.kegg.jp/entry/R01481



# Rxn Rules

| R26  | sulfate conjugation | only considers 6 C ring (based on ChatGPT examples), but Hao's result only considers 5C ring. |
| ---- | ---- | ---- |
|      |      |      |





