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



# Rxn Rules

| R26  | sulfate conjugation | only considers 6 C ring (based on ChatGPT examples), but Hao's result only considers 5C ring. |
| ---- | ---- | ---- |
|      |      |      |

