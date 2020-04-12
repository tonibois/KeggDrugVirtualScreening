# KeggDrugVirtualScreening
A software for analysis and fast virtual screen of KEGG_DRUG molecular database of accepted drugs classificated by ATC codes
Context
ATC is a classification of accepted medicines by therapeutic effects (https://en.wikipedia.org/wiki/Anatomical_Therapeutic_Chemical_Classification_System , https://www.whocc.no/atc_ddd_index/)

For examples of usage of this notebook and Virtual Screening results for COVID19 search of analogues, look at Jupyter Notebook in Github Repository (https://github.com/tonibois/KeggDrugVirtualScreening)

Content
This file contains information about ATC classificated drugs which structure is available via LigandBox repository (http://www.mypresto5.com/ligandbox/cgi-bin/lbox_download.cgi?LANG=en) and information about the relation of that compounds using KEGG identifiers and its classification in ATC (https://www.genome.jp/kegg/drug/).

The data has been processed using OpenBabel (http://openbabel.org/wiki/Main_Page) and Python programs to generate molecular identifiers and properties like fingerprints (FP2 and FP4, http://openbabel.org/wiki/Tutorial:Fingerprints), 3D-Electrosthatic molecular descriptors like PED (https://www.nature.com/articles/srep43738) and molecular properties (as Lipinski Rule of 5 parameters).

Note that this repository does not contain all ATC classificated drugs because many of them hasn't the structure available in the KEGG_DRUG repository and drug accepted list is changing over time, some molecules are deleted and other added.

Acknowledgements
Thanks to fundations as KEGG DRUG, LigandBox and OpenBabel to make public databases and free software tools.

Inspiration
This dataset is a compact collection of information about the compounds that are classified in ATC database as accepted drugs for medical care. Despite is not very extend dataset, the information that is inside has been validated by the long procedure that makes a pharmaceutical compound to the first experiments to the market. The big value of this database is that the compounds are classified by therapeutic action and then it can be easy to make groups of classes as for instance Antivirals or a subgroup ot them as Protease inhibitors (J05AE class).

This could help to find similar analogues by pharmaceutical industry to more efficient search for candidates in drug datasets or the generate therapeutic class kernels in order to improve succesful discovery rate chances.
