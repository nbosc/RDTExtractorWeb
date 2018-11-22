# DjangoAPI

## Download

- Download zip file 

or 

- `git clone https://github.com/phi-grib/DjangoAPI.git`

## Install

Install the enviroment

`conda env create -f environment.yml`

## Run
In the root folder execute:

`python manage.py runserver`

## Introduction
This tool is designed to extract data from the _in vivo_ repeat-dose toxicity (RDT) studies' database generated within the context of the [eTOX](http://www.etoxproject.eu/) project. These data are expanded using an histopathological observation and an anatomical entity ontologies. The [histopathological ontology](https://github.com/Novartis/hpath/blob/master/LICENSE.txt) is obtained from Novartis and can be used under the Apache License 2.0. The anatomical entities ontology is extracted from the following paper:
- [Hayamizu TF, Mangan M, Corradi JP, Kadin JA, Ringwald M. Genome Biol. 2005; 6(3): R29](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1088948/)

In order to be able to aggregate the data by parent compound, some pre-processing has to be done to data as they exist in the database. Each substance is standardised according to the following protocol:
- From [this repository](https://github.com/bet-gregori/standardiser) use the _process_smiles.std_ method to standardize, discard mixtures, discard compound with metal ions, and remove all salts. Also use the _neutralise.run_ method to neutralise all charges when possible.
- Using [molVS](https://molvs.readthedocs.io/en/latest/guide/intro.html), get the canonical tautomer.

This project is an extension of the work published in the following paper:
- [López-Massaguer O, Pinto-Gil K, Sanz F, Amberg A, Anger LT, Stolte M, Ravagli C, Marc P, Pastor M. Toxicol Sci. 2018 Mar; 162(1): 287–300.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5837688/)

## Manual
Exract studies' findings based on the given filtering and the organs' and
morphological changes' ontologies-based expansions of these findings.

### Output example
On clicking the 'Extract' button, two output files are generated, one with quantitative and the other with qualitative data. Both have a caption summarising the filtering criteria applied. After this caption, they both have a table with the data aggergated by parent compound. The table contains several fixed columns, namely 6 at the begining:
- inchi_key: Parent compound's InChIKey.
- study_count: Number of relevant studies (according to the current filtering scheme) in which the compound appears.
- dose_min: Minimum dose at which the compound has been tested among the relevant studies.
- dose_max: Maximum dose at which the compound has been tested among the relevant studies.
- min_observation_dose: Minimum dose for which a relevant finding (according to the current filtering scheme) has been reported for the compound.
- is_active: Boolean indicating whether the substance has been found to have any toxicity according to the current finding-related filtering criteria.
And two at the end:
- subst_id: All substance IDs corresponding to the parent compound.
- std_smiles: Smiles string corresponding to the standardised parent compound.

Between these two groups, there is a column for each relevant finding. In these columns a value is provided if the finding is reported for the given substance, and it is empty otherwise. The value will be the number of studies that report the finding in the qualitative file, and the minimum dose at which the finding is reported in the quantitative file.

This is an example of the qualitative output: 
![qualiative](https://raw.githubusercontent.com/phi-grib/DjangoAPI/master/img/qual.JPG)

This is an example of the quantitative output: 
![quantitative](https://raw.githubusercontent.com/phi-grib/DjangoAPI/master/img/quant.JPG)
