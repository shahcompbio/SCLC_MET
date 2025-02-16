# SCLC_MET

Code base related to the following paper:

```
FOXA2 promotes metastatic competence in small cell lung cancer
Kenta Kawasaki1,2, Sohrab Salehi3, Yingqian A. Zhan4, Kevin Chen5, Jun Ho Lee1, Eralda Salataj4, Hong Zhong2, Parvathy Manoj2, Dennis Kinyua2, Barbara P. Mello2, Harsha Sridhar2, Sam E. Tischfield6, Irina Linkov7, Nicholas Ceglia3, Matthew Zatzman3, Eli Havasov3, Neil Shah2, Fanli Meng6, Brian Loomis6, Umesh K. Bhanot7, Esther Redin2, Elisa de Stanchina5, Pierre-Jacques Hamard4, Richard P. Koche4, Andrew McPherson3, Álvaro Quintanal-Villalonga2, Sohrab P. Shah3, Joan Massagué1,8,*, Charles M. Rudin2,8,*
```


## Organization

The repository is organized as follows.

`src/Notebooks` contains the notebooks to reproduce figures in the study.
The directory contains the following files.

```
- common_utils.py
- preprocessing.ipynb
- FOXA2_compartments.ipynb
- integration_harmony.ipynb
- integration_scvi.ipynb
- tmb_comparison.ipynb
```

Here we describe each file. 

`common_utils.py`: Utility functions that are shared across the notebooks. 

`preprocessing.ipynb` Functions for dublet detection.

`FOXA2_compartments.ipynb`: Generates figures related to tissue compartment expressing FOXA2.

`integration_harmony.ipynb`: Generates figures related to batch correction using harmony.

`integration_scvi.ipynb`: Generates figures related to batch correction using scVI.

`tmb_comparison.ipynb`: Generates figures related to the study of the tumor mutation burden.