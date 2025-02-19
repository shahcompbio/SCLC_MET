
# FOXA2 Promotes Metastatic Competence in Small Cell Lung Cancer

**Authors:**  
Kenta Kawasaki<sup>1,2</sup>, Sohrab Salehi<sup>3</sup>, Yingqian A. Zhan<sup>4</sup>, Kevin Chen<sup>5</sup>, Jun Ho Lee<sup>1</sup>, Eralda Salataj<sup>4</sup>, Hong Zhong<sup>2</sup>, Parvathy Manoj<sup>2</sup>, Dennis Kinyua<sup>2</sup>, Barbara P. Mello<sup>2</sup>, Harsha Sridhar<sup>2</sup>, Sam E. Tischfield<sup>6</sup>, Irina Linkov<sup>7</sup>, Nicholas Ceglia<sup>3</sup>, Matthew Zatzman<sup>3</sup>, Eli Havasov<sup>3</sup>, Neil Shah<sup>2</sup>, Fanli Meng<sup>6</sup>, Brian Loomis<sup>6</sup>, Umesh K. Bhanot<sup>7</sup>, Esther Redin<sup>2</sup>, Elisa de Stanchina<sup>5</sup>, Pierre-Jacques Hamard<sup>4</sup>, Richard P. Koche<sup>4</sup>, Andrew McPherson<sup>3</sup>, Álvaro Quintanal-Villalonga<sup>2</sup>, Sohrab P. Shah<sup>3</sup>, Joan Massagué<sup>1,8</sup>,* Charles M. Rudin<sup>2,8</sup>,*

---
`
<sup>*</sup>Corresponding authors  
<sup>1</sup> [Institution 1]  
<sup>2</sup> [Institution 2]  
<sup>3</sup> [Institution 3]  
<sup>4</sup> [Institution 4]
<sup>5</sup> [Institution 5]  
<sup>6</sup> [Institution 6]  
<sup>7</sup> [Institution 7]  
<sup>8</sup> [Institution 8]`



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
- co_expression.ipynb
- fetal_scores.ipynb
- tmb_comparison.ipynb
```

Here we describe each file. 

`common_utils.py`: Utility functions that are shared across the notebooks. 

`preprocessing.ipynb` Functions for dublet detection.

`FOXA2_compartments.ipynb`: Generates figures related to tissue compartment expressing FOXA2.

`integration_harmony.ipynb`: Generates figures related to batch correction using harmony.

`integration_scvi.ipynb`: Generates figures related to batch correction using scVI.

`co_expression.ipynb`: Generates figures related to co-expression of FOXA2 and select TFs.

`fetal_scores.ipynb`: Generates figures related to fetal scores in FOXA2+ and FOXA2- cells.

`tmb_comparison.ipynb`: Generates figures related to the study of the tumor mutation burden.
