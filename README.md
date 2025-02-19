
# FOXA2 Promotes Metastatic Competence in Small Cell Lung Cancer

**Authors:**  
Kenta Kawasaki<sup>1,2</sup>, Sohrab Salehi<sup>3</sup>, Yingqian A. Zhan<sup>4</sup>, Kevin Chen<sup>5</sup>, Jun Ho Lee<sup>1</sup>, Eralda Salataj<sup>4</sup>, Hong Zhong<sup>2</sup>, Parvathy Manoj<sup>2</sup>, Dennis Kinyua<sup>2</sup>, Barbara P. Mello<sup>2</sup>, Harsha Sridhar<sup>2</sup>, Sam E. Tischfield<sup>6</sup>, Irina Linkov<sup>7</sup>, Nicholas Ceglia<sup>3</sup>, Matthew Zatzman<sup>3</sup>, Eli Havasov<sup>3</sup>, Neil Shah<sup>2</sup>, Fanli Meng<sup>6</sup>, Brian Loomis<sup>6</sup>, Umesh K. Bhanot<sup>7</sup>, Esther Redin<sup>2</sup>, Elisa de Stanchina<sup>5</sup>, Pierre-Jacques Hamard<sup>4</sup>, Richard P. Koche<sup>4</sup>, Andrew McPherson<sup>3</sup>, Álvaro Quintanal-Villalonga<sup>2</sup>, Sohrab P. Shah<sup>3</sup>, Joan Massagué<sup>1,8</sup>,* Charles M. Rudin<sup>2,8</sup>,*

---

<sup>*</sup>Corresponding authors  
<sup>1</sup> [Cancer Biology and Genetics Program, Sloan Kettering Institute, Memorial Sloan Kettering Cancer Center, New York, NY 10065, USA]  
<sup>2</sup> [Department of Medicine, Memorial Sloan Kettering Cancer Center, New York, NY, USA]  
<sup>3</sup> [Computational Oncology Program, Memorial Sloan Kettering Cancer Center, New York NY, USA]  
<sup>4</sup> [Center for Epigenetics Research, Memorial Sloan Kettering Cancer Center, New York, NY, USA]
<sup>5</sup> [Antitumor Assessment Core, Memorial Sloan Kettering Cancer Center, New York NY, USA]  
<sup>6</sup> [Marie-Josée and Henry R. Kravis Center for Molecular Oncology, Memorial Sloan Kettering Cancer Center, NY, USA]  
<sup>7</sup> [Pathology Core Facility, Memorial Sloan Kettering Cancer Center, New York, NY, USA.]  
<sup>8</sup> [Weill Cornell Medicine Graduate School of Medical Sciences, New York, NY, USA]




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
