Here we only provide two datasets with the limited space ( [Mammary epithelial cells (bct_raw.zip)](./bct_raw.zip) and  [macaque_retina (macaque_raw.zip)](./macaque_raw.zip)). And the full dataset will be available from [10.6084/m9.figshare.20499630](10.6084/m9.figshare.20499630)). Below is the description of all datasets: 


__Supplemental Table1__: Datasets analyzed in scDML 

| Species                | Tissue                   | Data Source                                                  | Batches                        | Dataset Dimensions       |
| ---------------------- | ------------------------ | ------------------------------------------------------------ | ------------------------------ | ----------------------- |
| Simluation1            | -                        | https://figshare.com/articles/dataset/Benchmarking_atlas-level_data_integration_in_single-cell_genomics_-_integration_task_datasets_Immune_and_pancreas_/12420968 (sim2_norm.h5ad)(Luecken and Theis, 2019) | 4 batches                      | 19138 cells <br> 10000 genes  |
| Simulation2            | -                        | https://figshare.com/articles/dataset/Benchmarking_atlas-level_data_integration_in_single-cell_genomics_-_integration_task_datasets_Immune_and_pancreas_/12420968 (sim1_1_norm.h5ad) (Luecken and Theis, 2019) | 6 batches                      | 12097 cells <br> 9979 genes   |
| Mouse                  | Mammary epithelial cells | https://github.com/NBISweden/single-cell_sib_scilifelab/tree/master/datasets/SCE_MammaryEpithelial_x3.rds; (PumbedID_30089273, PubmedID_29158510, PubmedID_29225342) | 3 batches                      | 9288 cells <br> 1222 genes    |
| Human                  | Pancreas                 | https://satijalab.org/seurat/archive/v3.2/integration.html.(standard workflow) <br> SeuratData::InstallData(“panc8”) | 8 batches                      | 14890 cells <br> 34363 genes  |
| Macaque                | Retina (Bipolar cells)    | GSE11848(Peng et al., 2019)                                  | 2 regions <br> 4 animals <br> 30 samples | 30302 cells <br> 36162 genes  |
| Mouse                  | Retina (Bipolar cells)    | GSE81905(Shekhar et al., 2016)                               | 6 batches | 23494 cells <br> 13166 genes  |
| Human and Mouse        | Lung                     | GSE133747(Raredon et al., 2019)                              | 2 batches  | 20760 cells <br> 62781 genes  |
| Mouse                  | Brain                    | GSE116470 <br> GSE110823 <br> http://scanorama.csail.mit.edu/data.tar.gz. | 2 batches  | 833206 cells <br> 17745 genes |
| Human           | Heart              | [GSE183852](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183852) (Koenig et al., 2022) | 45 batches              | 269794 cells <br>45068 genes |
| Human     | Heart            | [www.heartcellatlas.org](http://www.heartcellatlas.org/) (Litviňuková, et al. 2020) | 147 batches           | 486134 cells <br>33538 genes |
| Adam (Mouse) | Kidney        | https://drive.google.com/drive/folders/1BIZxZNbouPtGf_cyu7vM44G5EcbxECeu (Adam)(Adam et al., 2017) | 1 batch                                  | 3660 cells <br> 23797 genes |
| Mural (Human) | Pancreas | https://drive.google.com/drive/folders/1BIZxZNbouPtGf_cyu7vM44G5EcbxECeu (Muraro)(Muraro et al., 2016) | 1 batch | 2122 cells <br> 19046 genes |
| QX_Limb_Muscle (Mouse) | Limb Muscle | https://drive.google.com/drive/folders/1BIZxZNbouPtGf_cyu7vM44G5EcbxECeu (Quake_10x_Limb_Muscle)(The Tabula Muris Consortium et al., 2018) | 1 batch | 3909 cells <br> 23341 genes |


# Reference

Adam,M. **et al.** (2017) Psychrophilic proteases dramatically reduce single-cell RNA-seq artifacts: a molecular atlas of kidney development. **Dev. Camb. Engl.**, **144**, 3625–3632.

Luecken,M.D. and Theis,F.J. (2019) Current best practices in single‐cell RNA‐seq analysis: a tutorial. **Mol. Syst. Biol.**, **15**.

Muraro,M.J. **et al.** (2016) A Single-Cell Transcriptome Atlas of the Human Pancreas. **Cell Syst.**, **3**, 385-394.e3.

Peng,Y.-R. **et al.** (2019) Molecular Classification and Comparative Taxonomics of Foveal and Peripheral Cells in Primate Retina. **Cell**, **176**, 1222-1237.e22.

Raredon,M.S.B. **et al.** (2019) Single-cell connectomic analysis of adult mammalian lungs. **Sci. Adv.**, **5**, eaaw3851.

Shekhar,K. **et al.** (2016) Comprehensive Classification of Retinal Bipolar Neurons by Single-Cell Transcriptomics. **Cell**, **166**, 1308-1323.e30.

The Tabula Muris Consortium **et al.** (2018) Single-cell transcriptomics of 20 mouse organs creates a Tabula Muris. **Nature**, **562**, 367–372.

Koenig A.L., et al (2022). Single-cell transcriptomics reveals cell-type-specific diversification in human heart failure. **Nat Cardiovasc Res**. **1(3)**:263-280. doi: 10.1038/s44161-022-00028-6. 

Litviňuková M,et al. (2020). Cells of the adult human heart. **Nature** . **588**(7838):466-472.

