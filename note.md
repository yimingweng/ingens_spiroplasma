# 2021.08.01
1. According to [Weng et al.,](https://pubmed.ncbi.nlm.nih.gov/33938156/) the endosymbiotic *Spiroplasma* is present in the *Nebria ingens* complex, the alpine ground beetle species endemic to the Sierra Nevada of California.
2. For each of the 384 specimens of *N. ingens* complex, we have low-coverage whole genome sequence data available. 
3. The project is to detect and analyze the genomic data of the *Spiroplasma* from the *Nebria ingens* complex.
4. The work was done in [Schoville's Lab](https://molecularecology.russell.wisc.edu/) and the participants are Sean Schoville, Yi-Ming Weng, and Robert Hall.
5. The first step is to scan through the genomes of the beetle samples.
6. The working directory is in the Fuji computer under the path: `/media/Jade/YMW/spiroplasma/spiroscan`. In this repository, it is the scanning folder.
7. The result is summarized in the [distribution map](https://github.com/yimingweng/ingens_spiroplasma/blob/main/distribution_map/infection_rate.pdf)

<br />

# 2021.08.02
**Calling variants from the *Spiroplasma* bam files**
1. After extracting the *Spiroplasma* sequences from the scanning process, we are going to use this sequence data to call the variance.
2. This is for reconstructing the population structure of the *Spiroplasma*.
3. The scripts and the result are in the [varients_calling](https://github.com/yimingweng/ingens_spiroplasma/tree/main/varients_calling) folder. 

<br />

# 2021.08.04
**Population structure of the *Spiroplasma***
1. We use PCA and sNMF to visualize the gene structure (population structure) of the *Spiroplasma*.
3. The scripts and the result are in the folder sprio_structure. 