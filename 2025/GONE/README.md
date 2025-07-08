# GONE Practical

**Author:** Marianne Dehasque  
**Date:** June 30, 2025

## Introduction

### GONE
For this practical, we will use `GONE`, a software that infers recent effective population size from genetic data based on linkage disequilibrium. 

The original `GONE` paper can be found [here](https://academic.oup.com/mbe/article/37/12/3642/5869049). However, today we will use `GONE2`, which hasn't been officially published yet, but is faster and more user-friendly. Whereas the original GONE software depends on multiple scripts and a separate INPUT_PARAMETERS_FILE, GONE2 can be run directly from the command line with different parameters. 

The `GONE2` GitHub repository with installation and usage instructions can be found [here](https://github.com/esrud/GONE2). Note that for the practical today all necessary software and tools have already been installed.

### Dataset
We will use the killer whale dataset from [Kardos et al, 2023 (_Nature Ecology and Evolution_)](https://www.nature.com/articles/s41559-023-01995-0#Sec21). In this paper, Kardos and colleagues study the population dynamics of an endangered North Pacific killer whale population. The main author kindly provided a pre-filtered dataset consisting of genotype information from 151 individuals from five different killer whale populations that we can use for today's practical.

## Data types
In this practical, we will use two types of data: PLINK PED/MAP files, and PLINK binary BED files.

Given their compact size, PLINK binary files are most commonly used to store and filter data. Since `GONE2` requires PED/MAP files as input, a workflow will typically convert between these two data types. This can be done using the following commands:

```bash
# Convert PLINK binary files to PED/MAP file
plink --bfile <filename> --recode --out <filename>

# Convert PLINK PED/MAP to PLINK binary files
plink --file <filename> --make-bed --out <filename>
```

Below more detailed information on both data types is given. This is additional reference information. It is not needed to run the practical, so we can skip it now, but it can be useful if we want to poke into the data files.

---

### PLINK PED/MAP Files

PLINK PED/MAP files are used to store genotype data in a text-based format. They consist of two main files:

* `.ped` file: Contains genotype and sample information.
* `.map` file: Contains variant information.

#### .ped file format
The `.ped` file is a tab-delimited text file where each row corresponds to an individual sample. It contains both sample metadata and genotype data. The columns are structured as follows:

* `Family ID`: Identifier for the family (can be set to 0 if not applicable).
* `Individual ID`: Unique identifier for the individual.
* `Paternal ID`: Identifier for the father (0 if not available).
* `Maternal ID`: Identifier for the mother (0 if not available).
* `Sex`: Sex of the individual (1 = male, 2 = female, 0 = unknown).
* `Phenotype`: Phenotype information (1 = unaffected, 2 = affected, -9 = missing).
* `Genotype Data`: Genotypes for each variant, represented as pairs of alleles (e.g., A A, A B, B B). Missing genotypes are denoted as 0 0.

#### .map file format
The `.map` file is a tab-delimited text file where each row corresponds to a variant. It contains the following columns:

* `Chromosome`: Chromosome number or ID.
* `Variant ID`: Unique identifier for the variant (e.g., rsID).
* `Genetic Distance`: Genetic distance (can be set to 0 if not available).
* `Base-pair Position`: Physical position of the variant on the chromosome.

### PLINK Binary Files

PLINK binary files are used to store genotype data in a compact and efficient binary format. They typically consist of three main files:

* `.bed` file: Contains the binary genotype data.
* `.bim` file: Contains variant information.
* `.fam` file: Contains individual sample information.

#### .bed file format
The `.bed` file stores the genotype data in a compact binary format. It does not have a human-readable structure but is designed for efficient storage and access by PLINK and other compatible tools. 

**Genotype Encoding:**
* 00: Homozygous for the reference allele (AA)
* 01: Missing genotype
* 10: Heterozygous (AB)
* 11: Homozygous for the alternate allele (BB)

#### .bim file format
The `.bim` file is a text file that contains variant information. Each row corresponds to a variant, and it has the following columns:

* `chrom`: Chromosome number or ID
* `variant ID`: Unique identifier for the variant (e.g., rsID)
* `genetic distance`: Genetic distance (can be set to 0 if not available)
* `base-pair position`: Physical position of the variant on the chromosome
* `allele 1`: Reference allele (usually coded as the minor allele)
* `allele 2`: Alternate allele

#### .fam file format
The `.fam` file is a text file that contains information about each individual sample in the dataset. Each row corresponds to an individual and has the following columns:

* `Family ID`: Identifier for the family (can be set to 0 if not applicable)
* `Individual ID`: Unique identifier for the individual
* `Paternal ID`: Identifier for the father (0 if not available)
* `Maternal ID`: Identifier for the mother (0 if not available)
* `Sex`: Sex of the individual (1 = male, 2 = female, 0 = unknown)
* `Phenotype`: Phenotype information (1 = unaffected, 2 = affected, -9 = missing)

---

## Exploring the data

Open this link in your browser to set up the GONE training environment for this practical: [![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://classroom.github.com/a/ZG65xSNV)

We should be presented with a page where we can create a new GitHub Codespace. For this practical, the default environment with 2 cores should suffice.

In the training environment, we should see the following folders:

* **DATA**: contains the data for this practical
* **SCRIPTS**: contains (plotting) scripts you can use
* **GONE2**: contains the GONE2 installation

For this practical, we will start from PLINK binary files containing genotype information of 151 killer whales: `KW_all.bed`, `KW_all.bim` and `KW_all.fam`. The files are stored in the `./DATA` folder. The following filters have already been applied:

* Biallelic loci only
* No indels
* No loci with unusually high heterozygosity or read depth
* Minimum Allele Frequency (MAF) of 0.05
* Autosomal data only

> **Question:** Why do we filter our data? What's the relevance of each of these filters?

Additionally, the dataset has been reduced to chr1, chr10 and chr20. This is simply to speed up the analyses for today's practical.

Let's start by familiarizing ourselves with the dataset first. Killer whales in the Eastern North Pacific can be subdivided into genetically differentiated populations that differ in diet, behaviour, and distribution. In our dataset, individuals have been assigned to one of the following populations: Transient (TKW), offshore, Alaska residents (ARKW), Southern residents (SRKW), and Northern residents (NRKW). The geographic distribution of these populations is given in Figure 1 below.

**Figure 1:** [Geographic distribution map of killer whale populations to be inserted here]

To test if (1) we find this population structure and if (2) there are any outliers or mislabeled samples, we will first conduct a Principal Component Analysis (PCA). PCA is a statistical technique that reduces the dimensionality of our data by identifying the main sources of variation. In the field of genomics, it is typically used as a first step to explore the dataset.

### Setting up the workspace

First, create a new directory to store all results of today's practical. Assuming we are still in the login directory, we can type:

```bash
mkdir RESULTS/
cd RESULTS/
```

### Running PCA analysis

We will perform a PCA with the PLINK software. Before using a (new) software, it is always a good idea to read the software instructions first. For most software programs this can be simply done by using the `--help` flag.

To view PLINK's usage instructions, we can use:

```bash
plink --help
```

There is no need to completely read the manual now, but it is an important source of information if we want to better understand parameters and/or syntax. 

To run a PCA, we can use the following command:

```bash
plink --bfile ../DATA/KW_all --pca 151 --out KW_all
```

This should generate two files with `.eigenval` and `.eigenvec` file extensions. To plot the results, we can use the Python script provided under `./SCRIPTS`. You may have noticed that the PLINK output refers to "people". This is because PLINK was originally developed for human genetic data.

```bash
python ../SCRIPTS/plot_PCA.py 
```

Inspect the `pca_plot.png` file. 

> **Question:** Are there five genetically differentiated populations? Are there any outliers?

## Preparing data for GONE analysis

Since some of the populations in our dataset only have a few individuals, we will run the GONE analysis for the Southern residents (SRKW) and Alaska residents (ARKW) populations. Generally, a minimum sample size of 10-15 individuals is recommended for small populations, and more for larger populations to accurately infer population size. GONE generally doesn't require any hard filtering (i.e. the filters applied in the pre-filtered dataset are fine). However, the program does not allow any missing data. We will first run the analysis for the SRKW population.

### Subsetting by population

Let's first subdivide our main dataset into different populations. To achieve this, we can use the code below. 

```bash
plink --bfile ../DATA/KW_all --keep ../DATA/SRKW.txt --make-bed --out SRKW
```

The `--keep` flag accepts a space/tab-delimited text file with family IDs in the first column and within-family IDs in the second column, and removes all unlisted samples from the current analysis. To find the family IDs and within-family IDs, we can always check the `KW_all.fam` file.

### Checking for missing data

Next, we will inspect our population datafiles for individuals with excessive missing data.

```bash
plink --bfile SRKW  --missing --out SRKW
```

This will generate a `.imiss` file.

> **Question:** Inspect the `.imiss` file. How many SNPs were genotyped in the whole population? Which individual has the most missing data?

There are many reasons why missingness might vary between individuals, including depth of sequencing. In some cases, it could also be an indication of technical problems. Depending on the study, individuals with too much missing data can be removed for certain analyses. Our dataset looks fine, so we won't be removing any individuals from the dataset here.

### Converting to PED/MAP format

Finally, we will remove sites with missing data and transform our data into the `.map` `.ped` file format. We can use the following command:

```bash
plink --bfile SRKW --geno 0 --recode --out SRKW.geno0
```

Make sure to inspect the PLINK documentation if we are unsure about any of the flags in our command above.

## Running the GONE analysis

First, let's check the available options for GONE2:

```bash
gone2 -help
```

We will use parameters based on the original paper to run the GONE analysis. Additionally, to speed up the analysis, we will run the analysis on 10,000 subsampled SNPs:

```bash
# This will take ~1 minute
gone2 -g 0 -r 1 -u 0.02 -s 10000 SRKW.geno0.ped -o SRKW_geno0.u0.02.s10k
```

### Parameter explanation:
* `-g 0`: Type of genotyping data. 0:unphased diploids.
* `-r 1`: Constant recombination rate of 1 cM/Mb - This is a typical rate for large mammals
* `-u 0.02`: Upper bound of recombination rates to be considered (2 cM, default is 5 cM) - This parameter was set to 0.02 in the original paper to mitigate the possibility of bias arising from population substructure in recent population history.
* `-s 10000`: Subsample to 10,000 SNPs for faster analysis
* `-o`: Output file prefix

### Visualizing results

Once the analysis is complete, we can visualize the results:

```bash
python ../SCRIPTS/plot_GONE.py SRKW_geno0.u0.02.s10k_GONE2_Ne
```

## Exercises

### Exercise 1
Repeat the GONE2 analysis with the Alaska residents (ARKW) population. Start with generating the appropriate files as outlined under "Preparing data for GONE analysis". Compare the results of the Alaska residents to the Southern residents. 

### Exercise 2
If there is time left, we can start playing around with the different parameters in GONE. What happens if we change the recombination rate (`-r`) or upper bound of recombination rates (`-u`)?

### Exercise 3
If there is even more time left, we can also test the effect of population specification on demographic inference. What happens if we misspecify our population and run the analysis on SRKW and ARKW together? What if we reduce the number of individuals in the SRKW or ARKW populations?

## Additional Resources

- [GONE2 GitHub repository](https://github.com/esrud/GONE2)
- [Original GONE paper](https://academic.oup.com/mbe/article/37/12/3642/5869049)
- [PLINK documentation](https://www.cog-genomics.org/plink/)
- [Killer whale dataset reference](https://www.nature.com/articles/s41559-023-01995-0#Sec21)

## Repository Structure

```
.
├── DATA/
│   ├── KW_all.bed
│   ├── KW_all.bim
│   ├── KW_all.fam
│   ├── SRKW.txt
│   └── ARKW.txt
├── SCRIPTS/
│   ├── plot_PCA.py
│   └── plot_GONE.py
├── GONE2/
│   └── [GONE2 installation files]
└── RESULTS/
    └── [Generated during practical]
```
