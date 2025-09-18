# Code and Data for "Ribosome dynamics govern mRNA stability and translation across natural and synthetic systems"

This repository contains the code and data access instructions to reproduce the figures and analyses presented in our manuscript, "[Insert Paper Title Here]".

## Data Download

The raw sequencing data is publicly available from the China National Center for Bioinformation (CNCB) under the following accessions.

### Accession: HRA011266

This dataset includes TRiP-seq data for various conditions and cell lines.

* **Link:** [`https://ngdc.cncb.ac.cn/gsa-human/browse/HRA011266`](https://ngdc.cncb.ac.cn/gsa-human/browse/HRA011266)
* **Contents:**
    * **HEK293T cells:**
        * 4sU labeling at 30 min, 1 h, and 2 h
        * Puromycin treatment for 2 h
        * STM2457 (20 µM) treatment for 2 h
        * YTHDF1, YTHDF2, and YTHDF3 Knockout (KO) lines
    * **THP-1 cells (unstimulated):**
        * 4sU labeling at 30 min, 1 h, and 2 h
    * **THP-1 cells (LPS stimulated):**
        * 4sU labeling at 30 min, 1 h, and 2 h

### Accession: HRA013215

This dataset includes the Pulse-Chase experiment data for HEK293T cells.

* **Link:** [`https://ngdc.cncb.ac.cn/gsa-human/browse/HRA013215`](https://ngdc.cncb.ac.cn/gsa-human/browse/HRA013215)
* **Contents:**
    * **HEK293T Pulse-Chase:**
        * Pulse: 2 h
        * Chase: 1 h, 2 h, and 8 h

---

## Data Preprocessing

If you choose to start with the raw data from CNCB, follow the steps below to process it.

1.  **Run the Python preprocessing scripts** for each cell line to process the raw files:
    ```bash
    python preprocess_HEK293T.py
    python preprocess_THP1.py
    ```
2.  **Next, run the R scripts** to load the processed data into R objects for analysis:
    ```R
    "load_data_to_R_293_all_sample.R"
    "load_data_to_R_THP1_all_sample.R"
    ```
---

## Reproducing Figures

To reproduce all the figures in the paper, you can use the pre-processed data provided. This is the recommended and most straightforward approach.

1.  **Download the processed data** from Google Drive:
    > **https://drive.google.com/drive/folders/1rdiGYeeCrGBGh4Fh4qhOqvRbHUZpBayy?usp=sharing**

2.  Place the downloaded data into the project's working directory.

3.  **Run the analysis scripts** for each figure to generate the plots:
    * `Figure1.R`
    * `Figure2.R`
    * `Figure3.R`
    * `Figure4.R`
    * `Figure5.R`

*Note: Please ensure you have all the necessary R and Python packages installed. Check the script files for specific library dependencies.*
