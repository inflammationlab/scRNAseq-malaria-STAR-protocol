## scRNAseq malaria STAR protocol - data analysis

<br>

---

<br>

### Table of contents

   + [Overview](#overview)

   + [Raw Data](#raw-data)

   + [Data Analysis](#data-analysis)
   
   + [Authors](#authors)

<br>

---

<br>

### Overview

Documentation supporting the single-cell data analyses described in the following **STAR protocol** - _Single Cell RNA Sequencing and Analysis of Rodent Blood Stage Plasmodium_ - which refer to the published work in _Cell Metabolism_: [A hypometabolic defense strategy against malaria](https://doi.org/10.1016/j.cmet.2022.06.011). 

<br>

---

<br>

### Raw Data

The raw sequencing data has been submitted to ArrayExpress under the accession number [E-MTAB-10939](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-10939?query=E-MTAB-10939).

<br>

---

<br>

### Data Analysis

The G6pc1fl/lfl (or `gt1` in the STAR protocol) and G6pc1AlbΔ/Δ (`gt2`) samples described in the paper correspond to `Lox2` and `Cre3` in the code. 

Data analysis notebooks and scripts describing the following sections: 

   1. _QC, filtering and exploratory clustering, dimensional reduction and DGE_
   
      + _R markdown notebook_: 
      
         + `gt1`: [02_2012_miguel_elisa_lox2_processing.Rmd](https://github.com/inflammationlab/scRNAseq-malaria-STAR-protocol/blob/main/report/02_2012_miguel_elisa_lox2_processing.Rmd) 

         + `gt2`: [01_2012_miguel_elisa_cre3_processing.Rmd](https://github.com/inflammationlab/scRNAseq-malaria-STAR-protocol/blob/main/report/01_2012_miguel_elisa_cre3_processing.Rmd) 
      
      + _html report_: 
      
         + `gt1`: [02_2012_miguel_elisa_lox2_processing.html](https://inflammationlab.github.io/scRNAseq-malaria-STAR-protocol/reports/02_2012_miguel_elisa_lox2_processing.html) 

         + `gt2` (_too big to render - please download the file from GitHub_): [01_2012_miguel_elisa_cre3_processing.html](https://github.com/inflammationlab/scRNAseq-malaria-STAR-protocol/blob/main/docs/reports/01_2012_miguel_elisa_cre3_processing.html)

      
   2. _Integration, clustering, dimensional reduction, DGE and functional enrichment_
   
      + _R markdown notebook_: [07_2012_miguel_elisa_cre3_lox2_integration_w_new_plots.Rmd](https://github.com/inflammationlab/scRNAseq-malaria-STAR-protocol/blob/main/report/07_2012_miguel_elisa_cre3_lox2_integration_w_new_plots.Rmd) 
      
      + _html report_: [07_2012_miguel_elisa_cre3_lox2_integration_w_new_plots.html](https://inflammationlab.github.io/scRNAseq-malaria-STAR-protocol/reports/07_2012_miguel_elisa_cre3_lox2_integration_w_new_plots.html)


   3. _Velocity analysis_
   
      3.1. 
      
         + _R markdown notebook_: [09_2012_miguel_elisa_RNA_velocity_part1.Rmd](https://github.com/inflammationlab/scRNAseq-malaria-STAR-protocol/blob/main/report/09_2012_miguel_elisa_RNA_velocity_part1.Rmd) 

         + _html report_: [09_2012_miguel_elisa_RNA_velocity_part1.html](https://inflammationlab.github.io/scRNAseq-malaria-STAR-protocol/reports/09_2012_miguel_elisa_RNA_velocity_part1.html)

      3.2. 
      
         + _bash script_: [velocyto_script.sh](https://github.com/inflammationlab/scRNAseq-malaria-STAR-protocol/blob/main/scripts/velocyto_script.sh)

         + _command_: `./velocyto_script.sh &> velocyto_log.log`

   
      3.3. 
      
         + _ipython notebook_: [09_2012_miguel_elisa_RNA_velocity_part2.ipynb](https://github.com/inflammationlab/scRNAseq-malaria-STAR-protocol/blob/main/report/09_2012_miguel_elisa_RNA_velocity_part2.ipynb)

         + _html report_: [09_2012_miguel_elisa_RNA_velocity_part2.html](https://inflammationlab.github.io/scRNAseq-malaria-STAR-protocol/reports/09_2012_miguel_elisa_RNA_velocity_part2.html)

   4. Plotting markers and enriched functions (_not included in the STAR protocol_): 

      4.1. **Markers**: 

         + _R markdown notebook_: [10_2012_miguel_elisa_plot_markers.Rmd](https://github.com/inflammationlab/scRNAseq-malaria-STAR-protocol/blob/main/report/10_2012_miguel_elisa_plot_markers.Rmd) 

         + _html report_: [10_2012_miguel_elisa_plot_markers.html](https://inflammationlab.github.io/scRNAseq-malaria-STAR-protocol/reports/10_2012_miguel_elisa_plot_markers.html)

      4.2. **Enriched functions**:

         + _R markdown notebook_: [11_2012_miguel_elisa_func_enrich_plot_DEG_markers.Rmd](https://github.com/inflammationlab/scRNAseq-malaria-STAR-protocol/blob/main/report/11_2012_miguel_elisa_func_enrich_plot_DEG_markers.Rmd) 

         + _html report_: [11_2012_miguel_elisa_func_enrich_plot_DEG_markers.html](https://inflammationlab.github.io/scRNAseq-malaria-STAR-protocol/reports/11_2012_miguel_elisa_func_enrich_plot_DEG_markers.html)
   
<br>

---

<br>

### Authors

**Elisa Jentho1,2#`*`, António G. G. Sousa1#`*`, Susana Ramos1, Temitope W. Ademolue1, João Sobral1, João Costa1, Denise Brito1, Marta Manteiro1, Ricardo B. Leite1, Jingtao Lilue and Miguel P. Soares1,6##**

1 Instituto Gulbenkian de Ciência, Oeiras, Portugal. 

2 Jena University Hospital, Department of Anesthesiology and Intensive Care Medicine, Friedrich-Schiller-University, Jena, Germany. 

6 Lead contact

`*` These authors contributed equally

`#` Technical Correspondence: <ejentho@igc.gulbenkian.pt> and <aggode@utu.fi>

`##` Correspondence: <mpsoares@igc.gulbenkian.pt> 
   
<br>

---

<br>

<br>

<br>

<br>

<br>

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg

