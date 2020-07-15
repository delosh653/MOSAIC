# Welcome to MOSAIC!

<p align="center">
<img src="MOSAIC App/www/MOSAIC_Gene_Expression_Plot.png" width="500" />
</p>

Note: currently in beta testing. Since this app is currently in development, please check here for updates, especially for an updated MOSAIC citation.

This is the third step in the PAICE (Pipeline for Amplitude Integration of Circadian Exploration) Suite! This suite of tools provides high-throughput applications for circadian, ultradian, and infradian rhythms. The first step, ECHO, can be found [here](https://github.com/delosh653/ECHO), and the second step, ENCORE, can be found [here](https://github.com/delosh653/ENCORE). These applications can be run before or after running MOSAIC.

## README Outline

* Overview
* Use and First-Time Set-Up Instructions
* MOSAIC Features
* Data Format Example
* MOSAIC R Package
* Minimum Version Information
* Contact Information and Bug Reporting
* FAQ

## Overview

MOSAIC (Multi-Omics Selection with Amplitude Independent Criteria) is an R-powered application designed to find and visualize circadian and non-circadian trends in multi-omics time course data using model selection and joint modeling. To read more about this work and cite us, see [*MOSAIC: A Joint Modeling Methodology for Combined Circadian and Non-Circadian Analysis of Multi-Omics Data*](https://www.biorxiv.org/content/10.1101/2020.04.27.064147v1) by H. De los Santos, et al. (2020). To read more about the work one of the contributing models is based on, see [*ECHO: an Application for Detection and Analysis of Oscillators Identifies Metabolic Regulation on Genome-Wide Circadian Output*](https://doi.org/10.1093/bioinformatics/btz617) by H. De los Santos et al. (2019), published in *Bioinformatics*.

All images created by ECHO using data from [*Analysis of clock-regulated genes in Neurospora reveals widespread posttranscriptional control of metabolic potential*](https://www.ncbi.nlm.nih.gov/pubmed/25362047) by J. Hurley, et al. (2014) and [*Circadian Proteomic Analysis Uncovers Mechanisms of Post-Transcriptional Regulation in Metabolic Pathways*](https://doi.org/10.1016/j.cels.2018.10.014) by J. Hurley, et al. (2018).

## Use and First-Time Set-Up Instructions

Thank you for downloading MOSAIC (Multi-Omics Selection with Amplitude Independent Criteria)! MOSAIC is an app to find and visualize circadian and non-circadian trends in multi-omics time course data using model selection and joint modeling. This guide will lead you in first time set-up and use. Pictures have been provided for ease of use, using Windows 10, in the files MOSAIC README.docx and MOSAIC README.pdf, found above. A double asterisk indicates the step has an explanation below, and a tilde indicates the step is first-time set up only.

Steps: 
1.	** ~ Download [Firefox](https://www.mozilla.org/en-US/firefox/new/) or [Chrome](https://www.google.com/chrome/browser/desktop/index.html) and make it your default browser.
2.	~ [Download R](https://www.r-project.org/), if you do not already have it. 
3.	~ [Download RStudio](https://www.rstudio.com/products/rstudio/download/), if you do not already have it (RStudio Desktop is sufficient).
4.	Open RStudio.
5.	~ Copy and paste the following text into the console window (bottom left window of the RStudio Session), then press enter:

```r
install.packages("rstudioapi")
install.packages("shiny")
install.packages("ggplot2")
install.packages("VennDiagram")
install.packages("reshape2")
install.packages("minpack.lm")
install.packages("doParallel")
install.packages("foreach")
install.packages("iterators")
install.packages("doSNOW")
install.packages("colorRamps")
install.packages("ggplotify")
install.packages("gridExtra")
install.packages("dplyr")
```

This will install these packages (a set of functions that this application uses) onto your computer. This may ask for your input, so just say yes to the questions asked. If you run into errors saying “yes,” just say no instead. Note: this may take some time.

6.	Open mosaic_app.R, which should be included in the .zip file you downloaded and also contained README.pdf and README.docx. It should open in the top left window of your RStudio session.

7.	In the top right corner of the mosaic_app.R window, you should see the button, “Run App”. Click on the small downwards arrow next to it and choose “Run External”. 

8.	Now click “Run App”. This should open the ECHO application in your now default browser window (either Firefox or Chrome). The picture below is a representation in Firefox.

9.	Have fun!

** Why do I have to install either Firefox or Chrome, you ask? Why not Internet Explorer, or some other browser? Well, it is known there are problems downloading files when viewing shiny apps in Internet Explorer, so we definitely want to avoid that. However, I have not tested this app in browsers like Microsoft Edge, Safari, etc. If you can verify that these work, please let me know at delosh@rpi.edu.

## MOSAIC Features

MOSAIC's interface is divided into two sections: **Finding Trends** and **Visualizing Results**.

Within the **Finding Rhythms** tab, you can upload your RNA and protein data (.csv) and enter its information, such as time point range, resolution (in hours), and amount and type of replicates. You can then choose from a variety of preprocessing steps including smoothing, removing unexpressed genes, and normalization. You can then download your results as both a CSV (for viewing) and a .RData (for visualizations).

In the **Visualizing Results** tab, simply upload the .RData file from your results and choose from several visualization and gene subset exploration options. You can explore subsets of data under the "Gene List" tab and sort by the various output parameters, such as Period or P-Value. You can also choose from a host of automatically-generated visualizations, including summary visualizations, heat maps, gene expression plots (with or without replicates visualized), and parameter density graphs (examples displayed below).

<p align="center">
<img src="MOSAIC App/www/MOSAIC_Gene_Expression_Plot.png" width="200" /> <img src="MOSAIC App/www/MOSAIC_Summary_Viz.png" width="250" /> <img src="MOSAIC App/www/MOSAIC_Parameter_Density_Plot.png" width="200" /> <img src="MOSAIC App/www/MOSAIC_Heat_Map.png" width="250" /> <img src="MOSAIC App/www/MOSAIC_Heat_Map_Comparison.png" width="250" />
</p>

## Data Format Example

RNA and protein data should be .csvs (comma-separated values) with the first column being the expression names/labels, and the rest being numeric columns with expression data, ordered by time point, then by replicate. Missing data should be left blank. Expression names/labels should match *exactly* between corresponding RNA and protein expressions. Names/labels not found in either dataset will not be run. An example of this formatting is the following:

| RNA.Name |	TP2.1 | 	TP2.2| 	TP2.3	| TP4.1| 	TP4.2| 	TP4.3| 
| ------------- |-------------|-------------|-------------|-------------|-------------|-------------|
| Sample 1 |	1.633117905| 	| 	1.513810213| 	1.309553546 | 	1.302488129| 	|	
| Sample 2 | 	-0.630319173| 	| 	-0.510500938| 	| 	-0.543457041| 	-0.448383157|		
| Sample 3	| -0.780221402| 	| 	| 	0.178429468| 	0.306513019| 	1.376226634|

| Protein.Name |	TP2.1 | 	TP2.2| 	TP2.3	| TP4.1| 	TP4.2| 	TP4.3| 
| ------------- |-------------|-------------|-------------|-------------|-------------|-------------|
| Sample 1 |	0.465992299	| 		| 	| 		-0.037081559| 		-0.076233835 |-0.076233835 |
| Sample 2 |0.702637217	| 	0.640581226	| 	| 		1.719945272| 			 | 1.775870204|
| Sample 3	|	| 	0.885876804| 		| 		0.971783339| 		0.93201731 |0.93201731|

In this example, this is two hour resolution data taken from 2 to 4 hours, with three replicates, from 3 corresponding genes/proteins (as denoted by the same names). In RNA, the second replicate at time point 2 is entirely missing, and in protein the third replicate at time point 2 is missing. Each expression in RNA and protein has additional missing data at various time points and replicates.

Larger example datasets can be found in the folder you downloaded with MOSAIC, called "mosaic_example_data_rna.csv" and "mosaic_example_data_protein.csv", respectively. If you have unevenly sampled data, choose the smallest resolution and leave all missing column samples blank.

Note that MOSAIC is not recommended for 6 hour+ resolutions with only one replicate, as this results in high false discovery rates. At least 2 replicates for 6 hour+ resolutions are recommended.

## MOSAIC R Package

MOSAIC's methodology is now available as an R package on CRAN! To download and use MOSAIC as a package, enter the following in the R console:

```r
install.packages("mosaic.find")
library(mosaic.find)
```

With this, you then have access to finding rhythms in data with one function, mosaic_find(). For more information on how to use this package and its functionality, check out the [mosaic.find vignette](https://cran.r-project.org/web/packages/mosaic.find/vignettes/mosaic-vignette.html).

Note that using this package requires knowledge of coding in R. If you having no coding knowledge, we recommend that you download and use the app as directed above. Also note that this version of MOSAIC does not take advantage of parallelism that the MOSAIC app does and therefore takes longer to run (only using one core, rather than using all cores except one). Further, there is no console output to show a progress bar of how long the output will take. If you would prefer built in parallelism and a progress bar, we highly recommend that you use the app. However, we have thought of several workarounds -- if interested, feel free to contact us with the "Feedback" form below.

## Minimum Version Information

Minimum versions for packages and sytems used in ECHO are the following:

| Package        | Minimum Version |
| -------------: |-------------|
| R | >= 3.5.1 |
| rstudioapi | >= 0.8|
| shiny | >= 1.3.2 |
| ggplot2 | >= 3.1.0 |
| VennDiagram | >= 1.6.20 |
| reshape2 | >= 1.4.3 |
| minpack.lm | >= 1.2-1|
| doParallel | >= 1.0.14|
| foreach | >= 1.4.4|
| interators | >= 1.0.10|
| doSNOW | >= 1.0.16|
| colorRamps | >= 2.3|
| dplyr | >= 0.8.3|
| ggplotify | >= 0.0.4|
| gridExtra | >= 2.3|

## Contact Information and Bug Reporting

As you may have noticed, this is still in beta testing! Therefore, we would love to hear your feedback on the program. For general feedback, email hurlej2@rpi.edu with the subject line "MOSAIC Feedback".

If you run into any errors, please email hurlej2@rpi.edu with the following (subject line: "MOSAIC Error"): 
- a short desciption of your problem
- MOSAIC version number 
- your dataset/file(s) (this may be a sample of at least 50% of the data)
- your exact settings for the run (a screenshot will do) 
- your exact error from the console window (a screenshot will do)

However, *please* read the FAQ below and all given information (including instructions in the app, example data, etc.) before sending error reports.

Contact:
Jennifer Hurley /
email: hurlej2@rpi.edu /
Rensselaer Polytechnic Institute

## FAQ

**Q:** My dataset isn't working for some reason and I'm not using the latest MOSAIC version! Why?

**A:** Please update to the current MOSAIC version (you can see this at the top of the github repository by checking out what the latest commit was), since this may have been corrected in an update. If this problem persists, feel free to send another email!

---

**Q:** Does MOSAIC work with 24-hour time course length data to determine circadian rhythms?

**A:** Yes, it does! However, I would not categorize any rhythms convincingly "Damped", "Harmonic", or "Forced", since you need more than one cycle to determine this.

---

**Q:** Does MOSAIC work with other omic types of data (such as metabolomics, phosphoproteomics, etc.)?

**A:** Yes, it does! Upload these in the slots for "RNA" and "protein" data in the same format as transcriptomic and proteomic data. However, take note of which slots they are uploaded in, as they will be labeled accordingly as "RNA" and "protein" regardless. Further, one should upload the omics type earlier in the workflow in the "RNA" slot. For example, with proteomic and phosphoproteomic data, you should upload the proteomic data in the "RNA" slot, and the phosphoproteomic data in the "protein" slot.

---

**Q:** Some of my genes/proteins are not showing up in the MOSAIC results! Why?

**A:** Genes and proteins must be named *exactly* the same in each dataset to use this method. If your gene does not have a matching protein, or vice versa, it will not be run and will be automatically removed from the results.

---

**Q:** My data has starting points/ending points/resolution of less than an hour, or a fractional amount of hours! How do I run this through MOSAIC?

**A:** If you have resolution of less than an hour, please enter the fraction into the box, in the form: numerator/denominator. For example, if my resolution was every 10 minutes (or 6 times every hour), I would enter: 1/6. This fractional form extends to starting and ending points as well. You must enter the fraction, NOT a mixed number. For example, if my starting time was 16 hours, 10 minutes, my starting time would be: 97/6. (This stems from the following calculation: (6/6 x 16) + (1/6))

---

**Q:** I was running MOSAIC, and it suddenly went grey! What happened?

**A:** There was an error, the cause of which can be found in the console. Check through the FAQ to see if it has been addressed, or if it's an obvious error (such as not loading any data).

---

**Q:** I get the following error (or similar) when I try to view a Gene List in the Visualizations part of MOSAIC:

```r
DataTables warning: table id=DataTables_Table_0 - Requested unknown parameter '8' for row 0.
  For more information about this error, please see http://datatables.net/tn/4
```

**A:** This is a [bug](https://community.rstudio.com/t/data-table-issue-while-rendering-the-shiny-page-datatables-warning-table-id-datatables-table-0-requested-unknown-parameter/44016/3) with Data Tables in Shiny 1.4. To check your version of Shiny, enter the following in the console:
```r
packageVersion("shiny")
```
This should give you the Shiny version. If you have version 1.4, a quick fix is the following:

1. Open MOSAIC's mosaic_app.R file in RStudio.
2. Enter `install.packages("DT")` in the console, which will install the DT package.
3. Add `library(DT)` at the top of the mosaic_app.R script, on its own line.
4. Press ctrl/cmd+F, which will open the find and replace tool at the top of the script.
5. In the left box, which is the "find" box, enter `dataTableOutput`. In the right box, which is the "replace" box, enter `DT::dataTableOutput`. Then click the rightmost button, which says `All`.
6. After you have done step 3, in the left box, which is the "find" box, enter `renderDataTable`. In the right box, which is the "replace" box, enter `DT::renderDataTable`. Then click the rightmost button, which says `All`.
7. Press ctrl/cmd+S, which saves the mosaic_app.R file.

After you've completed all these steps, the problem should be fixed when you run it again!
