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
* Contact Information and Bug Reporting
* FAQ

## Overview

MOSAIC (Multi-Omics Selection with Amplitude Independent Criteria) is an R-powered application designed to find and visualize circadian and non-circadian trends in multi-omics time course data using model selection and joint modeling. To read more about this work and cite us, see [MOSAIC: A Joint Modeling Methodology for Combined Circadian and Non-Circadian Analysis of Multi-Omics Data](https://www.biorxiv.org/content/10.1101/2020.04.27.064147v1) by H. De los Santos, et al. (2020). To read more about this work one of the contributing models is based on, see [*ECHO: an Application for Detection and Analysis of Oscillators Identifies Metabolic Regulation on Genome-Wide Circadian Output*](https://doi.org/10.1093/bioinformatics/btz617) by H. De los Santos et al. (2019), published in *Bioinformatics*.

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

## ECHO Features

MOSAIC's interface is divided into two sections: **Finding Trends** and **Visualizing Results**.

Within the **Finding Rhythms** tab, you can upload your RNA and protein data (.csv) and enter its information, such as time point range, resolution (in hours), and amount and type of replicates. You can then choose from a variety of preprocessing steps including smoothing, removing unexpressed genes, and normalization. You can then download your results as both a CSV (for viewing) and a .RData (for visualizations).

In the **Visualizing Results** tab, simply upload the .RData file from your results and choose from several visualization and gene subset exploration options. You can explore subsets of data under the "Gene List" tab and sort by the various output parameters, such as Period or P-Value. You can also choose from a host of automatically-generated visualizations, including summary visualizations, heat maps, gene expression plots (with or without replicates visualized), and parameter density graphs (examples displayed below).

<p align="center">
<img src="MOSAIC App/www/MOSAIC_Gene_Expression_Plot.png" width="200" /> <img src="MOSAIC App/www/MOSAIC_Summary_Viz.png" width="250" /> <img src="MOSAIC App/www/MOSAIC_Parameter_Density_Plot.png" width="200" /> <img src="MOSAIC App/www/MOSAIC_Heat_Map.png" width="250" /> <img src="MOSAIC App/www/MOSAIC_Heat_Map_Comparison.png" width="250" />
</p>

## Contact Information and Bug Reporting

As you may have noticed, this is still in beta testing! Therefore, we would love to hear your feedback on the program. For general feedback, email delosh@rpi.edu with the subject line "MOSAIC Feedback".

If you run into any errors, please email delosh@rpi.edu with the following (subject line: "MOSAIC Error"): 
- a short desciption of your problem
- MOSAIC version number 
- your dataset/file(s) (this may be a sample of at least 50% of the data)
- your exact settings for the run (a screenshot will do) 
- your exact error from the console window (a screenshot will do)

However, *please* read the FAQ below before sending error reports.

Contact:
Hannah De los Santos /
email: delosh@rpi.edu /
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

**Q:** My data has starting points/ending points/resolution of less than an hour, or a fractional amount of hours! How do I run this through MOSAIC?

**A:** If you have resolution of less than an hour, please enter the fraction into the box, in the form: numerator/denominator. For example, if my resolution was every 10 minutes (or 6 times every hour), I would enter: 1/6. This fractional form extends to starting and ending points as well. You must enter the fraction, NOT a mixed number. For example, if my starting time was 16 hours, 10 minutes, my starting time would be: 97/6. (This stems from the following calculation: (1/6 x 16) + (1/6))
