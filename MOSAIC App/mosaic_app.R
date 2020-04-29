# MOSAIC App
# By Hannah De los Santos
# Originated on: 1/27/20

# clearing workspace ----

# clear workspace (sorry!)
rm(list=ls())

# version information ----

vers_mosaic <- "0.2.1"

# load libraries and outside functions ----

# increase input size
options(shiny.maxRequestSize=150*1024^2)

# libraries for overall
library(shiny)
library(rstudioapi)

# libraries for running mosaic
library(minpack.lm)
library(doParallel)
library(foreach)
library(iterators)
library(doSNOW)
library(dplyr)

# libraries for visualization
library(ggplot2)
library(VennDiagram)
library(ggplotify)
library(colorRamps)
library(reshape2)
library(gridExtra)

# set working directory - only works in RStudio (with rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# general data ----

# last number that is a parameter
end_num <- 33

# model information
mod_names <- c("Linear","Exponential","ECHO","ECHO Joint","ECHO Linear", "ECHO Linear Joint")
n_mod <- length(mod_names)

# parameter names
echo_lin_param_n <- c("Initial_Amplitude", "AC_Coefficient", "Oscillation_Type", "Radian_Frequency", "Period", "Phase_Shift", "Hours_Shifted", "Slope", "Equilibrium_Value")
echo_lin_joint_param_n <- c(paste0(echo_lin_param_n, "_RNA"), paste0(echo_lin_param_n,"_Protein"))
echo_param_n <- c("Initial_Amplitude", "AC_Coefficient", "Oscillation_Type", "Radian_Frequency", "Period", "Phase_Shift", "Hours_Shifted", "Equilibrium_Value")
echo_joint_param_n <- c(paste0(echo_param_n, "_RNA"), paste0(echo_param_n,"_Protein"))
exp_param_n <- c("Initial_Amplitude", "Growth_Rate", "Equilibrium_Value")
lin_param_n <- c("Slope", "Equilibrium_Value")

# put all the specific parameter names in a list
param_map <- list(
  "ECHO Linear" = echo_lin_param_n,
  "ECHO Linear Joint" = echo_lin_joint_param_n,
  "ECHO" = echo_param_n,
  "ECHO Joint" = echo_joint_param_n,
  "Exponential" = exp_param_n,
  "Linear" = lin_param_n#,
)

# color information
rna_blue_high <- "#268beb"
pro_red_high <- "#d61e1e"

rna_blue_low <- "#9ebfde"
pro_red_low <- "#d49d9d"

# making a function to create automatic palettes for rna and protein viz
rna_col_func <- colorRampPalette(c(rna_blue_low, rna_blue_high))
pro_col_func <- colorRampPalette(c(pro_red_low, pro_red_high))

# auxiliary functions -----

# function to get the order for the heat map
# inputs:
# fin_sub: subset of final.df, result of running MOSAIC
# omic_type: full name of the omic, RNA or Protein
# outputs:
# vector, gene names in order for heat map
order_heat_map <- function(fin_sub, omic_type){
  # we are building up to final order
  name_ord <- c()

  # take the echo ones first
  # get all the echo-type models
  sub_df <- fin_sub[fin_sub[,paste0("Best_Model_",omic_type)] %in% c("ECHO","ECHO Joint","ECHO Linear", "ECHO Linear Joint"),]
  if (nrow(sub_df) > 0){
    # order by hours shifted
    name_ord <- c(name_ord, sub_df$Gene_Name[order(sub_df[,paste0("Hours_Shifted_",omic_type)])])
  }

  # and now we order the exponential
  sub_df <- fin_sub[fin_sub[,paste0("Best_Model_",omic_type)]  == "Exponential",]
  if (nrow(sub_df) > 0){
    # subset order by growth rate, then amplitude (both negative to positive)
    # get all the subsets of growth rate and amplitudes
    ord_list <- list(
      "g_neg_a_neg" = sub_df[,paste0("Growth_Rate_",omic_type)] < 0 & sub_df[,paste0("Initial_Amplitude_",omic_type)] < 0 & !is.na(sub_df[,paste0("Growth_Rate_",omic_type)]),
      "g_neg_a_pos" = sub_df[,paste0("Growth_Rate_",omic_type)] < 0 & sub_df[,paste0("Initial_Amplitude_",omic_type)] >= 0 & !is.na(sub_df[,paste0("Growth_Rate_",omic_type)]),
      "g_pos_a_neg" = sub_df[,paste0("Growth_Rate_",omic_type)] >= 0 & sub_df[,paste0("Initial_Amplitude_",omic_type)] < 0 & !is.na(sub_df[,paste0("Growth_Rate_",omic_type)]),
      "g_pos_a_pos" = sub_df[,paste0("Growth_Rate_",omic_type)] >= 0 & sub_df[,paste0("Initial_Amplitude_",omic_type)] >= 0 & !is.na(sub_df[,paste0("Growth_Rate_",omic_type)])
    )

    # order each of the taken subsets by growth rates
    for (i in 1:length(ord_list)){
      sub_name <- sub_df$Gene_Name[ord_list[[i]]]
      name_ord <- c(name_ord, sub_name[order(sub_df[ord_list[[i]],paste0("Growth_Rate_",omic_type)])])
    }
  }

  # and now we order the linear
  sub_df <- fin_sub[fin_sub[,paste0("Best_Model_",omic_type)]  == "Linear",]
  if (nrow(sub_df) > 0){
    # order by slope
    name_ord <- c(name_ord, sub_df$Gene_Name[order(sub_df[,paste0("Slope_",omic_type)])])
  }

  # https://stackoverflow.com/questions/10827300/matching-up-two-vectors-in-r
  # match(x,y), the ith element of the output is the first index of y that matches x[i], unless x[i] doesn't appear in y, in which case it gives NA
  # reverse it for the heat map
  return(rev(match(name_ord, fin_sub$Gene_Name)))
}

# function to get the order for the heat map comparison
# inputs:
# fin_sub: subset of final.df, result of running MOSAIC
# omic_type: full name of the omic, RNA or Protein
# outputs:
# vector, gene names in order for heat map comparison
order_heat_map_comparison <- function(fin_sub, omic_type){
  # we are building up to final order
  name_ord <- c()

  # save the nas for later
  which_na <- is.na(fin_sub[,paste0("Best_Model_",omic_type)])

  # take the echo ones first
  sub_df <- fin_sub[fin_sub[,paste0("Best_Model_",omic_type)] %in% c("ECHO","ECHO Joint","ECHO Linear", "ECHO Linear Joint") & !which_na,]
  if (nrow(sub_df) > 0){
    # order by hours shifted
    name_ord <- c(name_ord, sub_df$Gene_Name[order(sub_df[,paste0("Hours_Shifted_",omic_type)])])
  }

  # and now we order the exponential
  sub_df <- fin_sub[fin_sub[,paste0("Best_Model_",omic_type)]  == "Exponential" & !which_na,]
  if (nrow(sub_df) > 0){
    # order by growth rate, then amplitude
    ord_list <- list(
      "g_neg_a_neg" = sub_df[,paste0("Growth_Rate_",omic_type)] < 0 & sub_df[,paste0("Initial_Amplitude_",omic_type)] < 0 & !is.na(sub_df[,paste0("Growth_Rate_",omic_type)]),
      "g_neg_a_pos" = sub_df[,paste0("Growth_Rate_",omic_type)] < 0 & sub_df[,paste0("Initial_Amplitude_",omic_type)] >= 0 & !is.na(sub_df[,paste0("Growth_Rate_",omic_type)]),
      "g_pos_a_neg" = sub_df[,paste0("Growth_Rate_",omic_type)] >= 0 & sub_df[,paste0("Initial_Amplitude_",omic_type)] < 0 & !is.na(sub_df[,paste0("Growth_Rate_",omic_type)]),
      "g_pos_a_pos" = sub_df[,paste0("Growth_Rate_",omic_type)] >= 0 & sub_df[,paste0("Initial_Amplitude_",omic_type)] >= 0 & !is.na(sub_df[,paste0("Growth_Rate_",omic_type)])
    )

    # order each taken subset by growth rate
    for (i in 1:length(ord_list)){
      sub_name <- sub_df$Gene_Name[ord_list[[i]]]
      name_ord <- c(name_ord, sub_name[order(sub_df[ord_list[[i]],paste0("Growth_Rate_",omic_type)])])
    }
  }

  # and now we order the linear
  sub_df <- fin_sub[fin_sub[,paste0("Best_Model_",omic_type)] == "Linear" & !which_na,]
  if (nrow(sub_df) > 0){
    # order by slope
    name_ord <- c(name_ord, sub_df$Gene_Name[order(sub_df[,paste0("Slope_",omic_type)])])
  }

  # we put the empty ones at the end
  sub_df <- fin_sub[which_na,]
  if (nrow(sub_df) > 0){
    name_ord <- c(name_ord, sub_df$Gene_Name)
  }

  # https://stackoverflow.com/questions/10827300/matching-up-two-vectors-in-r
  # match(x,y), the ith element of the output is the first index of y that matches x[i], unless x[i] doesn't appear in y, in which case it gives NA
  # reverse it for the heat map
  return(name_ord)
}

# function to make the heat map matrix to visualize
# inputs:
# dat: expression data
# ord: order for data
# num_reps: number of replicates
# heat_subset_rep: "all", or the number of the replicate you want to viz
# timen: time points
# outputs:
# heat map matrix, to be visualized
make_hm_matrix <- function(dat, ord, num_reps, heat_subset_rep, timen){

  #get matrix of just the relative expression over time
  hm_mat <- as.matrix(dat)

  if (heat_subset_rep == "all"){ # make an average of replicates
    #if there are replicates, average the relative expression for each replicate
    mtx_reps <- list() # to store actual matrix
    mtx_count <- list() # to store how many are NA
    for (i in 1:num_reps){
      mtx_reps[[i]] <- hm_mat[, seq(i,ncol(hm_mat), by=num_reps)]
      mtx_count[[i]] <- is.na(mtx_reps[[i]])
      mtx_reps[[i]][is.na(mtx_reps[[i]])] <- 0
    }
    repmtx <- matrix(0L,ncol = length(timen),nrow = nrow(hm_mat))+num_reps # to store how many we should divide by
    hm_mat <- matrix(0L,ncol = length(timen),nrow = nrow(hm_mat)) # to store the final result
    for (i in 1:num_reps){
      hm_mat <- hm_mat + mtx_reps[[i]] # sum the replicates
      repmtx <- repmtx - mtx_count[[i]] # how many replicates are available for each time point
    }
    repmtx[repmtx==0] <- NA # to avoid division by 0 and induce NAs if there are no time points available
    hm_mat <- hm_mat/repmtx

  } else { # display only one time point
    # subset the heat map
    rep_looking_at <- as.numeric(heat_subset_rep)
    hm_mat <- hm_mat[,seq(rep_looking_at,ncol(hm_mat),by=num_reps)]
  }

  # center rows around mean
  # vector of row means
  all_row_mean <- rowMeans(hm_mat, na.rm = TRUE)
  # center heat map
  hm_mat <- hm_mat - all_row_mean

  #normalize each row to be between -1 and 1
  if (length(ord) > 0){
    gene_max <- suppressWarnings(apply(abs(hm_mat), 1, function(x) max(x, na.rm = T)))
    hm_mat <- if(gene_max != 0){ hm_mat/gene_max } else {hm_mat}

    #sort by phase shift
    if (nrow(hm_mat) > 1){
      hm_mat <- hm_mat[ord,]
    }
  }

  return(hm_mat)
}

# function to make the heat map plot
# inputs:
# hm_mat: heat map matrix
# omic_type: the full name of the omic (RNA, Protein)
# subset_look: string, subset we're visualizing
# if_comp: if a comparison heat map
# comp_y: vector, labels for y axis
# output:
# heat map plot
make_hm_plot <- function(hm_mat, omic_type, subset_look, if_comp = F, comp_y = c()){
  # melt for ease with ggplot
  tmp <- melt(t(hm_mat))

  # make heat map plot
  p <- ggplot(tmp, aes(Var1, factor(Var2, levels = unique(Var2)), fill = value))+
    geom_tile()+
    scale_fill_gradientn("", limits = c(-1,1), breaks = c(-1,0,1),colors = blue2yellow(256), na.value = "white")+
    theme(
      text = element_text(size = 20),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(hjust = .5)
    )+
    xlab("Time (Hours)")+
    ylab("Genes")+
    ggtitle(paste0(omic_type, ", ", subset_look, " Subset"))+
    NULL

  # if it's not a heat map comparison plot
  if (!if_comp){
    p <- p +
      theme(
        axis.ticks = element_blank(),
        axis.text = element_blank(),
      ) +
      NULL
  } else { # heat map comparison plot
    p <- p +
      theme(
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
      ) +
      scale_y_discrete(labels = comp_y)+
      NULL
  }

  return(p)
}

# UI ----

ui <- navbarPage(
  "MOSAIC: Multi-Omics Selection with Amplitude Independent Criteria for Circadian Data",

  # ui for finding trends ----
  tabPanel("Finding Trends", {
    sidebarLayout(
      sidebarPanel(
        # instructions
        tags$p(HTML("<b>Note: Instructions can be found in the 'Instructions' tab.</b> Required fields, with the exception of example data, are marked with *. Required fields in all cases are marked with **.")),

        tags$p(paste("MOSAIC Version:", vers_mosaic)),

        hr(),

        # upload data
        div(style="display: inline-block; vertical-align: center; width: 250px;",
            fileInput("rna", "Upload RNA Data File (.csv) *", accept=c(".csv"), width = "250px")),

        div(style="display: inline-block; vertical-align: top; width: 5px;",
            actionButton("rna_file_help", icon("question", lib="font-awesome"))),
        uiOutput("Help_file_rna"),

        div(style="display: inline-block; vertical-align: center; width: 250px;",
            fileInput("pro", "Upload Protein Data File (.csv) *", accept=c(".csv"), width = "250px")),

        div(style="display: inline-block; vertical-align: top; width: 5px;",
            actionButton("pro_file_help", icon("question", lib="font-awesome"))),
        uiOutput("Help_file_pro"),

        checkboxInput("use_example", "Use Example?", value = FALSE, width = NULL),

        # dataset properties
        textInput("title", "Project Title **"),

        div(style="display: inline-block; width: 80px;",
            textInput("begin", "Start Time *")),
        div(style="display: inline-block; width: 80px;",
            textInput("end", "End Time *")),
        div(style="display: inline-block; width: 90px;",
            textInput("resol", "Resolution *")),
        div(style="display: inline-block; vertical-align: top; width: 5px;",
            actionButton("time_help", icon("question", lib="font-awesome"))),
        uiOutput("Help_time"),

        div(style="display: inline-block;",
            radioButtons("tied", "Type of replicates? *", choices = NULL, selected = NULL,
                         inline = TRUE, width = NULL, choiceNames = c("None","Paired","Unpaired"), choiceValues = c("none","TRUE","FALSE"))),
        # one replicate defaults to paired replicates (smoothed without average)

        div(style="display: inline-block; vertical-align:top;  width: 20px;",
            actionButton("tying_help", icon("question", lib="font-awesome"))),
        uiOutput("Help_tying"),

        div(style="display: inline-block; width: 180px;",
            numericInput("num_reps_res", "Number of replicates? *", 1, min = 1)),

        div(class="header", checked=NA,
            tags$b("Looking for rhythms between:")),

        div(style="display: inline-block; width: 70px;",
            textInput("low", "(lower:)")),

        div(style="display: inline-block; width: 70px;",
            textInput("high", "(upper:)")),

        div(style="display: inline-block; vertical-align:top;  width: 20px;",
            actionButton("limit_help", icon("question", lib="font-awesome"))),
        uiOutput("Help_limit"),

        # preprocessing and preference choices
        div(style="display: inline-block;",
            checkboxInput("smooth", "Smooth data?", value = FALSE, width = NULL)),
        div(style="display: inline-block; vertical-align:top;  width: 20px;",
            actionButton("smooth_help", icon("question", lib="font-awesome"))),
        uiOutput("Help_smooth"),

        div(style="display: inline-block;",
            checkboxInput("rem_unexpr", "Remove unexpressed genes?", value = FALSE, width = NULL)),
        div(style="display: inline-block; width: 70px;",
            numericInput("rem_unexpr_amt_below", "Cutoff?", value = 0, step = 1, width = NULL)),
        div(style="display: inline-block; width: 95px;",
            numericInput("rem_unexpr_amt", "Threshold %?", value = 70,min = 0, max = 100, step = 1, width = NULL)),
        div(style="display: inline-block; vertical-align:top;  width: 20px;",
            actionButton("remove_help", icon("question", lib="font-awesome"))),tags$br(),
        uiOutput("Help_remove"),

        div(style="display: inline-block;",
            checkboxInput("is_normal", "Normalize data?", value = FALSE, width = NULL)),
        div(style="display: inline-block; vertical-align:top;  width: 20px;",
            actionButton("normal_help", icon("question", lib="font-awesome"))),tags$br(),
        uiOutput("Help_normal"),

        div(class="header", checked=NA,
            tags$b("Set AC Coefficient cutoffs to:")),

        div(style="display: inline-block; width: 80px;",
            textInput("over_cut", "OE/RE Cut", value = "0.15")),

        div(style="display: inline-block; width: 80px;",
            textInput("harm_cut", "HA Cut", value = "0.03")),

        div(style="display: inline-block; vertical-align:top;  width: 20px;",
            actionButton("cut_help", icon("question", lib="font-awesome"))),
        uiOutput("Help_cut"),

        # action buttons for running echo code and downloading results
        actionButton("find_trends", "Find Trends!"),tags$br(),
        hr(),
        downloadButton("downloadMOSAIC", "Download MOSAIC CSV"),
        downloadButton("downloadViz", "Download Visualization RData")
      ),
      mainPanel(
        tabsetPanel(
          # welcome text
          tabPanel("Welcome to MOSAIC",
            verbatimTextOutput("finish"),
            HTML('<center><img src="MOSAIC_Gene_Expression_Plot.png" style="width:500px"></center><br>'),
            tags$p(tags$b("Welcome to MOSAIC!")),
            tags$p("MOSAIC (Multi-Omics Selection with Amplitude Independent Criteria) is an app designed to find and identify both circadian (ECHO, ECHO Linear) and non-circadian (Linear, Exponential) trends and perform joint modeling in data from multiple omics types (such as transcriptomics and proteomics). Just upload CSVs of your RNA and protein data and enter its parameters to get started."),
            tags$p("For more information on the format of your data, please click the corresponding '?' buttons. Examples of this data are available in the MOSAIC application download folder. For explanations for what specific terms mean, click the '?' buttons to their right. Recommendations are also made in the '?' text. Please note: in order to download results, you must open this app in the browser window before starting to run."),
            tags$p("If you'd just like to try this out or get a good idea of what data format is necessary, just check 'Run Example'. This will run the example CSVs that came with your download of MOSAIC. This example data is fabricated expression data from 2 to 48 hours with 2 hour resolution and 3 replicates. Models were generated for each gene across the available model sets. Random and systematic missing data is also included, left blank."),
            tags$p("When you run your data, a progress bar will display in the bottom left corner showing the stage of progress. Another progress bar will display in the console window to show directly how far one is in the MOSAIC progress. Upon finishing, the results will display above. You can then download the MOSAIC results CSV and the results for use in the visualization tab (.RData)."),
            tags$p(HTML("If you are using the results from this app or want to learn about its methods, please <a href='https://www.biorxiv.org/content/10.1101/2020.04.27.064147v1'>cite us</a>.")),
            ("If you run into any errors, please email delosh@rpi.edu with the following (subject line: MOSAIC Error):"),tags$br(),
            "- a short desciption of your problem" ,tags$br(),
            "- MOSAIC version number",tags$br(),
            "- your dataset/file(s)",tags$br(),
            "- your exact settings for the run (a screenshot will do)",tags$br(),
            "- your exact error from the console window (a screenshot will do)",tags$br(),tags$br(),
            "All images created by ECHO using data from:",tags$br(),
            "Hurley, J. et al. 2014. PNAS. 111 (48) 16995-17002. Analysis of clock-regulated genes in Neurospora reveals widespread posttranscriptional control of metabolic potential. doi:10.1073/pnas.1418963111 ",tags$br(),
            "Jennifer M. Hurley, Meaghan S. Jankowski, Hannah De los Santos, Alexander M. Crowell, Samuel B. Fordyce, Jeremy D. Zucker, Neeraj Kumar, Samuel O. Purvine, Errol W. Robinson, Anil Shukla, Erika Zink, William R. Cannon, Scott E. Baker, Jennifer J. Loros, Jay C. Dunlap, Circadian Proteomic Analysis Uncovers Mechanisms of Post-Transcriptional Regulation in Metabolic Pathways, Cell Systems, Volume 7, Issue 6, 2018, Pages 613-626.e5, ISSN 2405-4712, https://doi.org/10.1016/j.cels.2018.10.014.",
            tags$br(),tags$br(),
            tags$p(paste("MOSAIC Version", vers_mosaic))
          ),

          # results instructions
          tabPanel("Interpreting MOSAIC Result Files",
            tags$h4(HTML("<b><u>Interpreting MOSAIC Result Files</b></u>")),
            tags$p("Once you run MOSAIC, there are two results files available for download: a MOSAIC CSV and an .RData file (used for the Visualizations tab). Below, you can find an explanation of the contents of the CSV file. You can check our published paper for more in-depth information on these parameters."),

            ("The MOSAIC CSV has the following columns:"),tags$br(),
            ("- Gene_Name: The name of each expression, as entered in the original data file."),tags$br(),
            ("- Best_Model (RNA or protein): The selected model type for the RNA or protein data, from a choice of oscillatory (ECHO, ECHO Joint ECHO Linear, ECHO Linear Joint) and non-oscillatory (Linear, Exponential) models."),tags$br(),
            ("- P_Value (RNA or protein): Significance of MOSAIC fit, unadjusted."),tags$br(),
            ("- BH_Adj_P_Value  (RNA or protein): Significance of MOSAIC fit, adjusted using the Benjamini-Hochberg criterion. Corrects for multiple hypothesis testing."),tags$br(),
            ("- P_Value_Joint: Significance of MOSAIC joint fit if best model is ECHO Joint or ECHO Linear Joint, unadjusted."),tags$br(),
            ("- BH_Adj_P_Value_Joint: Significance of MOSAIC joint fit if best model is ECHO Joint or ECHO Linear Joint, adjusted using the Benjamini-Hochberg criterion. Corrects for multiple hypothesis testing."),tags$br(),
            ("- P_Value_Linear_Slope (RNA or protein): Significance of linear slope of MOSAIC fit if the best model is Linear, unadjusted."),tags$br(),
            ("- BH_Adj_P_Value_Linear_Slope  (RNA or protein): Significance of linear slope of MOSAIC fit if the best model is Linear, adjusted using the Benjamini-Hochberg criterion. Corrects for multiple hypothesis testing."),tags$br(),
            ("- AC_Coefficient  (RNA or protein): Amplitude Change Coefficient. Parameter which states the amount of amplitude change over time in the system. Used in ECHO, ECHO Joint, ECHO Linear, and ECHO Linear Joint models."),tags$br(),
            ("- Oscillation_Type  (RNA or protein): States the expression's category based on forcing coefficient (forced, damped, harmonic). Used in ECHO, ECHO Joint, ECHO Linear, and ECHO Linear Joint models."),tags$br(),
            ("- Initial_Amplitude (RNA or protein): Parameter describing initial amplitude of expression. Used in Exponential, ECHO, ECHO Joint, ECHO Linear, and ECHO Linear Joint models."),tags$br(),
            ("- Radian_Frequency (RNA or protein): Parameter describing frequency of oscillations, in radians. Used in ECHO, ECHO Joint, ECHO Linear, and ECHO Linear Joint models."),tags$br(),
            ("- Period (RNA or protein): States the time for one complete oscillation, assumed to be in hours. Used in ECHO, ECHO Joint, ECHO Linear, and ECHO Linear Joint models."),tags$br(),
            ("- Phase_Shift (RNA or protein): Parameter describing the amount the oscillator is shifted, in radians. Used in ECHO, ECHO Joint, ECHO Linear, and ECHO Linear Joint models."),tags$br(),
            ("- Hours_Shifted (RNA or protein): Desribes the amount the oscillator is shifted in hours, calculated from phase shift and fitted period. This is the time of the first peak of the oscillation, relative to 0 as determined by the time course entered by the user. Used in ECHO, ECHO Joint, ECHO Linear, and ECHO Linear Joint models."),tags$br(),
            "- Growth_Rate (RNA or protein): Parameter describing the exponential change in amplitude. Used in Exponential models.",tags$br(),
            "- Slope (RNA or protein): Parameter describing the linear slope. Used in Linear, ECHO Linear, and ECHO Linear Joint models.",tags$br(),
            ("- Equilibrium_Value (RNA or protein): Parameter describing the center, i.e. the y-intercept at time point 0, as determined by the user supplied time course. Used in Linear, Exponential, ECHO, ECHO Joint, ECHO Linear, and ECHO Linear Joint models."),tags$br(),
            ("- Processed_TPX.R (RNA or protein): Your original data, after any selected preprocessing, for time point (TP) X, and replicate R. Example: TP1, TP1.1, TP2, TP2.1, etc."),tags$br(),
            ("- Fitted_TPX.R (RNA or protein): MOSAIC's fitted data for time point (TP) X, and replicate R."),tags$br(),tags$br(),

            tags$p("The .RData file contains a series of R objects that are necessary for the automatic visualizations on the next tab. These objects include the MOSAIC output and user input information.")
          )
        )
      )
    )
  }),

  # ui for visualizing multiomics ----
  tabPanel("Visualize Results", {
    sidebarLayout(
      sidebarPanel(
        div(class="header", checked=NA,
            tags$b("Required fields in all cases are marked with *.")),

        # upload data and specify data
        fileInput("viz_file","Upload MOSAIC Visualization File (.RData) *", accept=c(".RData", ".Rds")),

        div(style="display: inline-block; width: 90px;",
            textInput("start_range", "Start Period")),
        div(style="display: inline-block; width: 90px;",
            textInput("end_range", "End Period")),
        div(style="display: inline-block; vertical-align: top; width: 5px;",
            actionButton("time_range_help", icon("question", lib="font-awesome"))),
        uiOutput("Help_time_range"),

        # plotting preferences
        selectInput("om_type",
                    "How should omics data types be displayed in plots?",
                    c("Both, Overlaid", "Both, Side-by-Side", "Only RNA", "Only Protein")),

        selectInput("viz",
                    "Choose Visualization Method",
                    c("Summary","Gene Expression","Gene Expression with Replicates","Heat Map","Heat Map Comparison","Parameter Density Graph (PDG)")),

        textInput("gene_name",
                  "Enter Gene to View (for Gene Expression with/without Replicates)"),

        selectInput("pval_cat",
                    "Enter P-Value Adjustment to View (for PDG, Heat Maps, Gene Lists)",
                    c("BH Adj P-Value" = "BH_Adj_P_Value",
                      "P-Value" = "P_Value")),

        numericInput("pval_cutoff",
                     "Enter P-Value Significance Cutoff (for PDG, Heat Maps, Gene Lists)",
                     min = 0, max = 1, step = .01, value = 0.05),

        selectInput("coeff",
                    "Choose Coefficient to View (for PDG)",
                    c("Amplitude Change Coefficient" = "AC_Coefficient",
                      "Period" = "Period",
                      "Initial Amplitude" = "Initial_Amplitude",
                      "Radian Frequency" = "Radian_Frequency",
                      "Phase Shift" = "Phase_Shift",
                      "Hours Shifted" = "Hours_Shifted",
                      "Equilibrium Value" = "Equilibrium_Value",
                      "Growth Rate" = "Growth_Rate",
                      "Slope" = "Slope",
                      "P-Value" = "P_Value",
                      "BH Adj P-Value" = "BH_Adj_P_Value",
                      "P-Value Joint" = "P_Value_Joint",
                      "BH Adj P-Value Joint" = "BH_Adj_P_Value_Joint",
                      "P-Value Linear Slope" = "P_Value_Linear_Slope",
                      "BH Adj P-Value Linear Slope" = "BH_Adj_P_Value_Linear_Slope")),

        div(style="display: inline-block;",
            selectInput("subset_look","Subset of Data to View (for Heat Maps, PDG, Gene Lists)",
                        c("None",
                          "MOSAIC",
                          "Linear",
                          "Exponential",
                          "ECHO",
                          "ECHO Joint",
                          "ECHO All",
                          "ECHO Linear",
                          "ECHO Linear Joint",
                          "ECHO Linear All",
                          "Oscillatory",
                          "Non-Oscillatory",
                          "Forced All Oscillatory",
                          "Damped All Oscillatory",
                          "Harmonic All Oscillatory",
                          "Forced ECHO All",
                          "Damped ECHO All",
                          "Harmonic ECHO All",
                          "Forced ECHO Linear All",
                          "Damped ECHO Linear All",
                          "Harmonic ECHO Linear All"#,
                          # "No Deviation"
                          ))),

        div(style="display: inline-block; vertical-align:top;  width: 20px;",
            actionButton("subset_help", icon("question", lib="font-awesome"))),
        uiOutput("Help_subset"),

        div(style="display: inline-block;",
            textInput("heat_subset_rep","Enter Replicate to View (for Heat Maps)", value = "all")),
        div(style="display: inline-block; vertical-align:top;  width: 20px;",
            actionButton("heat_subset_rep_help", icon("question", lib="font-awesome"))),
        uiOutput("Help_heat_subset_rep"),

        div(class="header", checked=NA,
            tags$b("Enter Comparison Genes (line separated, for Heat Map Comparison):")),
        div(style="display: inline-block;",
            textAreaInput("rna_focus","RNA:",width="150px",height = "100px")
        ),
        div(style="display: inline-block;",
            textAreaInput("pro_focus","Protein:",width="150px",height = "100px")
        ),
        div(style="display: inline-block; vertical-align:top;  width: 20px;",
            actionButton("focus_help", icon("question", lib="font-awesome"))),
        uiOutput("Help_focus"),

        radioButtons("focus_ord","Dominant Ordering (for Heat Map Comparison):",
                     choices = c("RNA", "Protein"), inline = TRUE),

        div(class="header", checked=NA,
            tags$b("Visualization window/download size (in pixels):")),

        div(style="display: inline-block; width: 70px;",
            textInput("viz_wid", "width:")),

        div(style="display: inline-block; width: 70px;",
            textInput("viz_height", "height:")),

        br(),
        # action buttons for visualization and downloading results
        actionButton("go","Update Visuzalization"),
        hr(),
        downloadButton('downloadPlot', 'Download Plot')
      ),
      mainPanel(
        tabsetPanel(
          # visualization tab
          tabPanel("Visualization",
            uiOutput("plot_ui"),
            tags$br(),
            fluidRow(style = "border: 1px #e3e3e3; border-style: solid; border-radius: 10px; background: #f5f5f5; padding: 10px;",
                     uiOutput("viz_text")
            )
          ),

          # gene list tab
          tabPanel("Gene List",
            tabsetPanel(
              tabPanel("Both, Side-by-Side",
                fluidRow(
                  column(6, h2("RNA"), downloadButton('download_rna_side', 'Download RNA Gene List'), dataTableOutput("gene_list_rna1")),
                  column(6, h2("Protein"), downloadButton('download_pro_side', 'Download Protein Gene List'), dataTableOutput("gene_list_pro1"))
                )
              ),
              tabPanel("Both, Overlaid", downloadButton('download_overlaid', 'Download Gene List'), dataTableOutput("gene_list_overlaid")),
              tabPanel("Only RNA", downloadButton('download_rna', 'Download Gene List'), dataTableOutput("gene_list_rna")),
              tabPanel("Only Protein", downloadButton('download_pro', 'Download Gene List'), dataTableOutput("gene_list_pro"))
            )
          ),

          # user inputs tab
          tabPanel("User Inputs",
            fluidRow(style = "border: 1px #e3e3e3; border-style: solid; border-radius: 10px; background: #f5f5f5; padding: 10px;",
                    uiOutput("user_input_text")
            )
          ),

          # visualization instructions
          tabPanel("Instructions",
            tags$h4(HTML("<b><u><center>Instructions:</center></b></u>")),
            tags$p("Once you have run your data through the Finding Trends part of MOSAIC, you can visualize and explore your results using the .RData file available for download after your run. To update visualizations or gene lists, select options from each of the drop down menus to the left and press Update Visualization. There are several types of visualizations/explorations currently available:"),

            HTML('<center>'),tags$b("Gene Lists:"),HTML('</center>'),
            tags$p("Gene lists display gene names and parameter information for specified subsets in RNA and Protein. This information is displayed side-by-side for both RNA and protein, overlaid, or for each individually. Overlaid gene lists indicate a subset where criteria is true for RNA or protein."),

            HTML('<center><img src="MOSAIC_Summary_Viz.png" style="width:400px"></center><br>'),
            HTML('<center>'),tags$b("Summary Visualization:"),HTML('</center>'),
            tags$p("In the Summary Visualization, users can compare resulting significant model distributions from transcriptomics (RNA) and proteomics (protein). At the top of the visualization, the left Venn Diagram shows the overlap of significant models in RNA and protein, while the right shows the overlap of significant oscillatory models. Below this, a comparative bar graphs shows the percentages of each significant model relative to the total number of significant models in RNA or protein, respectively. Counts of significant models are also displayed on the bar graph alongside model and omics type names. These counts are also reflected below the visualization, along with additional categories such as AC coefficient categories."),

            HTML('<center><img src="MOSAIC_Gene_Expression_Plot.png" style="width:300px"></center><br>'),
            HTML('<center>'),tags$b("Gene Expression Plots:"),HTML('</center>'),
            tags$p("Gene Expression Plots display processed gene expression values and fitted models for a selected genes. The fitted model is displayed as the dark corresponding color, while original data is displayed either as a lighter shaded line or, if replicates exist, filled in between the maximum and minimum replicate value at each time point. The ability to plot replicates, if they exist, is available as well. Expressions for RNA and protein can be plotted overlaid with each other, side-by-side in separate plots, or plotted with only one of the other on the plot. Below the visualization, relevant model parameters and p-values are listed."),

            HTML('<center><img src="MOSAIC_Heat_Map.png" style="width:400px"></center><br>'),
            HTML('<center>'),tags$b("Heat Maps:"),HTML('</center>'),
            tags$p("Heat Maps visualize averaged expression for a selected subset of genes. The expression for each gene in the heat map is mean-centered and normalized by the absolute maximum value of the averaged replicate expression, such that each row is on a scale of [-1,1]. Individual replicates can also be visualized. High expression values are colored yellow, while low values are colored blue. Available subsets for the heat map include the MOSAIC subset, each model subset, and AC coefficient subsets, among others. Heat maps are ordered first by model type, then by specified parameters, in the following manner: oscillatory genes, sorted by phase; exponential genes, subset by growth rate and then amplitude, ordered by growth rate; linear genes, ordered by slope. Heat maps for RNA and protein can be displayed side-by-side or individually. Total amounts of genes for each heat map are displayed below the plot."),

            HTML('<center><img src="MOSAIC_Heat_Map_Comparison.png" style="width:400px"></center><br>'),
            HTML('<center>'),tags$b("Heat Map Comparison:"),HTML('</center>'),
            tags$p("Heat Map Comparison plots are an extension of heat maps, which provide a method for direct comparison of user supplied lists of genes. Users enter genes for comparison in both RNA and Protein, which are displayed in heat maps side-by-side, with white reflecting genes that were only specified in the opposite omics type subset. These heat maps are ordered by a user-specified dominant omics type. This applies heat map ordering described above to the dominant omics type, then uses this exact ordering on both omics types. Below the plot, gene names for those that appear in both and only each omics types are displayed."),

            HTML('<center><img src="MOSAIC_Parameter_Density_Plot.png" style="width:300px"></center><br>'),
            HTML('<center>'),tags$b("Parameter Density Graphs (PDG):"),HTML('</center>'),
            tags$p("In Parameter Density Graphs, density graphs for a specified parameter and subset are displayed, showing the range and distribution of the specified parameter. RNA and protein graphs can be overlaid, plotted side-by-side, or plotted separately. Below the plot, quantile information for each distribution appears. Graphs will only be plotted if there are at least two non-missing values for the selected subset.")
          )
        )
      )
    )
  })
)

# SERVER ----

server <- function(input, output, session) {
  # help text mosaic run ----

  output$Help_file_rna=renderUI({ # csv help
    if(input$rna_file_help%%2){
      helpText("This RNA data must be in a .csv file with the following format: first row is column labels, first column has gene labels/names (which correspond exactly to names in the Protein file), and all other columns have expression data. This expression data must be ordered by time point then by replicate, and must have evenly spaced time points. Any missing data must have cells left blank.")
    }
    else{
      return()
    }
  })

  output$Help_file_pro=renderUI({ # csv help
    if(input$pro_file_help%%2){
      helpText("This Protein data must be in a .csv file with the following format: first row is column labels, first column has gene labels/names (which correspond exactly to names in the RNA file), and all other columns have expression data. This expression data must be ordered by time point then by replicate, and must have evenly spaced time points. Any missing data must have cells left blank.")
    }
    else{
      return()
    }
  })

  output$Help_time=renderUI({ # time inputs help
    if(input$time_help%%2){
      helpText("These numbers indicate the beginning, end, and resolution (difference between time points) of the time points in your data (in hours). For example, data that begins at 2 hours, ends at 48 hours, with a resolution of 2 hours. If your data has points in a fractional amount of hours, please enter the fraction into the box, in the form: numerator/denominator.")
    }
    else{
      return()
    }
  })

  output$Help_smooth=renderUI({ # data smoothing help
    if(input$smooth_help%%2){
      helpText("Indicates whether data should be smoothed, which also depends on type of data. If checked, the data is weighted smoothed over a rolling window of 3 points, with each of the points having weights of 1,2,1 respectively. If paired data, each replicate is smoothed independently. If unpaired data, each time point is smoothed by itself, centered, and the average expression per time point on either side. Smoothed data will be returned in output files. Note: this will increase running time.")
    }
    else{
      return()
    }
  })

  output$Help_limit=renderUI({ # upper and lower limits for rhythms help
    if(input$limit_help%%2){
      helpText("Upper and lower limits to look for rhythms, in hours. For example, one could look for hours between 20 and 26. If limits left blank, rhythms of any length within timecourse will be considered.")
    }
    else{
      return()
    }
  })

  # TODO: COME BACK TO THIS
  output$Help_remove=renderUI({ #remove unexpressed genes help
    if(input$remove_help%%2){
      helpText("If checked, marks and does not fit genes that have low detection/expression in the dataset. A gene is considered adequately expressed for modeling if a absolute value above the cutoff is available for at least the specified threshold percentage of the total time points. A 70% threshold is recommended. By default, genes with a zero value at all time points are marked and not fit to allow appropriate calculation of multiple hypothesis corrections. Recommended.")
    }
    else{
      return()
    }
  })

  output$Help_tying=renderUI({
    if((input$tying_help)%%2){ # paired or unpaired help
      helpText("Indicates whether replicates are paired (replicates are related within a time course) or unpaired (replicates are not related to any specific time course).")
    }
    else{
      return()
    }
  })

  output$Help_remove=renderUI({ #remove unexpressed genes help
    if(input$remove_help%%2){
      helpText("If checked, marks and does not fit genes that have low detection/expression in the dataset. A gene is considered adequately expressed for modeling if an absolute value above the cutoff is available for at least the specified threshold percentage of the total time points. A 70% threshold is recommended. Note that removal of an unexpressed gene from one dataset will remove them from both datasets. By default, genes with less than 3 nonmissing values are marked unexpressed. Recommended.") #  By default, genes with a zero value at all time points are marked and not fit to allow appropriate calculation of multiple hypothesis corrections.
    }
    else{
      return()
    }
  })

  output$Help_normal=renderUI({
    if(input$normal_help%%2){ # normalize help
      helpText("Normalizes data by row using the normal distribution (subtract each row by row mean and divide by row standard deviation). Normalized data is returned in results, rather than original data. Normalized data will be returned in output files. Recommended for un-normalized data. If you have normalized your data in any way, this option is strongly not recommended.")
    }
    else{
      return()
    }
  })

  output$Help_cut=renderUI({ # ac coeff specification help
    if(input$cut_help%%2){
      helpText("If you would like different than standard AC coefficient cutoffs for Overexpressed/Repressed (OE/RE) or Harmonic (HA) values, enter changes here. Since these are symmetric, these will be both positive and negative cutoffs. For example, defaults indicate that OE expressions are below -0.15, RE are above 0.15, and HA falls between -0.03 and 0.03. Damped and forced cutoffs fall between those values accordingly. This will affect both ECHO and ECHO Linear models.")
    }
    else{
      return()
    }
  })

  # running mosaic ----

  observeEvent(input$find_trends, {
    withProgress(message = "Finding MOSAIC Trends!", value = 0, {
      # add the example here
      incProgress(1/2,detail = paste("Running MOSAIC. Started on:",Sys.time()))

      # functions for mosaic
      source("mosaic_master.R", local = T)

      # use example data?
      if (input$use_example){
        # load data
        rna <- read.csv("mosaic_example_data_rna.csv", header = T, stringsAsFactors = F)
        pro <- read.csv("mosaic_example_data_protein.csv", header = T, stringsAsFactors = F)

        # creating times sequence used for the genes
        begin <- 2 # beginning
        end <- 48 # end
        resol <- 2 # resolution of data
        timen <- seq(begin,end,resol) # the times for cicadian rhythms

        num_reps <- 3 # number of replicates
        tied <- TRUE

      } else {
        # load data
        rna <- read.csv(input$rna$datapath, header = T, stringsAsFactors = F)
        pro <- read.csv(input$pro$datapath, header = T, stringsAsFactors = F)

        # creating times sequence used for the genes
        begin <- as.numeric(sapply(input$begin, function(x) eval(parse(text=x)))) # beginning
        end <- as.numeric(sapply(input$end, function(x) eval(parse(text=x)))) # end
        resol <- as.numeric(sapply(input$resol, function(x) eval(parse(text=x)))) # resolution of data
        timen <- seq(begin,end,resol) # the times for cicadian rhythms

        num_reps <- as.numeric(input$num_reps_res) # number of replicates

        if (input$tied=="none"){ # one replicate, default to true paired-ness
          tied <- TRUE
        } else {# more than one replicate
          tied <- as.logical(input$tied) # the type of replicate
        }

      }

      # subset to only run the fits in both
      both_n <- intersect(pro[,1], rna[,1])
      genes_rna <- rna[rna[,1] %in% both_n,]
      genes_pro <- pro[pro[,1] %in% both_n,]
      # order them the same
      genes_rna <- genes_rna[order(genes_rna[,1]),]
      genes_pro <- genes_pro[order(genes_pro[,1]),]

      if (nrow(genes_rna)==1){
        # creating a constant row and adding it to genes
        add_row <- data.frame(matrix(0L, 1, ncol(genes_rna)))
        add_row[1,1] <- "not considering"
        colnames(add_row) <- colnames(genes_rna)
        genes_rna <- rbind(genes_rna,add_row)

        colnames(add_row) <- colnames(genes_pro)
        genes_pro <- rbind(genes_pro,add_row)

        add_one <- TRUE # marker for appropriate displays for progress bar
      } else{
        add_one <- FALSE # marker for appropriate displays for progress bar
      }

      # figuring out whether a range is wanted, adjusting accordingly
      if (input$low ==""){ # empty low input, adjust to time series
        if (resol >= 1){
          low <- 2*pi/resol
          low_end <- resol
        }
        else{ # if the begining is >=0, smallest period available is 1
          low <- 2*pi/1
          low_end <- 1
        }
      } else{ # there is a low input
        low_input <- as.numeric(sapply(input$low, function(x) eval(parse(text=x))))
        low <- 2*pi/low_input
        low_end <- low_input
      }
      if (input$high ==""){ # empty high input, adjust to time series
        high <- 2*pi/(resol*length(timen))
        high_end <- (resol*length(timen))
      } else{ # there is a high input
        high_input <- as.numeric(sapply(input$high, function(x) eval(parse(text=x))))
        high <- 2*pi/high_input
        high_end <- high_input
      }

      # bootstrapping parameters are for when you have bootstrapping

      # harmonic and overexpressed cutoffs
      harm_cut <- abs(as.numeric(sapply(input$harm_cut, function(x) eval(parse(text=x)))))
      over_cut <- abs(as.numeric(sapply(input$over_cut, function(x) eval(parse(text=x)))))

      # preprocessing ----

      # removing unexpressed genes
      rem_unexpr <- input$rem_unexpr # indicator for removing unexpressed genes
      rem_unexpr_amt <- (input$rem_unexpr_amt)/100 # threshold for removing unexpressed genes, converted to a decimal
      rem_unexpr_amt_below <- abs(input$rem_unexpr_amt_below)
      # if yes, check for genes that are unexpressed before preprocessing
      if (rem_unexpr){
        rem_unexpr_vect_rna <- genes_unexpressed_all(genes_rna,rem_unexpr_amt,rem_unexpr_amt_below)
        rem_unexpr_vect_pro <- genes_unexpressed_all(genes_pro,rem_unexpr_amt,rem_unexpr_amt_below)
      } else{
        rem_unexpr_vect_rna <- rem_unexpr_vect_pro <- rep(FALSE,nrow(genes_rna))
      }

      rem_unexpr_combine <- (rem_unexpr_vect_rna | rem_unexpr_vect_pro)

      # normalize and store original data
      if (input$is_normal){
        norm_list <- normalize_all(genes_rna)
        genes_rna <- norm_list$dat

        norm_list <- normalize_all(genes_pro)
        genes_pro <- norm_list$dat
      }

      # make averages
      avg_genes_rna <- avg_all_rep(num_reps, genes_rna, timen)
      avg_genes_pro <- avg_all_rep(num_reps, genes_pro, timen)

      # smoothing, if necessary
      if (input$smooth){
        is_weighted <- T

        if (tied){ # if paired replicates
          genes_rna <- smoothing_all_tied(is_weighted, num_reps, genes_rna, avg_genes_rna, timen)
          genes_pro <- smoothing_all_tied(is_weighted, num_reps, genes_pro, avg_genes_pro, timen)
        } else{ # if unpaired replicates
          genes_rna <- smoothing_all_untied(is_weighted, num_reps, genes_rna, avg_genes_rna, timen)
          genes_pro <- smoothing_all_untied(is_weighted, num_reps, genes_pro, avg_genes_pro, timen)
        }

        # need to retake the averages because new smoothed data
        avg_genes_rna <- avg_all_rep(num_reps, genes_rna, timen)
        avg_genes_pro <- avg_all_rep(num_reps, genes_pro, timen)

      }

      # start the engines ----

      start <- Sys.time()

      # creating the final data frame

      all_param <- c("AC_Coefficient", "Oscillation_Type", "Period", "Hours_Shifted", "Initial_Amplitude", "Radian_Frequency", "Phase_Shift", "Growth_Rate", "Slope", "Equilibrium_Value")

      final_df <- data.frame(matrix(NA, nrow = nrow(genes_rna), ncol = 1+32+(ncol(genes_rna)-1)+(ncol(genes_pro)-1)+(length(timen))*2))
      colnames(final_df)[1] <- c("Gene_Name")
      colnames(final_df)[2:33] <-
        c("Best_Model_RNA", "Best_Model_Protein",
          "P_Value_RNA", "BH_Adj_P_Value_RNA",
          "P_Value_Protein", "BH_Adj_P_Value_Protein",
          "P_Value_Joint", "BH_Adj_P_Value_Joint",
          "P_Value_Linear_Slope_RNA", "BH_Adj_P_Value_Linear_Slope_RNA",
          "P_Value_Linear_Slope_Protein", "BH_Adj_P_Value_Linear_Slope_Protein",
          paste0(all_param,"_RNA"),
          paste0(all_param,"_Protein")
        )
      colnames(final_df)[34:ncol(final_df)] <- c(paste0("Processed_RNA_",colnames(genes_rna[-1])),
                                                 paste0("Fitted_RNA_",paste0("TP",timen)),
                                                 paste0("Processed_Protein_",colnames(genes_pro[-1])),
                                                 paste0("Fitted_Protein_",paste0("TP",timen)))
      # add the data we already know
      # gene names
      final_df[,1] <- genes_rna[,1]
      # data
      final_df[34:(33+(length(timen)*num_reps))] <- genes_rna[,-1]
      final_df[(33+(length(timen)*num_reps)+length(timen)+1):(33+(length(timen)*num_reps)+length(timen) + (length(timen)*num_reps))] <- genes_pro[,-1]

      # prepare for parallelism
      cores <- detectCores() # dectect how many processors
      cl <- makeCluster(cores[1]-1) # not to overload your computer, need one for OS
      registerDoSNOW(cl)

      # making a progress bar
      print(paste("Percentage finished out of",
                  if (!add_one){nrow(genes_rna)} else {1},
                  "expression(s):"))
      pb <- txtProgressBar(max = nrow(genes_rna), style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)

      # debug
      # source("mosaic_master.R")
      # fin <- list()
      # for (current_gene in 1:10){
      #   if (current_gene%%100 == 0 ){
      #     print(current_gene)
      #   }
      #
      #   fin[[current_gene]] <- multi_pipeline(current_gene, timen, resol, num_reps, final_df[current_gene,], genes_rna, genes_pro, avg_genes_rna, avg_genes_pro, rem_unexpr_combine, harm_cut, over_cut)
      # }

      # as.vector(lsf.str()) gets all the function names in the global environment
      fin <- suppressWarnings(
        foreach (current_gene=1:nrow(genes_rna), .export = as.vector(lsf.str()), .packages=c('minpack.lm'),.options.snow = opts) %dopar% {
          multi_pipeline(current_gene, timen, resol, num_reps, final_df[current_gene,], genes_rna, genes_pro,  avg_genes_rna, avg_genes_pro, rem_unexpr_combine, harm_cut, over_cut)
        }
      )
      final_df <- plyr::rbind.fill(lapply(fin, `[[`, 1))
      # AIC df is not included in forward facing app
      # aic_rna <- dplyr::bind_rows(lapply(fin, `[[`, 2))
      # aic_pro <- dplyr::bind_rows(lapply(fin, `[[`, 3))

      close(pb)

      stopCluster(cl) # stop using the clusters

      # change data frame values back to numeric
      final_df[,c(14,16:24,26:33)] <- sapply(final_df[,c(14,16:24,26:33)], as.numeric)

      # adjust p-values
      final_df$BH_Adj_P_Value_RNA <- p.adjust(final_df$P_Value_RNA, method = "BH")
      final_df$BH_Adj_P_Value_Protein <- p.adjust(final_df$P_Value_Protein, method = "BH")
      final_df$BH_Adj_P_Value_Joint <- p.adjust(final_df$P_Value_Joint, method = "BH")
      final_df$BH_Adj_P_Value_Linear_Slope_RNA <- p.adjust(final_df$P_Value_Linear_Slope_RNA , method = "BH")
      final_df$BH_Adj_P_Value_Linear_Slope_Protein <- p.adjust(final_df$P_Value_Linear_Slope_Protein , method = "BH")

      # debug
      # aic_rna$best_mod_pval_orig_adj <- p.adjust(aic_rna$best_mod_pval_orig, method = "BH")
      # aic_pro$best_mod_pval_orig_adj <- p.adjust(aic_pro$best_mod_pval_orig, method = "BH")
      #
      # aic_rna$best_mod_pval_after_adj <- p.adjust(aic_rna$best_mod_pval_after, method = "BH")
      # aic_pro$best_mod_pval_after_adj <- p.adjust(aic_pro$best_mod_pval_after, method = "BH")

      if (add_one){
        final_df <- final_df[-nrow(final_df),]
      }

      en <- Sys.time()
      print(en-start)

      time.taken <- difftime(en,start,units = "mins")

      output$finish <- renderPrint({
        cat(paste0("Done! Finished on: ",Sys.time(),"\n"))
        cat(paste("MOSAIC Time:",time.taken,"minutes\n"))
      })

      # download results ----

      # prep saving the user inputs
      user_input <- list(
        "MOSAIC_end_date" = en,
        "rna_file_name"= if (!is.null(input$rna$name)){input$rna$name} else {"Example"},
        "pro_file_name"= if (!is.null(input$pro$name)){input$pro$name} else {"Example"},
        "begin"=begin,
        "end"=end,
        "resol"=resol,
        "tied"=tied,
        "is_smooth"=input$smooth,
        "rem_unexpr"= input$rem_unexpr,
        "rem_unexpr_amt"= input$rem_unexpr_amt,
        "rem_unexpr_amt_below"= input$rem_unexpr_amt_below,
        "is_normal"=input$is_normal,
        "harm_cut"=input$harm_cut,
        "over_cut"=input$over_cut,
        "v_num"= vers_mosaic
      )

      # download mosaic csv
      output$downloadMOSAIC <- downloadHandler(
        filename = paste0(input$title,"_MOSAIC.csv"),
        content = function(file){
          write.csv(final_df, file, row.names = F, na = "")
        }
      )

      # download mosaic
      output$downloadViz <- downloadHandler(
        filename = paste0(input$title,"_MOSAIC_Visualization.RData"),
        content = function(file){
          save(file = file,
               list = c("final_df","num_reps","low_end","high_end","timen","user_input"),
               file
               )
        }
      )

    })
  })

  # help text viz ----

  output$Help_subset=renderUI({
    if(input$subset_help%%2){ # subset help
      helpText("Identifies which subset of data to view for specified visualizations. MOSAIC includes all model types. Non-Oscillatory models include Linear and Exponential Models. Oscillatory models include ECHO and ECHO Linear Models. For AC Coefficient categories, these display the chosen category as given by both oscillatory models.")
    }
    else{
      return()
    }
  })

  output$Help_heat_subset_rep=renderUI({
    if(input$heat_subset_rep_help%%2){ # heat map replicates help
      helpText("Identifies which replicate of data to view for heat map visualizations. If 'all' is entered, average of all replicates will be displayed. Otherwise, enter number of replicate to display (1 through total number of replicates). Note: visualization of specific replicates only recommended for paired data.")
    }
    else{
      return()
    }
  })

  output$Help_time_range=renderUI({ # time inputs help
    if(input$time_range_help%%2){
      helpText("These numbers indicate the beginning and end of time range you would like to view in visualizations and gene lists, applied to only ECHO and ECHO Linear models. If nothing is entered, the entire time range of rhythms that was searched for will be displayed. If your data has points in a fractional amount of hours, please enter the fraction into the box, in the form: numerator/denominator.")
    }
    else{
      return()
    }
  })

  output$Help_focus=renderUI({ # time inputs help
    if(input$focus_help%%2){
      helpText("Gene names for RNA and Protein subsets, which are line separated. Must match names in dataset exactly.")
    }
    else{
      return()
    }
  })

  # viz ----

  observeEvent(input$go, {
    # load MOSAIC visualization data
    load(input$viz_file$datapath)

    # Preprocessing ----

    # necessary overall parameters
    dat_rna <- final_df[,(end_num+1):((end_num+1)+((length(timen)*num_reps)-1))]
    fit_rna <- final_df[,(end_num+(length(timen)*num_reps)+1):((end_num+(length(timen)*num_reps)+1)+(length(timen)-1))]
    dat_pro <- final_df[,(end_num+(length(timen)*num_reps)+1+length(timen)):(end_num+(length(timen)*num_reps*2)+1+length(timen)-1)]
    fit_pro <- final_df[,(end_num+(length(timen)*num_reps*2)+1+length(timen)):ncol(final_df)]

    # create subsets
    # significance - no NAs here
    sig_rna <- final_df[,paste0(input$pval_cat,"_RNA")] < input$pval_cutoff &
      !is.na(final_df[,paste0(input$pval_cat,"_RNA")])
    sig_pro <- final_df[,paste0(input$pval_cat,"_Protein")] < input$pval_cutoff &
      !is.na(final_df[,paste0(input$pval_cat,"_Protein")])
    # models - no NAs
    lin_rna <- final_df$Best_Model_RNA == "Linear"
    lin_pro <- final_df$Best_Model_Protein == "Linear"
    exp_rna <- final_df$Best_Model_RNA == "Exponential"
    exp_pro <- final_df$Best_Model_Protein == "Exponential"
    echo_all_rna <- final_df$Best_Model_RNA %in% c("ECHO","ECHO Joint")
    echo_all_pro <- final_df$Best_Model_Protein %in% c("ECHO","ECHO Joint")
    echo_rna <- final_df$Best_Model_RNA %in% c("ECHO")
    echo_pro <- final_df$Best_Model_Protein %in% c("ECHO")
    echo_joint_rna <- final_df$Best_Model_RNA %in% c("ECHO Joint")
    echo_joint_pro <- final_df$Best_Model_Protein %in% c("ECHO Joint")
    echo_lin_all_rna <- final_df$Best_Model_RNA %in% c("ECHO Linear","ECHO Linear Joint")
    echo_lin_all_pro <- final_df$Best_Model_Protein %in% c("ECHO Linear","ECHO Linear Joint")
    echo_lin_rna <- final_df$Best_Model_RNA %in% c("ECHO Linear")
    echo_lin_pro <- final_df$Best_Model_Protein %in% c("ECHO Linear")
    echo_lin_joint_rna <- final_df$Best_Model_RNA %in% c("ECHO Linear Joint")
    echo_lin_joint_pro <- final_df$Best_Model_Protein %in% c("ECHO Linear Joint")
    # osc/non-osc - no NAs
    osc_rna <- echo_all_rna | echo_lin_all_rna
    osc_pro <- echo_all_pro | echo_lin_all_pro
    nonosc_rna <- lin_rna | exp_rna
    nonosc_pro <- lin_pro | exp_pro
    # ac coeff echo - no NAs in the combination
    forced_echo_rna <- echo_all_rna & final_df$Oscillation_Type_RNA == "Forced"
    damped_echo_rna <- echo_all_rna & final_df$Oscillation_Type_RNA == "Damped"
    harmonic_echo_rna <- echo_all_rna & final_df$Oscillation_Type_RNA == "Harmonic"
    forced_echo_pro <- echo_all_pro & final_df$Oscillation_Type_Protein == "Forced"
    damped_echo_pro <- echo_all_pro & final_df$Oscillation_Type_Protein == "Damped"
    harmonic_echo_pro <- echo_all_pro & final_df$Oscillation_Type_Protein == "Harmonic"
    # ac coeff echo lin - no NAs in the combination
    forced_echo_lin_rna <- echo_lin_all_rna & final_df$Oscillation_Type_RNA == "Forced"
    damped_echo_lin_rna <- echo_lin_all_rna & final_df$Oscillation_Type_RNA == "Damped"
    harmonic_echo_lin_rna <- echo_lin_all_rna & final_df$Oscillation_Type_RNA == "Harmonic"
    forced_echo_lin_pro <- echo_lin_all_pro & final_df$Oscillation_Type_Protein == "Forced"
    damped_echo_lin_pro <- echo_lin_all_pro & final_df$Oscillation_Type_Protein == "Damped"
    harmonic_echo_lin_pro <- echo_lin_all_pro & final_df$Oscillation_Type_Protein == "Harmonic"
    # comparisons between RNA and Protein
    only_rna <- sig_rna & !sig_pro
    only_pro <- !sig_rna & sig_pro
    rna_and_pro <- sig_rna & sig_pro
    neither_rna_pro <- !sig_rna & !sig_pro
    # unexpressed
    unexpr <- final_df$Best_Model_RNA == "Unexpressed" # RNA and protein are same for this

    # subset list - this combines significance, since it's implicit for these
    rna_sub_list <- list(
      "None" = rep(T, nrow(final_df)),
      "MOSAIC" = sig_rna,
      "Linear" = sig_rna & lin_rna,
      "Exponential" = sig_rna & exp_rna,
      "ECHO" = sig_rna & echo_rna,
      "ECHO Joint" = sig_rna & echo_joint_rna,
      "ECHO All" = sig_rna & echo_all_rna,
      "ECHO Linear" = sig_rna & echo_lin_rna,
      "ECHO Linear Joint" = sig_rna & echo_lin_joint_rna,
      "ECHO Linear All" = sig_rna & echo_lin_all_rna,
      "Oscillatory" = sig_rna & osc_rna,
      "Non-Oscillatory" = sig_rna & osc_rna,
      "Forced All Oscillatory" = sig_rna & (forced_echo_lin_rna | forced_echo_rna),
      "Damped All Oscillatory" = sig_rna & (damped_echo_lin_rna | damped_echo_rna),
      "Harmonic All Oscillatory" = sig_rna & (harmonic_echo_lin_rna | harmonic_echo_rna),
      "Forced ECHO All" = sig_rna & (forced_echo_rna),
      "Damped ECHO All" = sig_rna & (damped_echo_rna),
      "Harmonic ECHO All" = sig_rna & (harmonic_echo_rna),
      "Forced ECHO Linear All" = sig_rna & (forced_echo_lin_rna),
      "Damped ECHO Linear All" = sig_rna & (damped_echo_lin_rna),
      "Harmonic ECHO Linear All" = sig_rna & (harmonic_echo_lin_rna),
      "No Deviation"  = rep(F, nrow(final_df)) # TODO: NEED TO FIX
    )

    pro_sub_list <- list(
      "None" = rep(T, nrow(final_df)),
      "MOSAIC" = sig_pro,
      "Linear" = sig_pro & lin_pro,
      "Exponential" = sig_pro & exp_pro,
      "ECHO" = sig_pro & echo_pro,
      "ECHO Joint" = sig_pro & echo_joint_pro,
      "ECHO All" = sig_pro & echo_all_pro,
      "ECHO Linear" = sig_pro & echo_lin_pro,
      "ECHO Linear Joint" = sig_pro & echo_lin_joint_pro,
      "ECHO Linear All" = sig_pro & echo_lin_all_pro,
      "Oscillatory" = sig_pro & osc_pro,
      "Non-Oscillatory" = sig_pro & osc_pro,
      "Forced All Oscillatory" = sig_pro & (forced_echo_lin_pro | forced_echo_pro),
      "Damped All Oscillatory" = sig_pro & (damped_echo_lin_pro | damped_echo_pro),
      "Harmonic All Oscillatory" = sig_pro & (harmonic_echo_lin_pro | harmonic_echo_pro),
      "Forced ECHO All" = sig_pro & (forced_echo_pro),
      "Damped ECHO All" = sig_pro & (damped_echo_pro),
      "Harmonic ECHO All" = sig_pro & (harmonic_echo_pro),
      "Forced ECHO Linear All" = sig_pro & (forced_echo_lin_pro),
      "Damped ECHO Linear All" = sig_pro & (damped_echo_lin_pro),
      "Harmonic ECHO Linear All" = sig_pro & (harmonic_echo_lin_pro),
      "No Deviation"  = rep(F, nrow(final_df)) # TODO: NEED TO FIX
    )

    # then: subset to periods specified (if specified)
    criteria_start_rna <- criteria_start_pro <-
      criteria_end_rna <- criteria_end_pro <- rep(T, nrow(final_df))
    if (input$start_range != ""){
      strt <- as.numeric(sapply(input$start_range, function(x) eval(parse(text=x))))
      criteria_start_rna <- final_df[,"Period_RNA"] >= strt
      crtieria_start_rna[is.na(criteria_start_rna)] <- T

      criteria_start_pro <- final_df[,"Period_Protein"] >= strt
      crtieria_start_pro[is.na(criteria_start_pro)] <- T
    }

    if (input$end_range != ""){
      en <- as.numeric(sapply(input$end_range, function(x) eval(parse(text=x))))
      criteria_end_rna <- final_df[,"Period_RNA"] >= en
      crtieria_end_rna[is.na(criteria_end_rna)] <- T

      criteria_end_pro <- final_df[,"Period_Protein"] >= en
      crtieria_end_pro[is.na(criteria_end_pro)] <- T
    }

    # Gene List Creation ----

    # if it's true in either of them
    gene_list_overlaid <- final_df[(rna_sub_list[[input$subset_look]] | pro_sub_list[[input$subset_look]]) &
                                     (criteria_start_rna | criteria_start_pro) &
                                     (criteria_end_rna | criteria_end_pro),1:33]

    # only true in rna
    gene_list_rna <- final_df[(rna_sub_list[[input$subset_look]]) &
                                (criteria_start_rna) &
                                (criteria_end_pro),
                              c(1, which(grepl("RNA",colnames(final_df)[1:33])) )]

    # only true in protein
    gene_list_pro <- final_df[(pro_sub_list[[input$subset_look]]) &
                                (criteria_start_pro) &
                                (criteria_end_pro),
                              c(1, which(grepl("Protein",colnames(final_df)[1:33])) )]

    # render specified tables
    output$gene_list_overlaid <- renderDataTable({
      gene_list_overlaid
    }, options = list(scrollX = T))

    output$gene_list_rna1 <- output$gene_list_rna <- renderDataTable({
      gene_list_rna
    }, options = list(scrollX = T))

    output$gene_list_pro1 <- output$gene_list_pro <- renderDataTable({
      gene_list_pro
    }, options = list(scrollX = T))

    # download gene lists
    output$download_rna_side <- output$download_rna <- downloadHandler(
      filename = function() { paste("Gene_List_RNA_", input$subset_look, '_Subset_',input$om_type,'_MOSAIC.csv', sep='') },
      content = function(file) {
        write.csv(gene_list_rna, file, row.names = FALSE)
      }
    )

    output$download_pro_side <- output$download_pro <- downloadHandler(
      filename = function() { paste("Gene_List_Protein_", input$subset_look, '_Subset_',input$om_type,'_MOSAIC.csv', sep='') },
      content = function(file) {
        write.csv(gene_list_pro, file, row.names = FALSE)
      }
    )

    output$download_overlaid <- downloadHandler(
      filename = function() { paste("Gene_List_Overlaid_", input$subset_look, '_Subset_',input$om_type,'_MOSAIC.csv', sep='') },
      content = function(file) {
        write.csv(gene_list_overlaid, file, row.names = FALSE)
      }
    )
    # User Inputs ----

    # render user inputs for file
    output$user_input_text <- renderUI({
      HTML(
        paste0(
          "<h4><b><u>User Inputs:</h4></b></u>",
          "<b>- MOSAIC End Date and Time: </b>", user_input$MOSAIC_end_date,"</br>",
          "<b>- RNA File Name: </b>", user_input$rna_file_name,"</br>",
          "<b>- Protein File Name: </b>", user_input$pro_file_name,"</br>",
          "<b>- Time Course Start: </b>", user_input$begin,"</br>",
          "<b>- Time Course End: </b>", user_input$end,"</br>",
          "<b>- Time Course Resolution: </b>", user_input$resol,"</br>",
          "<b>- Number of Replicates: </b>", num_reps,"</br>",
          "<b>- Seeking Rhythms, Low End: </b>", low_end,"</br>",
          "<b>- Seeking Rhythms, High End: </b>", high_end,"</br>",
          "<b>- Paired Repicates: </b>", user_input$tied,"</br>",
          "<b>- Smoothing?: </b>", user_input$is_smooth,"</br>",
          "<b>- Remove unexpressed genes?: </b>", user_input$rem_unexpr,"</br>",
          "<b>- Remove unexpressed genes, percentage: </b>", user_input$rem_unexpr_amt,"</br>",
          "<b>- Remove unexpressed genes, cutoff: </b>", user_input$rem_unexpr_amt_below,"</br>",
          "<b>- Normalize data?: </b>", user_input$is_normal,"</br>",
          "<b>- Harmonic cutoff: </b>", user_input$harm_cut,"</br>",
          "<b>- Overexpressed/Repressed cutoff: </b>", user_input$over_cut,"</br>",
          "<b>- MOSAIC Version No.: </b>", user_input$v_num,"</br>"
        )
      )
    })

    # Visualizations ----

    if (input$viz == "Summary"){
      # Summary ----

      # visualization, above

      # venn diagram comparing rna and protein, only significant
      tp1 <- as.ggplot(grid.arrange(as.ggplot(grobTree(
        draw.pairwise.venn(area1 = sum(sig_rna),
                           area2 = sum(sig_pro),
                           cross.area = sum(sig_rna & sig_pro),
                           category = c("RNA","Protein"),
                           fill = c(rna_blue_high,pro_red_high),
                           lty = rep("blank",2),
                           alpha = rep(.8,2),
                           cat.pos = c(0, 0),
                           cat.dist=rep(.025,2),
                           cat.fontfamily = rep("sans",2),
                           fontfamily ="sans",
                           cex = rep(3, 3),
                           cat.cex = rep(3, 2),
                           margin = .05
                           )
      )), top = textGrob("Significant Trends in RNA and Protein", gp=gpar(fontface="bold"))))

      # venn diagram comparing rna and protein, oscillatory
      tp2 <- as.ggplot(grid.arrange(as.ggplot(grobTree(
        draw.pairwise.venn(area1 = sum(sig_rna & osc_rna),
                           area2 = sum(sig_pro & osc_pro),
                           cross.area = sum(sig_rna & sig_pro & osc_rna & osc_pro),
                           category = c("RNA","Protein"),
                           fill = c(rna_blue_high,pro_red_high),
                           lty = rep("blank",2),
                           alpha = rep(.8,2),
                           cat.pos = c(0, 0),
                           cat.dist=rep(.025,2),
                           cat.fontfamily = rep("sans",2),
                           fontfamily ="sans",
                           cex = rep(3, 3),
                           cat.cex = rep(3, 2),
                           margin = .05
        )
      )), top = textGrob("Significant Oscillatory Trends in RNA and Protein", gp=gpar(fontface="bold"))))

      # put the venn diagrams together to make the first row of the plot
      tp <- as.ggplot(grid.arrange(tp1,tp2, ncol = 2))

      # ggplot comparing model distributions in rna and protein
      # making the ggplot data frame
      gg.df <- data.frame(
        "count" = c(
            sum(lin_rna & sig_rna), sum(exp_rna & sig_rna), sum(echo_rna & sig_rna), sum(echo_joint_rna & sig_rna), sum(echo_lin_rna & sig_rna), sum(echo_lin_joint_rna & sig_rna),
            sum(lin_pro & sig_pro), sum(exp_pro & sig_pro), sum(echo_pro & sig_pro), sum(echo_joint_pro & sig_pro), sum(echo_lin_pro & sig_pro), sum(echo_lin_joint_pro & sig_pro)
          ),
        "perc" = c(
          sum(lin_rna & sig_rna)/sum(sig_rna), sum(exp_rna & sig_rna)/sum(sig_rna), sum(echo_rna & sig_rna)/sum(sig_rna), sum(echo_joint_rna & sig_rna)/sum(sig_rna), sum(echo_lin_rna & sig_rna)/sum(sig_rna), sum(echo_lin_joint_rna & sig_rna)/sum(sig_rna),
          sum(lin_pro & sig_pro)/sum(sig_pro), sum(exp_pro & sig_pro)/sum(sig_pro), sum(echo_pro & sig_pro)/sum(sig_pro), sum(echo_joint_pro & sig_pro)/sum(sig_pro), sum(echo_lin_pro & sig_pro)/sum(sig_pro), sum(echo_lin_joint_pro & sig_pro)/sum(sig_pro)
        )*100,
        "data" = c(rep(paste0("RNA (", sum(sig_rna),")"),n_mod), rep(paste0("Protein (", sum(sig_pro),")"),n_mod)),
        "model_type" = c(paste0(mod_names, "_RNA"),paste0(mod_names, "_Protein")),
        "model_name" = c(rep(mod_names,2)),
        stringsAsFactors = F
      )

      # create color palette for rna and protein
      col_pal <- c(rna_col_func(n_mod), pro_col_func(n_mod))
      names(col_pal) <- gg.df$model_type

      # create comparison plot, stacked bar -- x axis is percentage, counts are annotated
      # makes for a distribution of models
      bp <- ggplot(gg.df, aes(x = data, y= perc, fill = factor(model_type, levels = unique(model_type)), group = data))+
        geom_bar(stat = "identity")+
        geom_text(data = subset(gg.df, perc > 1),  aes(label = paste0(model_name, " (", count, ")")), size = 5, position = position_stack(vjust = .5), angle = 90)+
        theme_bw()+
        theme(
          text= element_text(size = 20),
          axis.text.y = element_text(angle = 90, hjust = .5),
          legend.position="none",
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(10,10,10,10), "pt")
        )+
        xlab("Model")+
        ylab("Percentage of Total Significant Terms")+
        scale_fill_manual(values = col_pal, name = "Model Type", breaks = names(col_pal))+
        ggtitle(paste0("Distribution of Significant Models"))+
        guides(fill = guide_legend(nrow = 2, byrow = T, title.position = "top", title.hjust = .5))+
        coord_flip()+
        scale_y_continuous(expand = c(0,0))+
        scale_x_discrete(expand = c(0,0))+
        NULL

      # make the final plot, venn diagrams on top, comparison bar plot on the bottom
      fin_plot <- as.ggplot(grid.arrange(tp,bp, heights = c(10,15)))

      # text, below: counts of different categories
      fin_text <-
        HTML(
          paste0(
            "<b>RNA:</b><br>",
            "<b>- Unexpressed: </b>",sum(unexpr),"<br>",
            "<b>- Significant: </b>",sum(sig_rna),"<br>",
            "<b>&emsp;  - Non-Oscillatory: </b>",sum(sig_rna & nonosc_rna),"<br>",
            "<b>&emsp;&emsp;    - Linear: </b>",sum(sig_rna & lin_rna),"<br>",
            "<b>&emsp;&emsp;    - Exponential: </b>",sum(sig_rna & exp_rna),"<br>",
            "<b>&emsp;  - Oscillatory: </b>",sum(sig_rna & osc_rna),"<br>",
            "<b>&emsp;&emsp;    - ECHO: </b>",sum(sig_rna & echo_all_rna)," (<b>Regular:</b> ", sum(sig_rna & echo_rna), ", <b>Joint: </b>", sum(sig_rna & echo_joint_rna), ")<br>",
            "<b>&emsp;&emsp;&emsp;      - Damped: </b>",sum(sig_rna & damped_echo_rna),"<br>",
            "<b>&emsp;&emsp;&emsp;      - Harmonic: </b>",sum(sig_rna & harmonic_echo_rna),"<br>",
            "<b>&emsp;&emsp;&emsp;      - Forced: </b>",sum(sig_rna & forced_echo_rna),"<br>",
            "<b>&emsp;&emsp;    - ECHO Linear: </b>",sum(sig_rna & echo_lin_all_rna)," (<b>Regular:</b> ", sum(sig_rna & echo_lin_rna), ", <b>Joint: </b>", sum(sig_rna & echo_lin_joint_rna), ")<br>",
            "<b>&emsp;&emsp;&emsp;      - Damped: </b>",sum(sig_rna & damped_echo_lin_rna),"<br>",
            "<b>&emsp;&emsp;&emsp;      - Harmonic: </b>",sum(sig_rna & harmonic_echo_lin_rna),"<br>",
            "<b>&emsp;&emsp;&emsp;      - Forced: </b>",sum(sig_rna & forced_echo_lin_rna),"<br>",
            "<br>",
            "<b>Protein:</b><br>",
            "<b>- Unexpressed: </b>",sum(unexpr),"<br>",
            "<b>- Significant: </b>",sum(sig_pro),"<br>",
            "<b>&emsp;  - Non-Oscillatory: </b>",sum(sig_pro & nonosc_pro),"<br>",
            "<b>&emsp;&emsp;    - Linear: </b>",sum(sig_pro & lin_pro),"<br>",
            "<b>&emsp;&emsp;    - Exponential: </b>",sum(sig_pro & exp_pro),"<br>",
            "<b>&emsp;  - Oscillatory: </b>",sum(sig_pro & osc_pro),"<br>",
            "<b>&emsp;&emsp;    - ECHO: </b>",sum(sig_pro & echo_all_pro)," (<b>Regular:</b> ", sum(sig_pro & echo_pro), ", <b>Joint: </b>", sum(sig_pro & echo_joint_pro), ")<br>",
            "<b>&emsp;&emsp;&emsp;      - Damped: </b>",sum(sig_pro & damped_echo_pro),"<br>",
            "<b>&emsp;&emsp;&emsp;      - Harmonic: </b>",sum(sig_pro & harmonic_echo_pro),"<br>",
            "<b>&emsp;&emsp;&emsp;      - Forced: </b>",sum(sig_pro & forced_echo_pro),"<br>",
            "<b>&emsp;&emsp;    - ECHO Linear: </b>",sum(sig_pro & echo_lin_all_pro)," (<b>Regular:</b> ", sum(sig_pro & echo_lin_pro), ", <b>Joint: </b>", sum(sig_pro & echo_lin_joint_pro), ")<br>",
            "<b>&emsp;&emsp;&emsp;      - Damped: </b>",sum(sig_pro & damped_echo_lin_pro),"<br>",
            "<b>&emsp;&emsp;&emsp;      - Harmonic: </b>",sum(sig_pro & harmonic_echo_lin_pro),"<br>",
            "<b>&emsp;&emsp;&emsp;      - Forced: </b>",sum(sig_pro & forced_echo_lin_pro),"<br>",
            "<br>",
            "<b>Comparing Significant RNA and Protein: </b><br>",
            "<b>- Only RNA: </b>",sum(only_rna),"<br>",
            "<b>- Only Protein: </b>",sum(only_pro),"<br>",
            "<b>- Both RNA and Protein: </b>",sum(rna_and_pro),"<br>",
            "<b>- Neither RNA nor Protein: </b>",sum(neither_rna_pro),"<br>"
          )
        )

    } else if (input$viz == "Gene Expression" | input$viz == "Gene Expression with Replicates"){
      # Gene Expression Plot ----

      if (input$gene_name == "" | !input$gene_name %in% final_df$Gene_Name){
        # if the gene name specified isn't in our dataset, don't plot/write anything
        fin_plot <- ggplot()+theme_bw()
        fin_text <- HTML("")
      } else {
        # basically we want to get the RNA and protein and plot the fit and the ribbon
        # get the logical of which row corresponds to the gene name
        gene_log <- final_df$Gene_Name == input$gene_name

        # form the data frame
        gg.df <- data.frame(
          "rna_fit" = as.numeric(fit_rna[gene_log,]),
          "pro_fit" = as.numeric(fit_pro[gene_log,]),
          "rna_max" = sapply(seq(1,ncol(dat_rna), by = num_reps), function(x) max(dat_rna[gene_log,x:(x+num_reps-1)], na.rm = T)),
          "pro_max" = sapply(seq(1,ncol(dat_pro), by = num_reps), function(x) max(dat_pro[gene_log,x:(x+num_reps-1)], na.rm = T)),
          "rna_min" = sapply(seq(1,ncol(dat_rna), by = num_reps), function(x) min(dat_rna[gene_log,x:(x+num_reps-1)], na.rm = T)),
          "pro_min" = sapply(seq(1,ncol(dat_pro), by = num_reps), function(x) min(dat_pro[gene_log,x:(x+num_reps-1)], na.rm = T)),
          "timen" = timen
        )

        # colors for lines
        col_vect <- c(
          "Original RNA" = rna_blue_low,
          "Original Protein" = pro_red_low,
          "Fit RNA" = rna_blue_high,
          "Fit Protein" = pro_red_high
        )
        # add replicate colors
        # add two to provide more variability
        rna_rep_col <- rna_col_func(num_reps+2)
        pro_rep_col <- pro_col_func(num_reps+2)

        # add replicates to df and color scale
        for (i in 1:num_reps){
          gg.df[,paste0("rna_rep",i)] <- as.numeric(dat_rna[gene_log,seq(i, ncol(dat_rna), by=num_reps)])
          gg.df[,paste0("pro_rep",i)] <- as.numeric(dat_pro[gene_log,seq(i, ncol(dat_pro), by=num_reps)])

          col_vect[paste0("RNA, Replicate ",i)] <- rna_rep_col[i+1]
          col_vect[paste0("Protein, Replicate ",i)] <- pro_rep_col[i+1]
        }

        # name breaks for the colors in the legend
        col_name_breaks <- c(names(col_vect)[grepl("RNA", names(col_vect), fixed = T)],
                             names(col_vect)[grepl("Protein", names(col_vect), fixed = T)])

        # base plot, filling in general information
        p <- ggplot(gg.df, aes(x = timen))+
          scale_fill_manual("", values = col_vect, breaks = col_name_breaks)+
          scale_color_manual("", values = col_vect, breaks = col_name_breaks)+
          ylab("Expression")+
          xlab("Time (Hours)")+
          theme_bw()+
          theme(
            text = element_text(size = 20),
            plot.title = element_text(hjust = .5),
            legend.position = "bottom",
            legend.direction = "horizontal"
          )+
          scale_x_continuous(expand = c(0, 0))+
          NULL

        if (input$om_type != "Both, Overlaid"){
          # if the plots are separate in some way for each omics type
          p_rna <- p_pro <- p

          if (num_reps > 1){
            # if there are replicates, add a ribbon to the plot
            p_rna <- p_rna +
              geom_ribbon(aes(ymax = rna_max, ymin = rna_min, fill = "Original RNA"), alpha = .5)+
              NULL

            p_pro <- p_pro +
              geom_ribbon(aes(ymax = pro_max, ymin = pro_min, fill = "Original Protein"), alpha = .5)+
              NULL
          } else {
            # if there are no replicates, add just a line to the plot
            p_rna <- p_rna +
              geom_line(aes(y = rna_rep1, color = "Original RNA"))+
              NULL

            p_pro <- p_pro +
              geom_line(aes(y = pro_rep1, color = "Original Protein"))+
              NULL
          }

          # if it's "gene expression w/replicates plot"
          if (grepl("Replicates", input$viz)){
            # add a line to the plot for each replicate
            for (i in 1:num_reps){
              p_rna <- p_rna +
                geom_line(aes_string(y = paste0("rna_rep",i), color = shQuote(paste0("RNA, Replicate ",i))))+
                NULL

              p_pro <- p_pro +
                geom_line(aes_string(y = paste0("pro_rep",i), color = shQuote(paste0("Protein, Replicate ",i))))+
                NULL
            }
          }

          # add the fitted data, for rna and protein separately
          p_rna <- p_rna +
            geom_line(aes(y = rna_fit, color = "Fit RNA"), size = 1)+
            ggtitle(paste0("RNA: ", final_df$Gene_Name[gene_log]))+
            guides(
              color = guide_legend(nrow = 2, byrow = F),
              fill = guide_legend(nrow = 2)
            )+
            NULL

          p_pro <- p_pro +
            geom_line(aes(y = pro_fit, color = "Fit Protein"), size = 1)+
            ggtitle(paste0("Protein: ", final_df$Gene_Name[gene_log]))+
            guides(
              color = guide_legend(nrow = 2, byrow = F),
              fill = guide_legend(nrow = 2)
            )+
            NULL

          # arrange the final plot based on the view type of omics
          fin_plot <-
            if (input$om_type == "Both, Side-by-Side"){
              as.ggplot(grid.arrange(p_rna, p_pro, ncol = 2))
            } else if (input$om_type == "Only RNA"){
              p_rna
            } else if (input$om_type == "Only Protein"){
              p_pro
            }
        } else {
          # otherwise, the omics types should be overlaid
          fin_plot <- p

          if (num_reps > 1){
            # if replicates, add a ribbon
            fin_plot <- fin_plot +
              geom_ribbon(aes(ymax = rna_max, ymin = rna_min, fill = "Original RNA"), alpha = .5)+
              geom_ribbon(aes(ymax = pro_max, ymin = pro_min, fill = "Original Protein"), alpha = .5)+
              NULL
          } else {
            # if no replicates, just put the original data in a line
            fin_plot <- fin_plot +
              geom_line(aes(y = rna_rep1, color = "Original RNA"))+
              geom_line(aes(y = pro_rep1, color = "Original Protein"))+
              NULL
          }

          # add the replicates to the plot, if specified
          if (grepl("Replicates", input$viz)){
            for (i in 1:num_reps){
              fin_plot <- fin_plot +
                geom_line(aes_string(y = paste0("rna_rep",i), color = shQuote(paste0("RNA, Replicate ",i))))+
                geom_line(aes_string(y = paste0("pro_rep",i), color = shQuote(paste0("Protein, Replicate ",i))))+
                NULL
            }
          }

          # add the fits to the plot
          fin_plot <- fin_plot +
            geom_line(aes(y = rna_fit, color = "Fit RNA"), size = 1)+
            geom_line(aes(y = pro_fit, color = "Fit Protein"), size = 1)+
            ggtitle(final_df$Gene_Name[gene_log])+
            guides(
              color = guide_legend(nrow = 2, byrow = T),
              fill = guide_legend(nrow = 2)
            )+
            NULL
        }

        # for the text below it, we put stats for each omics type
        text_rna <- paste0("<b>RNA Parameters:</b><br><b>- Best Model:</b> ",final_df$Best_Model_RNA[gene_log],"<br>")
        text_pro <- paste0("<b>Protein Parameters:</b><br><b>- Best Model:</b> ",final_df$Best_Model_Protein[gene_log],"<br>")

        # add each parameter for each model type
        if (!final_df$Best_Model_RNA[gene_log] %in% c("ECHO Joint", "ECHO Linear Joint")){
          for (p in 1:length(param_map[[final_df$Best_Model_RNA[gene_log]]])){
            text_rna <- paste0(text_rna, "<b>- ", param_map[[final_df$Best_Model_RNA[gene_log]]][p], ":</b> ",
                               final_df[gene_log, paste0(param_map[[final_df$Best_Model_RNA[gene_log]]][p], "_RNA")], "<br>")
          }
          for (p in 1:length(param_map[[final_df$Best_Model_Protein[gene_log]]])){
            text_pro <- paste0(text_pro, "<b>- ", param_map[[final_df$Best_Model_Protein[gene_log]]][p], ":</b> ",
                               final_df[gene_log, paste0(param_map[[final_df$Best_Model_Protein[gene_log]]][p], "_Protein")], "<br>")
          }
        } else {
          for (p in 1:length(param_map[[final_df$Best_Model_RNA[gene_log]]])){
            if (grepl("RNA", param_map[[final_df$Best_Model_RNA[gene_log]]][p])){
              text_rna <- paste0(text_rna, "<b>- ", param_map[[final_df$Best_Model_RNA[gene_log]]][p], ":</b> ",
                                 final_df[gene_log, param_map[[final_df$Best_Model_RNA[gene_log]]][p]], "<br>")
            } else {
              text_pro <- paste0(text_pro, "<b>- ", param_map[[final_df$Best_Model_Protein[gene_log]]][p], ":</b> ",
                                 final_df[gene_log, param_map[[final_df$Best_Model_Protein[gene_log]]][p]], "<br>")
            }
          }
        }

        # add the p_values, for the fit of each omics type
        text_rna <- paste0(text_rna,
          "<b>- P-Value:</b> ", final_df$P_Value_RNA[gene_log], "<br>",
          "<b>- BH Adj P-Value:</b> ", final_df$BH_Adj_P_Value_RNA[gene_log], "<br>"
        )
        text_pro <- paste0(text_pro,
          "<b>- P-Value:</b> ", final_df$P_Value_Protein[gene_log], "<br>",
          "<b>- BH Adj P-Value:</b> ", final_df$BH_Adj_P_Value_Protein[gene_log], "<br>"
        )

        # add the joint p-value, if joint
        if (final_df$Best_Model_RNA[gene_log] %in% c("ECHO Joint", "ECHO Linear Joint")){
          text_rna <- paste0(text_rna,
            "<b>- Joint P-Value:</b> ", final_df$P_Value_Joint[gene_log], "<br>",
            "<b>- BH Adj Joint P-Value:</b> ", final_df$BH_Adj_P_Value_Joint[gene_log], "<br>"
          )
          text_pro <- paste0(text_pro,
            "<b>- Joint P-Value:</b> ", final_df$P_Value_Joint[gene_log], "<br>",
            "<b>- BH Adj Joint P-Value:</b> ", final_df$BH_Adj_P_Value_Joint[gene_log], "<br>"
          )
        }

        # add the p-value of the slope, if it's a linear model
        if (final_df$Best_Model_RNA[gene_log] == "Linear"){
          text_rna <- paste0(
            text_rna,
            "<b>- Linear Slope P-Value:</b> ", final_df$P_Value_Linear_Slope_RNA[gene_log], "<br>",
            "<b>- BH Adj Linear Slope P-Value:</b> ", final_df$BH_Adj_P_Value_Linear_Slope_RNA[gene_log], "<br>"
          )
        }
        if (final_df$Best_Model_Protein[gene_log] == "Linear"){
          text_pro <- paste0(
            text_pro,
            "<b>- Linear Slope P-Value:</b> ", final_df$P_Value_Linear_Slope_Protein[gene_log], "<br>",
            "<b>- BH Adj Linear Slope P-Value:</b> ", final_df$BH_Adj_P_Value_Linear_Slope_Protein[gene_log], "<br>"
          )
        }

        # put the text together, depending on the omics types specified in the plot
        mid_text <-
          if (grepl("Both", input$om_type)){
            paste0(text_rna,"<br>",text_pro)
          } else if (input$om_type == "Only RNA") {
            text_rna
          } else {
            text_pro
          }

        # add the gene name to the text
        fin_text <- HTML(paste0("<b><u>", final_df$Gene_Name[gene_log], ":</b></u><br><br>", mid_text))
      }
    } else if (input$viz == "Heat Map" | input$viz == "Parameter Density Graph (PDG)") {
      # Heat Map ----

      # subset the heatmap
      rna_sub <- dat_rna[rna_sub_list[[input$subset_look]] & criteria_start_rna & criteria_end_rna,]
      pro_sub <- dat_pro[pro_sub_list[[input$subset_look]] & criteria_start_pro & criteria_end_pro,]

      # param subset
      fin_sub_rna <- final_df[rna_sub_list[[input$subset_look]] & criteria_start_rna & criteria_end_rna,]
      fin_sub_pro <- final_df[pro_sub_list[[input$subset_look]] & criteria_start_pro & criteria_end_pro,]

      if (input$viz == "Heat Map"){
        # get the gene order for each omics type
        ord_rna <- order_heat_map(fin_sub_rna, "RNA")
        ord_pro <- order_heat_map(fin_sub_pro, "Protein")

        # now make the heat map
        rna_map <- make_hm_matrix(rna_sub, ord_rna, num_reps, input$heat_subset_rep, timen)
        pro_map <- make_hm_matrix(pro_sub, ord_pro, num_reps, input$heat_subset_rep, timen)

        # you can't have an overlaid plot for heat maps -- just side by side
        if (grepl("Both", input$om_type)){
          rna_g <- make_hm_plot(rna_map, "RNA", input$subset_look)
          pro_g <- make_hm_plot(pro_map, "Protein", input$subset_look)

          # arrange and pass to the final plot
          fin_plot <- as.ggplot(grid.arrange(rna_g, pro_g, ncol = 2))
        } else if (input$om_type == "Only RNA"){ # it's an rna plot
          rna_g <- make_hm_plot(rna_map, "RNA", input$subset_look)

          # arrange and pass to the final plot
          fin_plot <- rna_g
        } else { # it's a protein plot
          pro_g <- make_hm_plot(pro_map, "Protein", input$subset_look)

          # arrange and pass to the final plot
          fin_plot <- as.ggplot(grid.arrange(rna_g, pro_g, ncol = 2))
        }

        # for the text below it, we put total genes
        text_rna <- paste0("<b>Total Genes, RNA:</b> ",nrow(fin_sub_rna),"<br>")
        text_pro <- paste0("<b>Total Genes, Protein:</b> ",nrow(fin_sub_pro),"<br>")

        # text to display based on the omics type
        mid_text <-
          if (grepl("Both", input$om_type)){
            paste0(text_rna,"<br>",text_pro)
          } else if (input$om_type == "Only RNA") {
            text_rna
          } else {
            text_pro
          }

        # output the final text
        fin_text <- HTML(paste0(mid_text))
      }

      # Parameter Density Graph (PDG) ----
      if (input$viz == "Parameter Density Graph (PDG)"){
        # param subset
        coeff_sub_rna <- final_df[rna_sub_list[[input$subset_look]] & criteria_start_rna & criteria_end_rna,paste0(input$coeff, "_RNA")]
        coeff_sub_pro <- final_df[pro_sub_list[[input$subset_look]] & criteria_start_pro & criteria_end_pro,paste0(input$coeff, "_Protein")]

        # remove all NAs
        coeff_sub_rna <- coeff_sub_rna[!is.na(coeff_sub_rna)]
        coeff_sub_pro <- coeff_sub_pro[!is.na(coeff_sub_pro)]

        # the data frame with the coefficients
        gg.df <- data.frame(
          "dens" = c(coeff_sub_rna,coeff_sub_pro),
          "label" = c(rep("RNA", length(coeff_sub_rna)), rep("Protein",length(coeff_sub_pro))),
          stringsAsFactors = F
        )

        # options to make the plot look pretty
        opts <- list(
          geom_density(alpha = .3),
          theme_bw(),
          scale_fill_manual("",values = c("RNA" = rna_blue_high, "Protein" = pro_red_high)),
          scale_color_manual("",values = c("RNA" = rna_blue_high, "Protein" = pro_red_high)),
          scale_y_continuous(expand = expand_scale(mult = c(0,.05))),
          scale_x_continuous(expand = expand_scale(mult = c(0,0))),
          theme(
            text = element_text(size = 20),
            legend.position = "bottom",
            legend.direction = "horizontal",
            plot.title = element_text(hjust = .5)
          ),
          xlab(input$coeff),
          ylab("Density")
        )

        # how to display the plot based on which omics type to view
        if (input$om_type == "Both, Overlaid"){
          # plot densities of both omics types on the same plot
          fin_plot <-
            if (length(coeff_sub_rna) > 1 & length(coeff_sub_pro) > 1){
              ggplot(gg.df, aes(dens, fill = label, color = label))+
                ggtitle(paste0("RNA and Protein: ", input$coeff, ", ", input$subset_look, " Subset"))+
                opts
            } else {
              ggplot()+theme_bw()
            }
        } else {
          # plot rna and protein densities separately
          p_rna <-
            if (length(coeff_sub_rna) > 1){
              ggplot(gg.df[gg.df$label == "RNA",], aes(dens, fill = label, color = label))+
                ggtitle(paste0("RNA: ", input$coeff, ", ", input$subset_look, " Subset"))
            } else {
              ggplot()+theme_bw()
            }
          p_pro <-
            if (length(coeff_sub_pro) > 1){
              ggplot(gg.df[gg.df$label == "Protein",], aes(dens, fill = label, color = label))+
                ggtitle(paste0("Protein: ", input$coeff, ", ", input$subset_look, " Subset"))
            } else {
              ggplot()+theme_bw()
            }

          # arrange final plot based on what omics type to view
          fin_plot <-
            if (input$om_type == "Both, Side-by-Side"){
              as.ggplot(grid.arrange(p_rna + opts, p_pro + opts, ncol = 2))
            } else if (input$om_type == "Only RNA"){
              p_rna + opts
            } else if (input$om_type == "Only Protein"){
              p_pro + opts
            }
        }
        # for the text below it, we put stats based on the quantile info
        text_rna <- paste0("<b>RNA Distribution (Total Genes: ",length(coeff_sub_rna),"):</b><br>")
        text_pro <- paste0("<b>Protein Distribution (Total Genes: ",length(coeff_sub_pro),"):</b><br>")

        # if there are any parameters for rna, we can get the summary
        if (length(coeff_sub_rna) > 0){
          sumy_rna <- summary(coeff_sub_rna)

          text_rna <- paste0(
            text_rna,
            "<b>- Min: </b>", sumy_rna[1],"<br>",
            "<b>- 1st Quantile: </b>", sumy_rna[2],"<br>",
            "<b>- Median: </b>", sumy_rna[3],"<br>",
            "<b>- Mean: </b>", sumy_rna[4],"<br>",
            "<b>- 3rd Quantile: </b>", sumy_rna[5],"<br>",
            "<b>- Max: </b>", sumy_rna[6],"<br>"
          )
        }
        # if there are any parameters for protein, we can get the summary
        if (length(coeff_sub_pro) > 0){
          sumy_pro <- summary(coeff_sub_pro)

          text_pro <- paste0(
            text_pro,
            "-<b>Min: </b>", sumy_pro[1],"<br>",
            "-<b>1st Quantile: </b>", sumy_pro[2],"<br>",
            "-<b>Median: </b>", sumy_pro[3],"<br>",
            "-<b>Mean: </b>", sumy_pro[4],"<br>",
            "-<b>3rd Quantile: </b>", sumy_pro[5],"<br>",
            "-<b>Max: </b>", sumy_pro[6],"<br>"
          )
        }

        # arrange text based on omics type plotted
        mid_text <-
          if (grepl("Both", input$om_type)){
            paste0(text_rna,"<br>",text_pro)
          } else if (input$om_type == "Only RNA") {
            text_rna
          } else {
            text_pro
          }

        # add title to final text
        fin_text <- HTML(paste0("<b><u>", input$coeff, ":</b></u><br><br>", mid_text))

      }

    } else if (input$viz == "Heat Map Comparison"){
      # Heat Map Comparison ----

      # start by getting the names for rna and protein
      rna_focus <- strsplit(input$rna_focus,"\n", fixed=T)[[1]]
      pro_focus <- strsplit(input$pro_focus,"\n", fixed=T)[[1]]

      # subset the heatmaps and data accordingly:
      # subset the heatmap
      rna_sub <- dat_rna[final_df$Gene_Name %in% rna_focus,]
      pro_sub <- dat_pro[final_df$Gene_Name %in% pro_focus,]

      # param subset
      fin_sub_rna <- final_df[final_df$Gene_Name %in% rna_focus,]
      fin_sub_pro <- final_df[final_df$Gene_Name %in% pro_focus,]

      # add empty lines for things appearing in the other
      in_both <- intersect(rna_focus, pro_focus)
      # rna additions
      if (sum(!pro_focus %in% in_both) > 0){
        rna_sub[(nrow(rna_sub)+1):(nrow(rna_sub)+sum(!pro_focus %in% in_both)),] <- NA
        fin_sub_rna[(nrow(fin_sub_rna)+1):(nrow(fin_sub_rna)+sum(!pro_focus %in% in_both)),] <- NA
        fin_sub_rna[(length(rna_focus)+1):(length(rna_focus)+sum(!pro_focus %in% in_both)),"Gene_Name"] <- pro_focus[!pro_focus %in% in_both]
      }
      # pro additions
      if (sum(!rna_focus %in% in_both) > 0){
        pro_sub[(nrow(pro_sub)+1):(nrow(pro_sub)+sum(!rna_focus %in% in_both)),] <- NA
        fin_sub_pro[(nrow(fin_sub_pro)+1):(nrow(fin_sub_pro)+sum(!rna_focus %in% in_both)),] <- NA
        fin_sub_pro[(length(pro_focus)+1):(length(pro_focus)+sum(!rna_focus %in% in_both)),"Gene_Name"] <- rna_focus[!rna_focus %in% in_both]
      }

      # now ordering based on the "dominant" comparison
      ord <-
        if (input$focus_ord == "RNA"){
          order_heat_map_comparison(fin_sub_rna, "RNA")
        } else {
          order_heat_map_comparison(fin_sub_pro, "Protein")
        }
      ord_rna <- rev(match(ord, fin_sub_rna$Gene_Name))
      ord_pro <- rev(match(ord, fin_sub_pro$Gene_Name))

      # now make the heat map
      rna_map <- make_hm_matrix(rna_sub, ord_rna, num_reps, input$heat_subset_rep, timen)
      pro_map <- make_hm_matrix(pro_sub, ord_pro, num_reps, input$heat_subset_rep, timen)

      # this always shows both side by side
      rna_g <- make_hm_plot(rna_map, "RNA", "Comparison", T, fin_sub_rna$Gene_Name[ord_rna])
      pro_g <- make_hm_plot(pro_map, "Protein", "Comparison", T, fin_sub_pro$Gene_Name[ord_pro])

      # arrange and pass to the final plot
      fin_plot <- as.ggplot(grid.arrange(rna_g, pro_g, ncol = 2))

      # text below is a list of genes in a venn diagram of the subsets
      fin_text <- HTML(
        paste0(
          "<b>Genes Appearing in Both: </b>", paste(in_both, collapse = " "), "<br>",
          "<b>Genes Appearing in RNA Only: </b>", paste(rna_focus[!rna_focus %in% in_both], collapse = " "), "<br>",
          "<b>Genes Appearing in Protein Only: </b>", paste(pro_focus[!pro_focus %in% in_both], collapse = " ")
        )
      )
    }

    # change the plot size, if specified as a number
    output$plot_ui <- renderUI({
      # width
      if (input$viz_wid == "" | grepl("\\D", input$viz_wid)){
        wid <- "100%"
      } else {
        wid <- paste0(input$viz_wid,"px")
      }

      # height
      if (input$viz_height == "" | grepl("\\D", input$viz_height)){
        height <- "800px"
      } else {
        height <- paste0(input$viz_height,"px")
      }

      plotOutput("plot_viz", width = wid, height = height)
    })

    # render final visualization plot
    output$plot_viz <- renderPlot({
      fin_plot
    })

    # render final visualization text
    output$viz_text <- renderUI({
      fin_text
    })

    # function to download png of plot
    output$downloadPlot <- downloadHandler(
      filename = function() { paste0(input$viz,"_",input$om_type, '_MOSAIC.png') },
      content = function(file) {
        if (input$viz_wid == "" | grepl("\\D", input$viz_wid)){
          wid <- 800
        } else {
          wid <- as.numeric(input$viz_wid)
        }

        if (input$viz_height == "" | grepl("\\D", input$viz_height)){
          height <- 800
        } else {
          height <- as.numeric(input$viz_height)
        }

        png(file, width = wid, height = height)
        print(fin_plot)
        dev.off()
      },
      contentType='image/png'
    )
  })

  # close server ----
}

# run app ----

shinyApp(ui, server)
