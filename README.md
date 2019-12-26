# Iterated Multi-Source Exchangeability Models (iMEMs)
This repository containst the data and code used to implement all simulation studies and real data anlayses for the manuscript entitled "Iterated Multi-Source Exchangeability Models for Individualized Inference with an Application to Mobile Sensor Data".  Briefly, Iterated Multi-source Exchangeability Models (iMEMs) are a method for efficiently estimating high-dimensional MEMs, which are a flexible Bayesian method for incorporating supplementary data into the analysis of a primary data source.  Code is licensed under the MIT License, and data under the CC-0 License.  

## "Code" Folder
This folder contains 7 different R scripts used to implement iMEMs and generate the tables and figures in the manuscript.  In each script (excepting **iMEM_functions.R**), a *root* variable is defined at the beginning of the script.  This variable indicates the location of the repository on the local machine.  Once this variable has been changed appropriately, all code should run properly.  Contents of the Code folder are as follows:

  - **iMEM_functions.R**:  This script contains a set of functions for implementing iMEMs and executing the simulation for use in the other scripts.  This script is sourced from all the other R scripts, so should never need to be separately run.
  - **Similarity_Score_Heatmap.R**: This script generates the similarity score heatmap, Supplementary Figure 1 in the manuscript. 
  - **Asymptotics_Plot.R**: This script generates the asymptotic weights plot, Supplementary Figure 2 in the manuscript.
  - **Primary_Simulation_Code.R**: This script implements the simulation study described in Section 5 of the iMEM manuscript.  Simulation computation time is on the order of weeks, so simulation results are stored in the "Data/Simulation Results" subfolder.  
  - **Secondary_Simulation_Code.R**:  This script implements the secondary simulation study described in Appendix B of the supplementary materials.  Simulation results are stored in the "Data/Simulation Results" subfolder.  
  - **Simulation_Plots.R**: This script imports simulation results generates all simulation-related plots in the manuscript and supplementary materials: manuscript Figures 3 and 4, and Supplementary Figures 3-7, along with other plots from the simulation.  Imports data from the "Data/Simulation Results" subfolder. Text search "Manuscript Figure" or "Supplementary Figure" for the code snippets specifically corresponding to each figure.  
  - **Real_Data_Analysis.R**: This script implements the real data analysis described in Section 5 of the manuscript, and generates manuscript Figure 5 and Supplementary Figure 8.  Text search "Manuscript Figure" or "Supplementary Figure" for the code snippets specifically corresponding to each figure.  Imports **imem_realdata.Rdata** from the "Data" folder. 

## "Data" Folder
This folder contains the data imported in the scripts **Simulation_Plots.R** and **Real_Data_Analysis.R**.  Contents are as follows:

 - **imem_realdata.Rdata**:  Real trip and activity data with accompanying self-reported emotional state status.  Collected using the Daynamica mobile application in a Minneapolis-area study.  Variable names with descriptions are below.
 
    | Variable Name | Description | Variable Type |
    | ------ | ------ | ------ | 
    | surveycode | User ID | character string |
    | type | Indicator for activity or trip | character string |
    | primary_mode | Trip mode | character string |
    | happy | self-reported happiness intensity, 1-7 scale with 7 being highest | Integer: 1-7 |
    | tired | self-reported tiredness intensity, 1-7 scale with 7 being highest | Integer: 1-7 |
    | stressful | self-reported stress intensity, 1-7 scale with 7 being highest | Integer: 1-7 |
    | sad | self-reported sadness intensity, 1-7 scale with 7 being highest | Integer: 1-7 |
    | meaningful | self-reported meaningful intensity, 1-7 scale with 7 being highest | Integer: 1-7 |  
    | pain | self-reported pain intensity, 1-7 scale with 7 being highest |   Integer: 1-7 |

 - **"Data/Simulation Results"**: Subfolder containing .Rdata files with the simulation results from **Primary_Simulation_Code.R** and **Secondary_Simulation_Code.R**.  Data are stored as multi-dimensional arrays which are processed and used to generate plots in **Simulation_Plots.R**.  
