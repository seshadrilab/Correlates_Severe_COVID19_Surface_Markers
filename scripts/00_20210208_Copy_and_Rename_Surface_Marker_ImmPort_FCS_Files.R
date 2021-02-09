# This script copies the downloaded ImmPort FCS files to the appropriate folders for data analysis and renames them to the original FCS file names. 
# This is necessary in order to perform the rest of the downstream analyses.
# Before you run this script, download the result files for study SDY1680 from ImmPort (https://www.immport.org/home).
# As of 20210208, these FCS files are on ImmPort under `Browse Shared Data > SDY1680 > ResultFiles (569 files)` https://browser.immport.org/browser?path=SDY1680.
# This will download both the ICS and surface marker panel FCS files, including compensation files.
# The download process using Aspera will take a while (18.3 GB, took me 2 hrs and 42 min)

library(here)
library(purrr)

##################################
# MODIFY THIS NEXT PATH TO YOUR PERSONAL DOWNLOAD FOLDER
ImmPort_dl_path <- "/home/malisa/Downloads/ResultFiles"
##################################

surface_marker_file_map <- read.table(here::here("data/Surface_Markers_ImmPort_FCS_FileMap.tsv"), sep = "\t", header = T,
                                      colClasses = c("character", "character", "character", "character", "numeric", "numeric"))

# Loop through each FCS file required for the Surface Marker analysis and copy it to the analysis folder, and then rename it
quietly(pmap(surface_marker_file_map %>%
       dplyr::select("Original_File_Name", "ImmPort_File_Name",
                     "Destination_Folder_Path", "ImmPort_Subfolder"),
     function(Original_File_Name, ImmPort_File_Name, Destination_Folder_Path, ImmPort_Subfolder) {
  # Create the destination folder if it doesn't already exist
  if(!dir.exists(here::here(Destination_Folder_Path))) {
    cat(sprintf("Creating folder %s\n", here::here(Destination_Folder_Path)))
    dir.create(here::here(Destination_Folder_Path), recursive=T)
  }
  
  # Copy the FCS file to the destination folder
  cat(sprintf("Copying %s to %s\n", paste0(ImmPort_Subfolder, "/", ImmPort_File_Name), Destination_Folder_Path))
  file.copy(file.path(ImmPort_dl_path, ImmPort_Subfolder, ImmPort_File_Name),
            here::here(Destination_Folder_Path),
            overwrite = F)
  # Rename the copied file to the original FCS file name
  file.rename(here::here(Destination_Folder_Path, ImmPort_File_Name),
              here::here(Destination_Folder_Path, Original_File_Name))
}))
