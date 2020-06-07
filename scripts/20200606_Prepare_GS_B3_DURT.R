library(openCyto) # 1.24.0
library(CytoML) # 1.12.0
library(flowCore) # required for description()
library(flowWorkspace) # required for gh_pop_get_data()
library(ggcyto) # devtools::install_github("RGLab/ggcyto", ref="ggplot3.3") for update_theme()
library(here)
library(tidyverse)
library(readxl)

# Read in Batch 3 workspace and prepare GatingSet. (Second version of workspace)

date <- 20200606

xml_path_b3 <- here::here("data/20200604 HAARVI DURT B3/20200605 HAARVI DURT B3_KY.xml")
fcs_subfolder_b3 <- here::here("data/20200604 HAARVI DURT B3")

ws_b3 <- open_flowjo_xml(xml_path_b3)
merge(fj_ws_get_sample_groups(ws_b3), fj_ws_get_samples(ws_b3), by = "sampleID")
fj_ws_get_keywords(ws_b3, 37)
names(fj_ws_get_keywords(ws_b3, 37))
keywords2import <- c("TUBE NAME", "EXPERIMENT NAME", "$DATE", "SAMPLE ID")

sampleGroup <- "Samples"
gs_b3 <- flowjo_to_gatingset(ws_b3, name=sampleGroup, keywords=keywords2import,
                                  path=fcs_subfolder_b3)
pop_lists <- lapply(gs_b3, gh_get_pop_paths)
# Check that there is one gating tree for all samples
unique(pop_lists)

pData(gs_b3)$filename <- sapply(rownames(pData(gs_b3)), function(x) {
  gsub("/.+/", "", description(gh_pop_get_data(gs_b3[[x]]))$FILENAME)
}, USE.NAMES = F)
pData(gs_b3)$rowname <- rownames(pData(gs_b3))
pData(gs_b3)

# Read in the patient manifest
manifest <- readxl::read_excel(here::here("data/Seshadri manifest 22May2020.xlsx"), range = "A6:T80") %>% 
  dplyr::rename("Comments" = `...20`)
# Add metadata to pData
pData_tmp <- pData(gs_b3) %>% 
  mutate(`SAMPLE ID` = toupper(`SAMPLE ID`)) %>% 
  left_join(manifest %>%
              dplyr::select(`Record ID`, `Sample ID`, Cohort, Age, Sex, Race, `Hispanic?`, `Days symptom onset to visit 1`, `Pair ID`, Comments),
            by = c("SAMPLE ID" = "Record ID")) %>% 
  mutate(Batch = 3)
rownames(pData_tmp) <- rownames(pData(gs_b3))
pData(gs_b3) <- pData_tmp

pData(parameters(gh_pop_get_data(gs_b3[[1]])))[,c(1, 2)]

# Add names to all channels
dput(unname(pData(parameters(gh_pop_get_data(gs_b3[[1]])))[,2]))
markernames_b3<- c("Time", "FSC-A", "FSC-H", "SSC-A", "SSC-H", "CD8a", "L/D", "GD", "CD56", 
                   "CD3", "aGC", "CD4", "Delta2", "R660-A", "CD19", "CCR7", "CD14", "CD38", 
                   "V570-A", "V510-A", "MR1", "CD45RA", "HLADR")
names(markernames_b3) <- pData(parameters(gh_pop_get_data(gs_b3[[1]])))[,1]
# markernames_b3
markernames(gs_b3) <- markernames_b3
pData(parameters(gh_pop_get_data(gs_b3[[1]])))[,c(1, 2)]
# name   desc
# $P1      Time   Time
# $P2     FSC-A  FSC-A
# $P3     FSC-H  FSC-H
# $P4     SSC-A  SSC-A
# $P5     SSC-H  SSC-H
# $P6  <B710-A>   CD8a
# $P7  <B515-A>    L/D
# $P8  <G780-A>     GD
# $P9  <G660-A>   CD56
# $P10 <G610-A>    CD3
# $P11 <G575-A>    aGC
# $P12 <R780-A>    CD4
# $P13 <R710-A> Delta2
# $P14   R660-A R660-A
# $P15 <V780-A>   CD19
# $P16 <V710-A>   CCR7
# $P17 <V655-A>   CD14
# $P18 <V610-A>   CD38
# $P19   V570-A V570-A
# $P20   V510-A V510-A
# $P21 <V450-A>    MR1
# $P22 <U730-A> CD45RA
# $P23 <U395-A>  HLADR

png(here::here(sprintf("out/QC/B3_GatingTree_%s.png", date)), width = 7, height = 5, units = "in", res = 300)
(plot(gs_b3, fontsize=15, bool=T))
dev.off()

save_gs(gs_b3, here::here("out/GatingSets/20200606_HAARVI_DURT_GatingSet_B3"))

#####################################################################

# gs_b3 <- load_gs(here::here("out/GatingSets/20200606_HAARVI_DURT_GatingSet_B3"))

dput(gh_get_pop_paths(gs_b3))

# Perform QC on CD3 counts
cd3_path <- "/Time/S/Live/3+"
cd3_counts_b3 <- pData(gs_b3) %>% 
  left_join(gs_pop_get_count_fast(gs_b3, subpopulations = c(cd3_path)) %>% 
              pivot_wider(id_cols = name, names_from = "Population", values_from = "Count") %>% 
              rename(CD3 = !!cd3_path),
            by = c("rowname" = "name")) %>% 
  dplyr::select("SAMPLE ID", "Cohort", "Age", "Sex", "Race", "Days symptom onset to visit 1", "Batch", "CD3")

png(here::here(sprintf("out/QC/B3_CD3_Counts_%s.png", date)),
    width = 10, height = 6, units="in", res=300)
ggplot(cd3_counts_b3 %>% 
         mutate(Color = ifelse(CD3 < 10000, "red", "black")),
       aes(`SAMPLE ID`, CD3, color = Color)) +
  geom_point() +
  geom_hline(aes(yintercept = 10000, color="red"), linetype="dashed") +
  scale_color_identity() +
  ggtitle("Batch 3 CD3 Counts")
dev.off()

#####################################################################

# Plot gates using ggcyto

plotter <- function(myGS, myGates) {
  firstGate <- myGates[[1]]
  currentGateBoundaries <- attributes(gh_pop_get_gate(myGS[[1]], firstGate))$boundaries
  currentXaxis <- colnames(currentGateBoundaries)[[1]]
  currentYaxis <- colnames(currentGateBoundaries)[[2]]
  parentGate <- sub("\\/[^\\/]+$", "", firstGate)
  
  ggcyto(myGS,
         aes(!!currentXaxis, !!currentYaxis),
         subset = if(parentGate == "") "root" else parentGate,
         filter = marginalFilter) +
    geom_hex(bins=128) +
    geom_gate(myGates) +
    axis_x_inverse_trans() + axis_y_inverse_trans() +
    ggcyto_par_set(limits = "instrument") +
    facet_wrap(. ~ `SAMPLE ID`) +
    theme_bw(base_size=20) +
    theme(#plot.title = element_blank(),
      legend.position = "none",
      strip.text.x = element_text(size = 14, margin = margin(0,0,0,0, "cm")),
      panel.grid.major = ggplot2::element_blank()) +
    theme_bw(base_size=28) +
    geom_stats(size=8, alpha=0.4) +
    labs(title=myGates[[1]], caption = "Batch 3")
}

gates2draw <- list("/Time", "/Time/S", "/Time/S/Live", "/Time/S/Live/3+",
                   c("/Time/S/Live/3+/4+", "/Time/S/Live/3+/8+"), 
                   c("/Time/S/Live/3+/CD38+", "/Time/S/Live/3+/HLADR+"), 
                   "/Time/S/Live/3+/GD-",
                   c("/Time/S/Live/14+", "/Time/S/Live/19+"),
                   "/Time/S/Live/56+",
                   "/Time/S/Live/iNKT", 
                   "/Time/S/Live/MAIT")

for(currentGates in gates2draw) {
  png(filename = file.path(here::here("out/QC/FACS_Plots"),
                           sprintf("%s_B3_%s.png",
                                   sub(".*\\/([^\\/]+$)", "\\1", currentGates[[1]]),
                                   format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
      width = 1900, height = 1030)
  print(plotter(gs_b3, currentGates))
  dev.off()
}

# Memory gates and Vd2 were drawn on histograms and won't work with the above plotter function
plotter_hist <- function(myGS, myGates, yAxis="GD") {
  firstGate <- myGates[[1]]
  firstGateAttributes <- attributes(gh_pop_get_gate(myGS[[1]], firstGate))
  currentXaxis <- names(firstGateAttributes$min)
  
  currentYaxis <- if(length(myGates) == 2) {
    secondGate <- myGates[[2]]
    secondGateAttributes <- attributes(gh_pop_get_gate(myGS[[1]], secondGate))
    names(secondGateAttributes$min)
  } else {
    # Defaults to GD. Works for now.
    "GD"
  }
  
  parentGate <- sub("\\/[^\\/]+$", "", firstGate)
  
  ggcyto(myGS,
         aes(!!currentXaxis, !!currentYaxis),
         subset = if(parentGate == "") "root" else parentGate,
         filter = marginalFilter) +
    geom_hex(bins=128) +
    geom_gate(myGates) +
    axis_x_inverse_trans() + axis_y_inverse_trans() +
    ggcyto_par_set(limits = "instrument") +
    facet_wrap(. ~ `SAMPLE ID`) +
    theme_grey(base_size = 20) +
    theme(#plot.title = element_blank(),
      legend.position = "none",
      strip.text.x = element_text(size = 14, margin = margin(0,0,0,0, "cm")),
      panel.grid.major = ggplot2::element_blank()) +
    theme_bw(base_size=28) +
    geom_stats(size=8, alpha=0.4) +
    labs(title=myGates[[1]], caption = "Batch 3")
}

cd4_mem_gates <- c("/Time/S/Live/3+/4+/CCR7+", "/Time/S/Live/3+/4+/CD45RA+")
cd8_mem_gates <- c("/Time/S/Live/3+/8+/CCR7+", "/Time/S/Live/3+/8+/CD45RA+")

for(currentGates in list(cd4_mem_gates, cd8_mem_gates, "/Time/S/Live/3+/GD+/Vd2+")) {
  png(filename = file.path(here::here("out/QC/FACS_Plots"),
                           sprintf("%s_B3_%s.png",
                                   sub(".*\\/([^\\/]+$)", "\\1", currentGates[[1]]),
                                   format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
      width = 1900, height = 1030)
  print(plotter_hist(gs_b3, currentGates))
  dev.off()
}
