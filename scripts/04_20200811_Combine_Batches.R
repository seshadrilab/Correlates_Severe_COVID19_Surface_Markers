library(openCyto) # 1.24.0
library(CytoML) # 1.12.0
library(flowCore) # required for description()
library(flowWorkspace) # required for gh_pop_get_data()
library(ggcyto) # devtools::install_github("RGLab/ggcyto", ref="ggplot3.3") for update_theme()
library(here)
library(tidyverse)

date <- 20200812

gs1 <- load_gs(here::here("out/GatingSets/20200811_HAARVI_DURT_GatingSet_B1"))
gs2 <- load_gs(here::here("out/GatingSets/20200811_HAARVI_DURT_GatingSet_B2"))
gs3 <- load_gs(here::here("out/GatingSets/20200811_HAARVI_DURT_GatingSet_B3"))

# Make sure nodes and markers are consistent between the three batches
setdiff(sort(gh_get_pop_paths(gs1)), sort(gh_get_pop_paths(gs2)))
setdiff(sort(gh_get_pop_paths(gs1)), sort(gh_get_pop_paths(gs3)))
all(sort(gh_get_pop_paths(gs1)) == sort(gh_get_pop_paths(gs2)))
all(sort(gh_get_pop_paths(gs1)) == sort(gh_get_pop_paths(gs3)))
all(colnames(pData(gs1)) == colnames(pData(gs2)))
all(colnames(pData(gs1)) == colnames(pData(gs3)))
all(markernames(gs1) == markernames(gs2))
all(markernames(gs1) == markernames(gs3))

pData(parameters(gh_pop_get_data(gs1[[1]])))[,c(1, 2)]

# # What about the transformation?
# gs1_trans_list <- gh_get_transformations(gs1[[1]])
# gs2_trans_list <- gh_get_transformations(gs2[[1]])
# gs3_trans_list <- gh_get_transformations(gs3[[1]])

# Combine Batches into a GatingSetList
gs <- GatingSetList(list(gs1, gs2, gs3))

# Clean up the Cohort column by changing any remaining "NA" or "Healthy control" to "Healthy"
pData(gs)$Cohort <- ifelse(pData(gs)$Cohort %in% c(NA, "Healthy control"), "Healthy", pData(gs)$Cohort)

# Lower the lower boundary of the /Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/iNKT gate for just patient 109C
ggcyto(subset(gs, `SAMPLE ID` == "109C"), aes(CD3, aGC), subset = "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+",
       filter=marginalFilter) + geom_hex(bins=128) +
  # axis_x_inverse_trans() + axis_y_inverse_trans() +
  ggcyto_par_set(limits = "instrument") +
  facet_wrap(. ~ `SAMPLE ID`) +
  theme_bw(base_size = 28) +
  geom_gate("/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/iNKT") +
  theme(
         legend.position = "none",
         strip.text.x = element_text(margin = margin(0,0,0,0, "cm")),
         panel.grid.major = ggplot2::element_blank()) +
  geom_stats(size=8, alpha=0.4) +
  geom_hline(aes(yintercept = 1650))
inkt_gate_109C <- gh_pop_get_gate(subset(gs, `SAMPLE ID` == "109C")[[1]], "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/iNKT")
attributes(inkt_gate_109C)$boundaries[c(1,4), "<G575-A>"] <- 1650
gh_pop_set_gate(subset(gs, `SAMPLE ID` == "109C")[[1]], "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/iNKT", inkt_gate_109C)
recompute(subset(gs, `SAMPLE ID` == "109C"), "/Time/S/Live")
ggcyto(subset(gs, `SAMPLE ID` == "109C"), aes(CD3, aGC), subset = "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+",
       filter=marginalFilter) + geom_hex(bins=128) +
  # axis_x_inverse_trans() + axis_y_inverse_trans() +
  ggcyto_par_set(limits = "instrument") +
  facet_wrap(. ~ `SAMPLE ID`) +
  theme_bw(base_size = 28) +
  geom_gate("/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/iNKT") +
  theme(
    legend.position = "none",
    strip.text.x = element_text(margin = margin(0,0,0,0, "cm")),
    panel.grid.major = ggplot2::element_blank()) +
  geom_stats(size=8, alpha=0.4)

###########################################################

# Extract cell count data to take a look at batch effect. Make boxplots.
cd3_path <- "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+"
cd4_path <- "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/4+"
qc_counts <- gs_pop_get_count_with_meta(gs, subpopulations = c(cd3_path, cd4_path)) %>% 
  dplyr::select(Population, Count, "SAMPLE ID", "Cohort", "Age", "Sex", "Race", "Hispanic?", 
                "Days symptom onset to visit 1", "Batch") %>% 
  pivot_wider(names_from = Population, values_from = Count) %>% 
  dplyr::rename(CD3_Count = !!cd3_path,
                CD4_Count = !!cd4_path)

unique(table(qc_counts$`SAMPLE ID`)) # each patient appears once

figwidth <- 5
figheight <- 3.5

cd3_counts_vs_batch_plot <- ggplot(qc_counts, aes(x=Batch, y=CD3_Count, group=Batch)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.15, height = 0) +
  theme_bw(base_size = 22) +
  theme(axis.title = element_text(size=22),
        axis.text = element_text(color="black", size=18),
        plot.title = element_blank(),
        legend.position = "none")
ggsave(filename = here::here("out/QC/QC_Counts/AllBatches_CD3_Counts_vs_Batch.png"),
       plot = cd3_counts_vs_batch_plot,
       width = figwidth,
       height = figheight,
       units = "in")

cd4_percent_cd3_vs_batch_plot <- ggplot(qc_counts %>% 
         mutate(CD4_prop_CD3 = CD4_Count / CD3_Count),
       aes(x=Batch, y=CD4_prop_CD3, group=Batch)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.15, height = 0) +
  theme_bw(base_size = 22) +
  theme(axis.title = element_text(size=22),
        axis.text = element_text(color="black", size=18),
        plot.title = element_blank(),
        legend.position = "none") +
  scale_y_continuous(labels = function(x) paste0(x*100)) +
  labs(y = "% CD4 of CD3 cells")
ggsave(filename = here::here("out/QC/QC_Counts/AllBatches_CD4_Percent_CD3_vs_Batch.png"),
       plot = cd4_percent_cd3_vs_batch_plot,
       width = figwidth,
       height = figheight,
       units = "in")

##########################################

# Now add in boolean gates as necessary to define cell subpopulations of interest

# Cell Populations (DURT panel)
# 
# CD14-CD19-CD3+ as % parent (T cells)
# CD56-CD19+CD14-CD3- as % parent (B cells)
# CD56-CD19-CD14+CD3- as % parent (Monocytes)
# CD56+CD19-CD14-CD3- as % parent (NK cells)
# 
# Monocyte:Lymphocytes (CD14+ : CD3+ U CD19+)
# 
# CD4 as % CD3
# CD8 as % CD3
# 
# Activated CD4 (HLA-DR+CD38+) as % CD4
# Memory CD4 (Naive, TCM, TEM, TEMRA) as % CD4
# 
# Activated CD8 (HLA-DR+CD38+) as % CD8
# Memory CD8 (Naive, TCM, TEM, TEMRA) as % CD8
# 
# DURTs
# 
# Total GD as % CD3
# Vd2+GD+ as % CD3
# Vd2-GD+ as % CD3
# 
# Activated GD (HLA-DR+CD38+) as % GD
# Activated Vd2+GD+ (HLA-DR+CD38+) as % Vd2+GD+
# Activated Vd2-GD+ (HLA-DR+CD38+) as % Vd2-GD+
# 
# Total MR1 as % CD3
# Activated MR1 (HLA-DR+CD38+) as % MR1
# 
# Total iNKT as % CD3
# Activated iNKT (HLA-DR+CD38+) as % iNKT

dput(gh_get_pop_paths(gs[[1]]))

# Re-define PBMC subpopulations to be mutually exclusive where necessary:
# The latest gating strategy should already take care of this for all but the NK cells, which should be defined as CD3-
# And remember some T cells are CD56+ (NK-like T cells, often CD1d-specific)
live_path <- "/Time/S/Live"
cd3_path <- "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+" # T cells
cd14_path <- "/Time/S/Live/14+" # Monocytes
cd19_path <- "/Time/S/Live/19+" # B cells
cd56_path <- "/Time/S/Live/CD14-CD19-/LD/Singlet/56+" # NK cells and NK-like T cells

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", cd56_path,
                                                         "&!", cd3_path))))),
           parent = live_path, name = "NK_cells")

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", cd3_path,
                                                         "|", cd19_path))))),
           parent = live_path, name = "Lymphocytes")

# Add "Not Vd2" gate under GD
gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol("!/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/GD+/Vd2+")))),
           parent = "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/GD+", name = "Vd2-")

# Add memory gates to CD4 and CD8
cd4_ccr7_path <- "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/4+/CCR7+"
cd4_cd45ra_path <- "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/4+/CD45RA+"
cd8_ccr7_path <- "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/8+/CCR7+"
cd8_cd45ra_path <- "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/8+/CD45RA+"

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", cd4_ccr7_path,
                                                         "&", cd4_cd45ra_path))))),
           parent = "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/4+", name = "Naive")
gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", cd4_ccr7_path,
                                                         "&!", cd4_cd45ra_path))))),
           parent = "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/4+", name = "TCM")
gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("!", cd4_ccr7_path,
                                                         "&", cd4_cd45ra_path))))),
           parent = "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/4+", name = "TEMRA")
gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("!", cd4_ccr7_path,
                                                         "&!", cd4_cd45ra_path))))),
           parent = "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/4+", name = "TEM")

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", cd8_ccr7_path,
                                                         "&", cd8_cd45ra_path))))),
           parent = "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/8+", name = "Naive")
gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", cd8_ccr7_path,
                                                         "&!", cd8_cd45ra_path))))),
           parent = "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/8+", name = "TCM")
gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("!", cd8_ccr7_path,
                                                         "&", cd8_cd45ra_path))))),
           parent = "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/8+", name = "TEMRA")
gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("!", cd8_ccr7_path,
                                                         "&!", cd8_cd45ra_path))))),
           parent = "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/8+", name = "TEM")

# Add activation gate (HLA-DR+CD38+) under CD4, CD8, GD, Vd2+GD+, Vd2-GD+, MR1, and iNKT
# First, copy activation gates up to the Live gate.
cd3_hladr_path <- "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/HLADR+"
cd3_cd38_path <- "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/CD38+"
gs_pop_add(gs, lapply(gs, gh_pop_get_gate, y=cd3_hladr_path),
           parent = "/Time/S/Live")
gs_pop_add(gs, lapply(gs, gh_pop_get_gate, y=cd3_cd38_path),
           parent = "/Time/S/Live")
live_hladr_path <- "/Time/S/Live/HLADR+"
live_cd38_path <- "/Time/S/Live/CD38+"

for (parent_path in c("/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/4+", "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/8+",
                      "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/GD+", "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/GD+/Vd2+", "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/GD+/Vd2-",
                      "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/MAIT", "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/iNKT")) {
  gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                                 list(v = as.symbol(paste0("", cd3_hladr_path,
                                                           "&", cd3_cd38_path))))),
             parent = parent_path, name = "HLADR+CD38+")
}

png(here::here(sprintf("out/QC/Final_GatingTree_All_Batches_%s.png", date)), width = 7, height = 5, units = "in", res = 300)
plot(gs[[1]], bool=T, fontsize=15)
dev.off()

recompute(gs, live_path)

# Combine into a single GatingSet to facilitate extracting data later for t-SNE/UMAP
gs2 <- gslist_to_gs(gs) # This takes a while

save_gs(gs2, here::here("out/GatingSets/20200812_HAARVI_DURT_GatingSet_AllBatches"))
# gs <- load_gs(here::here("out/GatingSets/20200812_HAARVI_DURT_GatingSet_AllBatches"))

######################################

# Make ggcyto plots for CD38 and HLADR gates on the Live node
# recompute(gs, "/Time/S/Live")

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
    theme_bw(base_size=28) +
    theme(
      legend.position = "none",
      strip.text.x = element_text(margin = margin(0,0,0,0, "cm")),
      panel.grid.major = ggplot2::element_blank()) +
    geom_stats(size=8, alpha=0.4) +
    labs(title=myGates[[1]])
}

for(b in 1:3) {
  currentGates <- c("/Time/S/Live/CD38+", "/Time/S/Live/HLADR+")
  png(filename = file.path(here::here("out/QC/FACS_Plots"),
                           sprintf("%s_on_Live_B%s_%s.png",
                                   sub(".*\\/([^\\/]+$)", "\\1", currentGates[[1]]),
                                   b, format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
      width = 1900, height = 1030)
  print(plotter(subset(gs, Batch == b),
                currentGates) +
          labs(caption = sprintf("Batch %s", b)))
  dev.off()
}

####################

# Re-plot the iNKT gate for B2
png(filename = file.path(here::here("out/QC/FACS_Plots"),
                         sprintf("%s_B2_Updated_%s.png",
                                 sub(".*\\/([^\\/]+$)", "\\1", "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/iNKT"),
                                 format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
    width = 1900, height = 1030)
print(plotter(subset(gs, Batch == 2), "/Time/S/Live/CD14-CD19-/LD/Singlet/CD3+/iNKT"))
dev.off()
