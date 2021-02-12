

# load libraries

library(plyr)
library(limma)
library(UpSetR)
library(reshape2)
library(dplyr)
library(tidyr)
library(viridis)
library(RColorBrewer)
library(ggplot2)
library(ggsci)
library(rlang)
library(biomaRt)

upset_plot_up <- read.csv(file = "final_upset_data/Finalized/upset_input_A_trial.csv",
                          stringsAsFactors = FALSE,
                          header = TRUE)

gene_up<- upset_plot_up %>%mutate(value=1) %>% spread(sample, value, fill = 0)

metadata <- read.csv(file = "final_upset_data/Finalized/metadata.csv", stringsAsFactors = FALSE, header = TRUE)

upset(gene_up, nsets = 15, boxplot.summary = c("gene"),
      queries = list(list(query = intersects, params = list("Bone.4", "Bone.3"),
                          color = "darkorchid2", active = T),
                     list(query = intersects, params = list("Bone.3", "Bone.2"),
                          color = "darkorchid2", active = T),
                     list(query = intersects, params = list("Bone.4", "Bone.1"),
                          color = "darkorchid2", active = T),
                     list(query = intersects, params = list("Bone.4", "Bone.2"),
                          color = "darkorchid2", active = T),
                     list(query = intersects, params = list("Bone.3", "Bone.1"),
                          color = "darkorchid2", active = T),
                     list(query = intersects, params = list("Adrenal.2", "Adrenal.3"),
                          color = "orange", active = T),
                     list(query = intersects, params = list("Adrenal.2", "Adrenal.1"),
                          color = "orange", active = T),
                     list(query = intersects, params = list("Adrenal.3", "Adrenal.1"),
                          color = "orange", active = T),
                     list(query = intersects, params = list("Adrenal.1", "Adrenal.2", "Adrenal.3"),
                          color = "orange", active = T),
                     list(query = intersects, params = list("Brain.1", "Brain.3"),
                          color = "forestgreen", active = T),
                     list(query = intersects, params = list("Brain.1", "Brain.2"),
                          color = "forestgreen", active = T),
                     list(query = intersects, params = list("Brain.3", "Brain.2", "Brain.1"),
                          color = "forestgreen", active = T),
                     list(query = intersects, params = list("Brain.3", "Brain.2"),
                          color = "forestgreen", active = T)),
      order.by = "degree",
      mainbar.y.label = "Intersecting Gene Count", 
      sets.x.label = "Upregulated DEG count",
      sets = c("Bone.6", "Bone.5",
               "Kidney.1",
               "Brain.4",
               "Bone.4", "Bone.3", "Bone.2", "Bone.1",
               "Brain.3","Brain.2", "Brain.1",
               "Adrenal.3","Adrenal.2","Adrenal.1"),
      set_size.show = TRUE,
      text.scale = c(1.5, 1, 1.5, 1.5, 1.5, 1.5),
      # c(1: intersection size title, 2: intersection size
      # 3: tick labels, 4: set size title, set size tick labels, set names, numbers above bars)
      mb.ratio = c(0.6,0.4),
      point.size = 3,
      line.size = 1.1,
      matrix.color = "slategray",
      sets.bar.color = c("darkorchid2", "darkorchid2",
                         "black", 
                         "forestgreen",
                         "darkorchid2", "darkorchid2", "darkorchid2", "darkorchid2",
                         "forestgreen", "forestgreen", "forestgreen",
                         "orange", "orange", "orange"),
      main.bar.color = "slategray",
      shade.color = "slategray",
      shade.alpha = 0.1,
      matrix.dot.alpha = 0.3,
      color.pal = 1,
      keep.order = TRUE,
      set.metadata = list(data=metadata, plots=list(
        list(type="text", column = "Cell_line", assign = 14),
        list(type="heat", column = "E2_status", assign = 14,
             colors = c(ED = "lightblue", E2 = "lavenderblush2"))
      ))
)
