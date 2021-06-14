
# libraries ----
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

# load complete DEGs & metadata ----
  complete_DEGs <- read.csv(file = "3. Differential expression analysis/Upset_plots/Data_input/complete_DEGs.csv", 
                          stringsAsFactors = FALSE, header = TRUE)
  metadata <- read.csv(file = "3. Differential expression analysis/Upset_plots/Data_input/metadata.csv", stringsAsFactors = FALSE, header = TRUE)

  head(complete_DEGs)

# subset upregulated DEG's ----
  up <- complete_DEGs[complete_DEGs$logFC >= 0, ]
  head(up)
  up<-up[,-c(1,3:7)]

# format df 
  upset_degs_up<- up %>%mutate(value=1) %>% spread(sample, value, fill = 0)

# upset plot
    upset(upset_degs_up, nsets = 15, boxplot.summary = c("gene"),
        queries = list(list(query = intersects, params = list("Bone.4", "Bone.3"), #
                            color = "darkorchid2", active = T),
                       list(query = intersects, params = list("Bone.3", "Bone.2"), #
                            color = "darkorchid2", active = T),
                       list(query = intersects, params = list("Bone.4", "Bone.1"), #
                            color = "darkorchid2", active = T),
                       list(query = intersects, params = list("Bone.4", "Bone.2"), #
                            color = "darkorchid2", active = T),
                       list(query = intersects, params = list("Bone.3", "Bone.1"), #
                            color = "darkorchid2", active = T),
                       list(query = intersects, params = list("Bone.1", "Bone.2"), #
                            color = "darkorchid2", active = T),
                       list(query = intersects, params = list("Adrenal.2", "Adrenal.3"), #
                            color = "orange", active = T),
                       list(query = intersects, params = list("Adrenal.2", "Adrenal.1"), #
                            color = "orange", active = T),
                       list(query = intersects, params = list("Adrenal.3", "Adrenal.1"), #
                            color = "orange", active = T),
                       list(query = intersects, params = list("Adrenal.1", "Adrenal.2", "Adrenal.3"), #
                            color = "orange", active = T),
                       list(query = intersects, params = list("Brain.1", "Brain.3"), #
                            color = "forestgreen", active = T),
                       list(query = intersects, params = list("Brain.1", "Brain.2"), #
                            color = "forestgreen", active = T),
                       list(query = intersects, params = list("Brain.3", "Brain.2", "Brain.1"), #
                            color = "forestgreen", active = T),
                       list(query = intersects, params = list("Brain.3", "Brain.2"), #
                            color = "forestgreen", active = T)),
        order.by = "degree",
        mainbar.y.label = "Intersecting Gene Count",
        sets.x.label = "Unique Gene Count",
        sets = c("Bone.6", "Bone.5",
                 #"Kidney.1",
                 "Brain.4",
                 "Bone.4", "Bone.3", "Bone.2", "Bone.1",
                 "Brain.3","Brain.2", "Brain.1",
                 "Adrenal.3","Adrenal.2","Adrenal.1"),
        # set_size.show = TRUE,
        text.scale = c(1.8, 1.8, 1.8, 1.8, 1.8, 1.8),
        # c(1: intersection size title, 2: intersection size
        # 3: tick labels, 4: set size title, set size tick labels, set names, numbers above bars)
        mb.ratio = c(0.55,0.45),
        point.size = 3,
        line.size = 1.1,
        matrix.color = "black",
        sets.bar.color = c("darkorchid2", "darkorchid2",
                           #"black", 
                           "forestgreen",
                           "darkorchid2", "darkorchid2", "darkorchid2", "darkorchid2",
                           "forestgreen", "forestgreen", "forestgreen",
                           "orange", "orange", "orange"),
        main.bar.color = "black",
        shade.color = "black",
        shade.alpha = 0.1,
        matrix.dot.alpha = 0.3,
        color.pal = 1,
        keep.order = TRUE,
        set.metadata = list(data=metadata, plots=list(
          # list(type="text", column = "Cell_line", assign = 14),
          list(type="heat", column = "Estrogen_status", assign = 10,
               colors = c(ED = "lightblue", E2 = "lavenderblush2")))
        )
    )

  
# subset downregulated DEG's ----
  down <- complete_DEGs[complete_DEGs$logFC <= 0, ]
  head(down)
  down<-down[,-c(1,3:7)]

# format df 
  down_trans<- down %>%mutate(value=1) %>% spread(sample, value, fill = 0)
  
# upset plot
    upset(down_trans, nsets = 15, boxplot.summary = c("gene"),
          queries = list(list(query = intersects, params = list("Bone.4", "Bone.3"), #
                              color = "darkorchid2", active = T),
                         list(query = intersects, params = list("Bone.3", "Bone.2"), #
                              color = "darkorchid2", active = T),
                         list(query = intersects, params = list("Bone.4", "Bone.2"), #
                              color = "darkorchid2", active = T),
                         list(query = intersects, params = list("Bone.3", "Bone.1"), #
                              color = "darkorchid2", active = T),
                         list(query = intersects, params = list("Adrenal.2", "Adrenal.3"), #
                              color = "orange", active = T),
                         list(query = intersects, params = list("Adrenal.2", "Adrenal.1"), #
                              color = "orange", active = T),
                         list(query = intersects, params = list("Adrenal.3", "Adrenal.1"), #
                              color = "orange", active = T),
                         list(query = intersects, params = list("Adrenal.1", "Adrenal.2", "Adrenal.3"), #
                              color = "orange", active = T),
                         list(query = intersects, params = list("Brain.1", "Brain.3"), #
                              color = "forestgreen", active = T),
                         list(query = intersects, params = list("Brain.3", "Brain.2"), #
                              color = "forestgreen", active = T)),
          order.by = "degree",
          mainbar.y.label = "Intersecting Gene Count", 
          sets.x.label = "Downregulated DEG's",
          sets = c("Bone.6", "Bone.5",
                   "Brain.4",
                   "Bone.4", "Bone.3", "Bone.2", "Bone.1",
                   "Brain.3","Brain.2", "Brain.1",
                   "Adrenal.3","Adrenal.2","Adrenal.1"),
          # set_size.show = TRUE,
          text.scale = c(1.8, 1.8, 1.8, 1.8, 1.8, 1.8),
          mb.ratio = c(0.55,0.45),
          point.size = 3,
          line.size = 1.1,
          matrix.color = "slategray",
          sets.bar.color = c("darkorchid2", "darkorchid2",
                             "forestgreen",
                             "darkorchid2", "darkorchid2", "darkorchid2", "darkorchid2",
                             "forestgreen", "forestgreen", "forestgreen",
                             "orange", "orange", "orange"),
          main.bar.color = "black",
          shade.color = "black",
          shade.alpha = 0.1,
          matrix.dot.alpha = 0.3,
          color.pal = 1,
          keep.order = TRUE
          # set.metadata = list(data=metadata, plots=list(
          #   # list(type="text", column = "Cell_line", assign = 14),
          #   list(type="heat", column = "Estrogen_status", assign = 10,
          #        colors = c(ED = "lightblue", E2 = "lavenderblush2")))
          # )
    )
    
      
    
    
  
  