# Assembly  of Figure 1 a) points distribution map b) three predicted abundance maps 
# in forest and pasture

library(magick)
library(ggplot2)
library(cowplot)

setwd("F:/Doctorado/Tesis/figura_1")

fig_a <- image_read("./fig_1A.jpeg")
fig_b <- image_read("./fig_1B.jpg")
fig_c <- image_read("./fig_1C.jpeg")
fig_d <- image_read("./fig_1D.jpg")
fig_A <- ggdraw() + draw_image(fig_a)
fig_B <- ggdraw() + draw_image(fig_b)
fig_C <- ggdraw() + draw_image(fig_c)
fig_D <- ggdraw() + draw_image(fig_d)

final_plot_cowplot <- plot_grid(
  plot_grid(fig_A, fig_B, fig_C, fig_D, ncol = 2, labels = c("a", "b", "c", "d"), 
            label_size = 12, label_fontface = "bold"))

ggsave("fig_1_rev.jpg", final_plot_cowplotA, width = 8.5, height = 10, units = "in", dpi = 300)

# enssambly

final_plot <- ggdraw() +
  draw_plot(fig_A, x = 0, y = 0.5, width = 2/3, height = 0.5) +
  draw_label("a", x = 0.02, y = 0.97, fontface = "bold", size = 14) +
  draw_plot(fig_B, x = 2/3, y = 0.5, width = 1/3, height = 0.5) +
  draw_label("b", x = 0.68, y = 0.97, fontface = "bold", size = 14) +
  draw_plot(fig_C, x = 0.1, y = 0, width = 0.9, height = 0.35) +
  draw_label("c", x = 0.02, y = 0.47, fontface = "bold", size = 14)

ggsave("draft1.jpg", final_plot, width = 8.5, height = 8, units = "in", dpi = 300)
