library(tidyverse, warn.conflicts = FALSE)
library(ggrastr)
library(ggthemes)
library(patchwork)

idx <- commandArgs(trailingOnly=TRUE)
n_int_AB <- 5
n_int_ABC <- 5

print(idx)

df_p1 <- read_csv(paste0('tmp/test_', idx, '_p1.csv'), show_col_types = FALSE) 

df_mod_p1 <- df_p1 %>%
    filter(position != -9) %>%
    pivot_longer(-position) %>%
    separate_wider_delim(name, ', ', names = c('topology', 'int_1', 'int_2')) %>%
    mutate(
        topology = as.integer(topology),
        int_1 = as.integer(int_1),
        int_2 = as.integer(int_2),
        is_V0 = topology == 0
    )

df_p2 <- read_csv(paste0('tmp/test_', idx, '_p2.csv'), show_col_types = FALSE)

df_mod_p2 <- df_p2 %>%
    filter(position != -9) %>%
    pivot_longer(-position) %>%
    separate_wider_delim(name, ', ', names = c('topology', 'int_1', 'int_2')) %>%
    mutate(
        topology = as.integer(topology),
        int_1 = as.integer(int_1),
        int_2 = as.integer(int_2),
        is_V0 = topology == 0
    )
    
p1 <- df_mod_p1 %>%
    group_by(position, is_V0, int_1) %>%
    summarize(prob = sum(value)) %>%
    ggplot() +
    rasterize(geom_tile(aes(position, int_1+(!is_V0)*(n_int_AB), fill = prob, color = prob)), dpi = 300) +
    geom_hline(aes(yintercept = n_int_AB-1+0.5), color = 'white') +
    theme_few() +
    scale_fill_viridis_c(name = 'Posterior\nprobability', 
                            limits = c(0, 1),
                            option="inferno"
                        ) +
    scale_color_viridis_c(name = 'Posterior\nprobability', 
                            limits = c(0, 1),
                            option="inferno"
                            ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(
        breaks = c(0:(n_int_AB-1), ((n_int_AB):(n_int_AB+n_int_ABC-1))), 
        # labels = c(0:(n_int_AB-1), 0:(n_int_ABC_1-1)),
        labels = paste0('F', c(0:(n_int_AB+n_int_ABC-1))),
        expand = c(0, 0),
        sec.axis = dup_axis(
            # breaks = c(0:(n_int_AB-1), ((n_int_AB):(n_int_AB+n_int_ABC_1)))-0.5,
            # labels = sprintf("%.0f", round(c(cut_AB[0:(length(cut_AB)-1)], cut_ABC_new)))
            # labels = rep('', length(round(c(cut_AB[0:(length(cut_AB)-1)], cut_ABC_new))))
        )
    ) +
    theme(
        axis.title.y.right = element_blank()
    ) +
    labs(y = 'First coalescent', x = 'Position')

p2 <- df_mod_p2 %>%
    group_by(position, int_2) %>%
    summarize(prob = sum(value)) %>%
    ggplot() +
    rasterize(geom_tile(aes(position, int_2, fill = prob, color = prob)), dpi = 300) +
    theme_few() +
    scale_fill_viridis_c(name = 'Posterior\nprobability', 
                            limits = c(0, 1),
                            option="inferno") +
    scale_color_viridis_c(name = 'Posterior\nprobability', 
                            limits = c(0, 1),
                            option="inferno") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(
        breaks = 0:(n_int_ABC-1), 
        # labels = 0:(n_int_ABC_2-1),
        labels = paste0('S', 0:(n_int_ABC-1)),
        expand = c(0, 0),
        sec.axis = dup_axis(
            # breaks = (0:n_int_ABC_2)-0.5,
            # labels = sprintf("%.0f", round(cut_ABC))
            # breaks = seq(-0.5, n_int_ABC_2-0.5, 0.5),
            # labels = c(rbind(rep('', length(cut_ABC)-1), paste0('S', 0:(n_int_ABC_2-1))), '')
        )
    ) +
    labs(y = 'Second coalescent') +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y.right = element_blank()
    )

p3 <- df_mod_p1 %>%
    group_by(position, topology) %>%
    summarize(prob = sum(value)) %>%
    ggplot() +
    rasterize(geom_tile(aes(position, topology, fill = prob, color = prob)), dpi = 300) +
    theme_few() +
    scale_fill_viridis_c(name = 'Posterior\nprobability', 
                            limits = c(0, 1),
                            option="inferno") +
    scale_color_viridis_c(name = 'Posterior\nprobability', 
                            limits = c(0, 1),
                            option="inferno") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(
        breaks = c(0, 1, 2, 3), 
        labels = c('V0', 'V1', 'V2', 'V3'),
        expand = c(0, 0),
        sec.axis = dup_axis()
    ) +
    labs(y = 'Topology') +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y.right = element_blank()
    )

fac <- 0.7 

plt_tot <- (p3+p2+p1 + 
    plot_layout(ncol = 1, heights = c(4, n_int_ABC, n_int_AB+n_int_ABC), guides = 'collect')) & theme(legend.position = 'none')

ggsave(paste0('./tmp/tmp_', idx, '.pdf'), plt_tot, 
        width = 14*fac, height = 9*fac)