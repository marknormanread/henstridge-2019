library(stats)
library(mixOmics)
library(phyloseq)
library(vegan)  # Adonis function.  


#### Number of sequences per group, with instances (samples) explicitly shown and bar graphs indicating distribution. 
sample_sequencing_depth_by_groups = function(
  # Each of these are 
  sample_names, # Vector of string, one item per instance/sample.
  groups,  # Vector of class for each instance. Tested with type factor. Other types might work too. 
  num_reads,  # Vector of numerical values, one item per insteance. 
  data_path,  # Where to write data. Omit extentions, they are added here. 
  col_per_group  # Vector of named items, c('group1' = '#E9A023', ...)
  ) 
{
  sample_depth_df = data.frame(sample_names, groups, num_reads)
  p = ggplot(sample_depth_df, aes(x=groups, y=num_reads, color=groups)) +
      geom_boxplot(outlier.shape=NA) +  # Turn off outliers, as we plot each datapoint below. 
      geom_jitter(shape=16, position=position_jitter(0.35), size=2) +
      ylab('Number of sequences') +
      xlab('Groups') +
      # scale_colour_brewer(palette=rainbow(n=9)) +
      scale_colour_manual(values=col_per_group) +
      theme_gray() +
      theme(text = element_text(size=9)) +
      # theme(axis.text.x=element_blank()) + # Turn of x axis tick labels  
      theme(legend.position="none")
  print(p)
  
  ggsave(paste(data_path, '.pdf', sep=''), width=8, height=5.5, units='cm')
  
  # Stats on these data using PERMANOVA.
  print(adonis(num_reads ~ groups, data=sample_depth_df, method="euclidean"))  # Print to the notebook. 
  # Same again, but save to the file system.
  sink(paste(data_path, '.txt', sep=''))
  print(adonis(num_reads ~ groups, data=sample_depth_df, method="euclidean"))
  sink()
  
  return(p)
}


sample_richness_by_group = function(
  phylo,  # PhyloSeq object (requires sample-by-taxa table, taxonomy map and metadata - this can be an arse to prepare without errors)
  group_name,  # String, EXACTLY as it appears in the metadata file. 
  data_path,  # Where to write data. Omit extentions, they are added here. 
  col_per_group,  # Vector of named items, c('group1' = '#E9A023', ...)
  plot_means = FALSE
)
{
  # Plot richness
  # Ignore the warning about singletons - pretty sure this is because we use DADA2 which doesn't produce as many 
  # singletons as previous technologioes. 
  p = plot_richness(phylo, x = group_name, measures = c("Observed","Shannon","InvSimpson"), color = group_name)
  # Change the order of the samples on the x axis to be listed by experimental group. This was manually determined in the mapping file: "SortOrder"
  p$data$samples = as.character(p$data$samples)  # "Sample" here is data-set specific. It's the dataframe column being plotted on X
  p$data$samples = factor(p$data$samples, levels=levels(sample_data(phylo)$Sample))
  p = p + geom_point(size=1.5) +
      theme(legend.position="none") + 
      theme(text = element_text(size=9)) +
      theme(axis.title.x = element_blank()) + # Turn off x axis tick labels
      # theme(axis.text.x=element_blank()) + # Turn off x axis tick labels
      # scale_color_brewer(palette = 'rainbow')
      scale_colour_manual(values=col_per_group)
  if (plot_means) {
    p = p + stat_summary(  # Plot means. 
      geom = "point",
      fun.y = "mean",
      col = "black",
      size = 2,
      shape = 1,
      fill = NA
    )
  }
  print(p)  # Display graph. 
  ggsave(paste(data_path, '.pdf', sep=''), width=8, height=5.5, units='cm')  # Adjust the width of the saved plot, stop sample IDs being squeezed. 

  richness = estimate_richness(phylo)  # Actually used in the plot_richness function above. Same data. 
  richness[group_name] = sample_data(phylo)[, group_name]
  
  # Perform a PERMANOVA to seek statistically significant differences in groups, in terms of alpha diversity. 
  # adonis(eval(parse(text=paste('Observed ~ ', group_name, sep=''))), data=richness, method="euclidean")
  print(adonis(eval(parse(text=paste('Observed ~ ', group_name, sep=''))), data=richness, method="euclidean"))  # Print to the notebook.
  print(adonis(eval(parse(text=paste('Shannon ~ ', group_name, sep=''))), data=richness, method="euclidean"))  # Print to the notebook.
  print(adonis(eval(parse(text=paste('InvSimpson ~ ', group_name, sep=''))), data=richness, method="euclidean"))  # Print to the notebook.
  # Repeat, but save to file system.
  sink(paste(data_path, '_permanova.txt', sep=''))
  cat('Observed taxa\n')
  print(adonis(eval(parse(text=paste('Observed ~ ', group_name, sep=''))), data=richness, method="euclidean"))  # Print to the notebook.
  cat('\n\nShannon taxa\n')
  print(adonis(eval(parse(text=paste('Shannon ~ ', group_name, sep=''))), data=richness, method="euclidean"))  # Print to the notebook.
  cat('\n\nInverse Simpson taxa\n')
  print(adonis(eval(parse(text=paste('InvSimpson ~ ', group_name, sep=''))), data=richness, method="euclidean"))  # Print to the notebook.
  sink()
}

#### How well separated are the samples from two different experimental groups?
# Two groups, A and B. 
# Calculates all pairwise distances between samples in A with themselves. Same with samples from B. 
# These compbine to give the 'within group' distribution. 
# Then calculate the pairwise distances between samples from A with those from B: 'between group' distribution. 
# These are plotted as probability density functions. 
# Kolmogorov smirnov test conducted on the distributions to gauge statistical significance of any separations. 
#
# Mark N. Read, 2019. 
within_between_group_distances = function(
  # Euclidean distance between CLR-transformed data recommended, as it copes well with compositional data. 
  distance_matrix,  # n x n matrix of distances between items. 
  group_vector,  # Eg column from metadata file. Order of items must match those in distance matrix.
  plotsave_prefix='',  # To save to a different directory, provide path through this prefix.
  # If TRUE, plot separate distributions for within groups A and B distances. Otherwise, combine them.
  separate_within_distributions=FALSE,
  titles=TRUE,
  xlimits=NULL,  # or c(lower_limit, upper_limit)
  plot_size = c(34, 22)  # mm
)
{
  dir.create(plotsave_prefix, showWarnings = FALSE)  # Creat if it doesn't already exist. 
  
  group_set = as.vector(unique(group_vector))  # Contrast only the groups in the supplied vector
  # Store pairwise group comparisons in terms of kolmogorov smirnov D and associated  p values in these nxn dataframes.
  group_differences_ks = data.frame(matrix(nrow=length(group_set), ncol=length(group_set)))
  group_differences_p = data.frame(matrix(nrow=length(group_set), ncol=length(group_set)))
  
  rownames(group_differences_ks) = group_set
  colnames(group_differences_ks) = group_set
  rownames(group_differences_p) = group_set
  colnames(group_differences_p) = group_set
  
  pairwise_combinations = combn(unique(group_set), 2)
  num_comparisons = ncol(pairwise_combinations)
  
  for (comparison in 1:num_comparisons)
  {
    group_a = pairwise_combinations[1, comparison]
    group_b = pairwise_combinations[2, comparison]
    
    a_vs_a = distance_matrix[group_vector == group_a, # Each of these produces a boolean vector.
                           group_vector == group_a]
    mask = lower.tri(a_vs_a)  # Retruns a boolean mask of dim(a_vs_a). Removes the diagonal too.
    a_vs_a_vect = a_vs_a[mask]  # Transforms to a vector.
    
    b_vs_b = distance_matrix[group_vector == group_b, group_vector == group_b]
    b_vs_b_vect = b_vs_b[lower.tri(b_vs_b)]  
    
    a_vs_b = distance_matrix[group_vector == group_a, group_vector == group_b]
    # DO NOT take the lower triangle only, all distances here need to be considered - no duplicate values. 
    between_group_distances = a_vs_b[TRUE]  # Transforms to a vector. 
    within_group_distances = c(a_vs_a_vect, b_vs_b_vect)  # Concatenate vectors
    
    ks = ks.test(x=within_group_distances, y=between_group_distances, 
                 alternative=c('greater'))  # Hypothesise that within (x) not greater than between (y)
    # Adjust for multiple comparisons. 
    p_value = p.adjust(ks$p.value, method='bonferroni', n=num_comparisons)
    
    if (separate_within_distributions == FALSE) {
      # ggplot needs data in dataframe format
      df = data.frame(Distance=c(within_group_distances, between_group_distances), 
                      Condition=c(rep('Within', length(within_group_distances)),
                                  rep('Between', length(between_group_distances))))
      p = ggplot(df, aes(x=Distance, fill=Condition)) + 
          # There are other options; this shows subtble differences in distribution.  
          # http://www.cookbook-r.com/Graphs/Plotting_distributions_(ggplot2)/
          geom_density(alpha=.3) +
          theme(plot.background = element_rect(fill = "transparent", color = NA)) + # Transparent plot background. 
          coord_cartesian(xlim = xlimits)
    } else {
      # ggplot needs data in dataframe format
      df = data.frame(Distance=c(within_group_distances, between_group_distances), 
                      Condition=c(rep(paste('Within ', group_a, sep=''), length(a_vs_a_vect)),
                                  rep(paste('Within ', group_b, sep=''), length(b_vs_b_vect)),
                                  rep('Between', length(between_group_distances))))
      p = ggplot(df, aes(x=Distance, color=Condition)) + geom_density(alpha=.3) +
          theme(plot.background = element_rect(fill = "transparent", color = NA))  # Transparent plot background. 
    }
    
    if (! is.null(xlimits)) {
      p + xlim(xlimits[1], xlimits[2])
    }
    if (titles) {
      p + ylab('Probability density') +
            xlab('Distance') +
            ggtitle(paste('Groups: ', group_a, ' and ', group_b,
                          # Include stats on title.
                          '; D=', round(ks$statistic, digits=2), ', p=', round(p_value, digits=2), sep=''))
      ggsave(paste(plotsave_prefix, 'group_distances_', group_a, '_vs_', group_b, '.pdf', sep=''))
    } else {
      # This completely removes ALL space and items around the polot itself. 
      p + theme(legend.position="none",
                axis.title = element_blank(),  # Turn off axis titles.
                axis.text = element_blank(),  # Turn off tick label titles.
                axis.ticks = element_blank(),  # Turn off black tick marks.
                text = element_blank(),
                panel.border = element_blank(),
                plot.margin=grid::unit(c(0,0,0,0),"cm"),
                axis.ticks.length = unit(0, "pt")
                )
      ggsave(paste(plotsave_prefix, 'group_distances_', group_a, '_vs_', group_b, '.pdf', sep=''),
             width=plot_size[1], height=plot_size[2], units='mm')
    }

    
    # Write to dataframe of results
    group_differences_ks[group_b, group_a] = ks$statistic
    group_differences_p[group_b, group_a] = p_value
    # Put in the upper triangle too. It's redundant data, but easier for look up. 
    group_differences_ks[group_a, group_b] = group_differences_ks[group_b, group_a]
    group_differences_p[group_a, group_b] = group_differences_p[group_b, group_a]
  }
  return(list(group_differences_ks=group_differences_ks, 
              group_differences_p=group_differences_p))
}


#### Compositional data-robust distance calculations. 
# Conventional distance metrics (Bray Curtis, UniFrac, etc) have been criticised for introducing artefacts under 
# compisitional data, as microbiome sequencing data is. 
# Gloor et al advocate the 'Aitchison' distance metric instead.
# (https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224/full) 
# This is Euclidean distance in centred logratio transformed data.
aitchison_distance = function(taxa_matrix,  # Taxa/Features as columns, samples as rows
                              lr='CLR'  # Can also be 'ILR' for Isometric Logratio Transform. 
                              ) 
{
  taxa_clr = mixOmics::logratio.transfo(taxa_matrix,  # Samples as rows, Taxa/Features as columns. 
                                        logratio=lr, offset=0)
  
  clr_dist_raw = dist(taxa_clr, method='euclidean', diag=FALSE)  # Returns some weird 'dist' class ...
  # ... turn it into something better behaved, with functioning colnames and rownames. 
  clr_dist = as.matrix(clr_dist_raw) 
  return (clr_dist)
}
