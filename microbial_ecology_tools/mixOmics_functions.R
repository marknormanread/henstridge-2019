library(matrixStats)
library(ggbeeswarm)

#### Filter low abundance features out of the dataset. 
# Adapted from http://mixomics.org/mixmc/pre-processing/
# Function to perform pre-filtering. 
# To be retained, a feature must comprise at least the given percentage of all counts found across all samples (pooled). 
low_count_removal = function(
  data,  # feature count data frame of size n (sample, rows) x p (features, cols). Features as columns. 
  percent=0.01)  # Cutoff chosen
{
  keep_features = which(colSums(data) * 100 / (sum(colSums(data))) > percent)
  data_filter = data[, keep_features]
  return(list(data_filter = data_filter, keep_features = keep_features))
}


#### Total sum scaling (TSS) normalisation
# Each variable read count is divided by the total number of read counts. 
# Use on single columns at a time, NOT the whole table. 
# E.g. table_relab = t(apply(table, 1, TSS_normalisation)) where tables place features as columns. 
TSS_normalisation = function(x) { x/sum(x) }


#### Useful to run before choosing filter cutoff. 
# What is the distribution of percentages of total reads captured by the OTUs?
plot_feature_abundances = function(feature_table,  # Features as columns, samples as rows. Counts (not rel. abundance).
                                   xlab = 'Features, sorted by % total counts',
                                   ylab = 'Percent of total counts (%)',
                                   plot_log = TRUE,
                                   plot_ggplot = FALSE,
                                   filter_relative_abundance = -Inf,  # Exclude features with percentate less than this. 
                                   filter_n_abundant = Inf  # Include the top n.
                                   )
{
  feature_percent_all_counts = (100 * colSums(feature_table)) / sum(colSums(feature_table))
  feature_percent_all_counts = sort(feature_percent_all_counts, decreasing=TRUE)
  # Filter out ASV capturing less than this given relative abundance
  feature_percent_all_counts = feature_percent_all_counts[feature_percent_all_counts >= filter_relative_abundance]
  # Show only top n ASV by relative abundance. 
  feature_percent_all_counts = feature_percent_all_counts[1: min(length(feature_percent_all_counts), filter_n_abundant)]
  if (plot_ggplot)
  {
    df = data.frame(feature_percent_all_counts, rank=1:length(feature_percent_all_counts))
    p = ggplot(df, aes(x=rank, y=feature_percent_all_counts)) +
        geom_point() +
        ylab(ylab) +
        xlab(xlab)
    return(p)
  }
  else {
    plot(feature_percent_all_counts, xlab=xlab, ylab=ylab)
    if (plot_log)
      plot(feature_percent_all_counts, xlab=xlab, ylab=ylab, log='y')
  }
}

#### Helps in selecting cutoff threshold for filtering feature-instance tables. 
# Each feature is a datapoint. 
# Plots feature relative abundance (where you would be filtering) against cumulative count of reads captured as a %.
# Hence, you can see (for example) that retaining features with >0.1% relative abundance of whole dataset 
#   captures 90% of all reads in the dataset. 
plot_feature_cumulative_abundance = function(
  feature_table,  # Features as columns, samples as rows. Counts (not rel. abundance).
  xlab = 'Feature relative abundance',
  ylab = 'Total counts captured (%)',
  plot_log = 'x',
  percent_cutoff = NULL
  ) 
{
  feature_ids = rownames(feature_table)
  
  total_reads = sum(feature_table)
  feature_read_count = c()  # Named vector (Dictionary), [Feature ID] = # reads feature captures. 
  for (f_id in feature_ids) {
    feature_read_count[f_id] = sum(feature_table[f_id, ])  # Counts captured by this feature across all samples. 
  }  
  feature_read_count = rev(sort(feature_read_count))
  feature_read_count_relab = 100 * feature_read_count / total_reads
  feature_read_count_cumulative = 100 * cumsum(feature_read_count)  # Cumulative count
  feature_read_count_cumulative = feature_read_count_cumulative / total_reads  # Turn into proportion of all counts observed. 
  plot(x = feature_read_count_relab, y = feature_read_count_cumulative, 
       xlab = xlab, 
       ylab = ylab,
       log = plot_log)
  if (! is.null(percent_cutoff)) {
    abline(v = percent_cutoff)
  }
}

#### Transforms a raw feature-sample (as rows-cols) table into transposed relative abundance (total sum scaling).
# mixOmics adheres to traditional format of features as columns, however many microbial ecology bioinformatics 
# pipelines place them as rows. 
# Additionally, this prepares a pseudocount, needed for centred/isometric log ratio transforms, and filters. 
# Consider if filtering is important before performing it. 
# In targetted applications (looking at specific bugs or pathways) filtering may not be needed. 
# For exploratory analyses with many features it can be beneficial in disuading overfitting. 
prepare_mixomics = function(
  table,  # Supply as samples across columns and taxa/genes/predictors/features as rows. 
  pseudocount = 1,  # Select this to be small w.r.t. the magnitude of data in `table`. 
  filter_percent = 0.01,  # Remove features representing < this value of all counts.
  # Downstream analysis needs observations/taxa/features as columns and samples as rows. 
  # Transpose the input table if it is not already in this format. 
  transpose_table = TRUE  
  )
{
  # Notes. I have wondered if it is best to do the filtering before or after the relative abundance calculation. 
  # It is possible to remove quite a lot of features through stringent filtering. For datasets where there are not
  # many samples there can be an arguement for doing so, too many features relative to the number of instances
  # can lead to overfitting. Under stringent filtering, the relative abundance is then calculated a comparatively 
  # small subset of the total counts (e.g., 60% is readily achievable). Is this still really relative abundance?
  # CLR transforms (and I imagine ILR transforms, though I haven't explicitly checked that) are dependent on all the 
  # features included in the presented dataset. Post relative abundance calculations, a full and filtered dataset 
  # will give different CLR values. Hence, to ensure that CLR values are not altered when filtering features we would
  # need to calculate the CLR-transformed data BEFORE applying any filtering. This also means that the mixOmics 
  # tools are given pre-CLR-transformed data, and not instructed to perform the CLR themselves on filtered data. 
  table_mixomics = table
  if (transpose_table) {
    table_mixomics = t(table)  # Transpose such that samples = rows (mixomics format, also more conventional)
  }
  table_pseudocount = table_mixomics + pseudocount  # Can't log(0), which is used in CLR transform: add pseudocount.
  
  filter_result = low_count_removal(table_pseudocount, percent=filter_percent)
  table_filtered = filter_result$data_filter
  
  table_relab = t(apply(table_filtered, 1, TSS_normalisation))
  return(table_relab)
}


#### Genereric function for plotting a PCA/sPLS-DA ordination. 
ordinate = function(projection,  # A PCA or sPLS-DA
                    group,  # Factor vector descirbing the group of each sample included in `projection`.
                    plot_title,
                    filename = NULL,  # Location to save plot. Remember to append '.pdf'. 
                    name_samples = FALSE,  # If TRUE, write name of each sample on the plot.
                    sample_plot_characters = 16,  # Single char, else vector of n=len(samples). Used if name_samples=FALSE
                    ellipse = FALSE,
                    background = NULL,
                    distance = "mahalanobis.dist", 
                    components = c(1, 2),
                    col_per_group_map = NULL,  # Named vector. `group` names used as look up for colors. 
                    color_scheme = NULL,  # e.g., use 'rainbow'
                    show_legend = TRUE,
                    dimensions_mm = c(NA, NA),
                    fontsize = NULL,
                    title_fontsize = 18,  # Overriden by fontsize, if it is specified. 
                    spartan_plot = FALSE,  # Disable axis tick marks, ticks lines,  
                    point_size = 3,
                    cex_1D = 5.0  # In case switch to 1D mode. 
                    )
{
  if (is.null(col_per_group_map))
  {
    if (length(levels(group)) > 10 | identical(color_scheme, 'rainbow'))
      col_per_group = rainbow(n=length(levels(group)))
    else
      col_per_group = color.mixo(1:length(levels(group)))
  } else {
    col_per_group = col_per_group_map[levels(group)]  # Extract the colours pertinent to the groups. 
    if (! is.factor(group))
      cat('\n\nWarning!! Supplied class labels are not factors, this is likely to cause an error.\n')
  }

  if (projection$ncomp == 1)
  {
    cat('Warning! Model has only 1 dimension, executing ordinate_1d plotting procedure instead.\n')
    return(ordinate_1d(splsda_model = projection, 
                       shape_vector = factor(sample_plot_characters), 
                       col_per_group = col_per_group_map,
                       graph_path = filename,
                       cex_1D = cex_1D))
  }
  
  if (! is.null(fontsize))
    title_fontsize = fontsize
  
  if (identical(background, TRUE) && identical(components, c(1, 2))) {
    # Not clear how to generate a background for anything but the first two components. 
    background = background.predict(projection, comp.predicted=2, dist=distance)
  }
  
  if (name_samples)
  {
    res = plotIndiv(projection, 
                    comp = components, # Select which components to plot on each axis. 
                    # This needs to be absent from the call to plot sample names. 
                    # pch = 16,  # Plot character to use
                    cex = point_size,
                    ind.names = name_samples,  # Do not name each item on the plot
                    group = group, 
                    col.per.group = col_per_group,
                    legend = show_legend,
                    ellipse = ellipse, 
                    background = background,
                    title = plot_title,
                    size.title = title_fontsize,
                    style = 'ggplot2'  # Set to return ggplot object as res$graph. Supposed to be default but doens't work. 
                    )
  } else {
    res = plotIndiv(projection,
                    comp = components, # Select which components to plot on each axis.
                    pch = sample_plot_characters,  # Plot character to use
                    cex = point_size,
                    ind.names = name_samples,  # Do not name each item on the plot
                    group = group,
                    col.per.group = col_per_group,
                    legend = show_legend,
                    ellipse = ellipse,
                    background = background,
                    title = plot_title,
                    size.title = title_fontsize,
                    style = 'ggplot2'  # Set to return ggplot object as res$graph. Supposed to be default but doens't work.
                    )
  }
  # MixOmics labels axes "X-variate". Much more meaningful to say "Component". 
  res$graph$labels$x = gsub('X-variate', 'Component', res$graph$labels$x)
  res$graph$labels$y = gsub('X-variate', 'Component', res$graph$labels$y)
  
  if (spartan_plot)
  {
    res$graph = res$graph + theme(legend.position="none",
                  text = element_text(size=fontsize),
                  axis.text = element_blank(),  # Turn off tick label titles.
                  axis.ticks = element_blank()  # Turn off black tick marks.
                  # plot.title = element_text(size=fontsize)
                  # Haven't been able to turn off the plot title completely. Stuck with the grey bar at top.
                  )
  }

  if (! is.null(filename)) {
    ggsave(filename, width=dimensions_mm[1], height=dimensions_mm[2], units='mm')
  }
  return(res$graph)
}


ordinate_1d = function(splsda_model, 
                       shape_vector = NULL,
                       col_per_group = NULL,
                       y_label = 'Groups\nseparated',
                       graph_path = NULL,  # Location on filesystem. Remember to provide file extension. 
                       width_mm = 60, height_mm = 30,
                       cex_1D = 5.0  # Spacing within groups on the swarm plot. 
                       )
{
  data_df = data.frame(comp = splsda_model$variates$X, component=rep('1', length(splsda_model$variates$X)),
                       group_labels = splsda_model$Y, 
                       shape = shape_vector)
  names(data_df)[which(names(data_df) == 'comp.1')] = 'value'  # Rename column. 
  if (! is.null(shape_vector)) {
    # There is shape information; separate groups by this too. 
    p = ggplot(data_df, aes(x=value, y=component, color=group_labels, shape=shape)) +
      # cex adjusts point spacing within group; dodge.width adjusts spacing between groups. 
      geom_beeswarm(cex = cex_1D, groupOnX = FALSE, size = 1, dodge.width=1.3, aes(color=group_labels, shape=shape))
  } else {
    p = ggplot(data_df, aes(x=value, y=component, group=group_labels)) +
      geom_beeswarm(cex = 1.2, groupOnX = FALSE, size = 1, dodge.width=1.5, aes(color=group_labels, shape=group_labels))
  }
  p = p +
    theme(text = element_text(size=8)) +
    theme(legend.position="none") +
    ylab(y_label) +
    theme(axis.title.y = element_text(size=7)) +
    theme(axis.text.y = element_blank()) +  # Turn off tick label titles.
    xlab(paste('PC1: ', round(100 * splsda_model$explained_variance$X), '% expl. var', sep='')) +
    theme(axis.title.x = element_text(size=7)) +
    theme(axis.text.x = element_blank()) +  # Turn off tick label titles.
    theme(axis.ticks = element_blank())  # Turn off black tick marks.
  if (! is.null(col_per_group)) {
    p = p + scale_color_manual(values=col_per_group)
  }

  if (! is.null(graph_path)) {
    ggsave(graph_path, width=width_mm, height=height_mm, units='mm')
  }
  return(p)
}


#### Gauge the statistical significance of a given real sPLS-DA model's performance by permutation testing. 
# The real sample-class assignments are permuted and the model is trained again. 
# Use the supplied arguments in exactly the same manner as the real model was built on (otherwise it is not a fair comparison)
statisitcal_significance_permutation_test = function(seed=123, repetitions=50,
                                                     feature_table,  # Features as columns, instances as rows.
                                                     select_max_ncomp, 
                                                     select_distance,  # e.g. 'mahalanobis.dist'
                                                     select_test_keepX,  # c(<integer values>)
                                                     select_error_mode,  # 'BER' or 'overall'
                                                     select_validation,  # e.g. 'loo'
                                                     select_logratio = 'CLR',
                                                     select_cpus = 8, 
                                                     data_write_path = NULL,  # Do not include extensions e.g. '.pdf'
                                                     graph_title = '',
                                                     real_model_performance,
                                                     real_model_classes  # vector of classes used for the real model.
                                                     )
{
  set.seed(seed)  # Ensure reproducibility
  repetitions = repetitions
  lowest_errors = rep(NA, repetitions)
  for (repetition in 1:repetitions)
  {
    randomised_categories = factor(sample(real_model_classes))  # Shuffle the class labels assigned to our samples.
    tune_splsda = tune.splsda(X=feature_table, 
                              Y=randomised_categories,
                              ncomp=select_max_ncomp, 
                              logratio=select_logratio, 
                              test.keepX=select_test_keepX,  
                              validation=select_validation,
                              dist=c(select_distance),  
                              measure=select_error_mode,  # Our classes are unbalanced, this is a better measure.  
                              scale=TRUE,
                              near.zero.var=TRUE,
                              progressBar=FALSE,
                              cpus=select_cpus
    )
    lowest_errors[repetition] = min(tune_splsda$error.rate)
  }

  #### Make graph for publication
  # Paper is written in terms of accuracy, not error. 
  real_model_best_accuracy = 100 * (1 - 
            # Gets the last item, which is be the best performance.  
            tail(real_model_performance$error.rate[[select_error_mode]][, paste(select_distance, '.dist', sep='')], 
            n=1))
  highest_accuracies = (1 - lowest_errors) * 100
  
  df = data.frame(accuracies = highest_accuracies)
  p = ggplot(df, aes(accuracies)) + stat_ecdf(geom = "step") + 
    # Represent where the real model's accuacy lay. 
    geom_vline(xintercept=real_model_best_accuracy, color='red') +
    ggtitle(graph_title) +
    xlab('Highest accuracies following tuning (%)') + 
    ylab('Cumulative distribution') +
    xlim(c(0, 100)) + ylim(c(0, 1)) +
    theme(text = element_text(size=8))
  if (! is.null(data_write_path)) {
    ggsave(paste(data_write_path, '.pdf', sep = ''), width=7, height=5, units='cm')
  }
  
  beaten_by_random = sum(sort(highest_accuracies) >= real_model_best_accuracy)
  p_val_resolution = 1 / repetitions
  p_value_upper = (beaten_by_random + 1) * p_val_resolution
  p_value_lower = beaten_by_random * p_val_resolution
  results_text = paste(
    'Highest accuracy from the real model was ', real_model_best_accuracy, '% (using ', select_error_mode, ' error calculations)\n',
    'Random group re-assignments beat the performance of the real data (accuracy) ', beaten_by_random, ' times\n',
    'This gives a p-value in range: ', 
    p_value_lower, ' < p <= ', p_value_upper, 
    ' (', repetitions, ' repetitions, giving p-resolution of ', p_val_resolution, ')\n', sep = '')
  cat(results_text)
  # Write to file. 
  sink(paste(data_write_path, '.txt', sep = ''))
  cat(results_text)  
  sink()
  return(invisible(list(plot=p, highest_accuracies=sort(highest_accuracies), lowest_errors=sort(lowest_errors))))
}


#### Wrapper that handles a select number of components for a splsda model. 
extract_important_features_all_components = function(
  model,  # an splsda model
  model_perf,  # Same `model``, but run through the `perf` function.
  max_components,  # How many components to explore. 
  csv_file_prefix,  # Specify directory and prefix (to identify model and groups explored) for csv files. 
  plot_loadings = TRUE,
  # 1D matrix, OTU/feature indexes as row names, the single column containing
  # a human-readable string (column named "name").
  feature_map = NULL)
{
  all_component_features = list()  # Populate in following loop.
  # Analyse each component in turn. Plot loadings, and output table that also describes stability in selection under LOO.
  for (component in 1:max_components) {
    contributions = plotLoadings(model, comp=component, 
                                 method='median',  # Advised for microbiome data
                                 contrib='max', 
                                 plot=FALSE  # Allows user to extract contribution matrix. Extraction disabled under TRUE. 
    )
    comp_features = extract_important_features(
      model_perf, contributions,
      component=component,
      filename=paste(csv_file_prefix, '_loadingsDim', component, '_selectedFeatures.csv', sep=''),
      feature_map=feature_map)
    all_component_features[[paste('component', component, sep='')]] = comp_features
    if (plot_loadings)
      pl = plotLoadings(model, comp=component, method='median',  # Advised for microbiome data
                        contrib='max')  
  }
  # How many features used in each component, and for each group. 
  for (component in 1:max_components) {
    used_features = all_component_features[[component]][complete.cases(all_component_features[[component]]), ]
    cat('\nComponent', component, 'employed', dim(used_features)[1], 'features\n')
    groups = unique(used_features$GroupContrib)
    for (group in groups) {
      group_features = used_features[used_features$GroupContrib == group, ]
      cat(' ', dim(group_features)[1], 'feature(s) were associated with group', group, '\n')
    }
  }
  
  return(invisible(all_component_features))
}

#### Writes a table, captures features used in training a splsda model, but also queries the performance object 
# that evaluates it to capture how frequently each feature was used in training a model under cross validation. 
extract_important_features = function(splsda_perf, 
                                      # As generated by the plotLoadings function (w plot=FALSE)
                                      feature_component_contributions,  
                                      # 1D matrix, OTU/feature indexes as row names, the single column containing
                                      # a human-readable string (column named "name"). 
                                      feature_map = NULL, 
                                      component, filename)
{
  # How often (in cross validation) was each feature selected in model building?
  # Behaves differently if there is 1 vs >1 items in `splsda_perf$features$stable[[component]]`.
  # Try to handle both cases with the if statement conditional on length. 
  stability = as.data.frame(splsda_perf$features$stable[[component]])  
  if (length(splsda_perf$features$stable[[component]]) == 1)
    colnames(stability) = c('Freq')  # This isn't named correctly for length = 1
  else
    stability$Var1 = NULL  # Delete this column if length > 1
  rownames(stability) = names(splsda_perf$features$stable[[component]])
  
  # Adds new rows (containing NA values) for features used in model building, but not in the final model.
  feature_component_contributions = merge(feature_component_contributions, stability, 
                                          by.x='row.names', by.y='row.names', all=TRUE)
  rownames(feature_component_contributions) = feature_component_contributions$Row.names
  feature_component_contributions$Row.names = NULL
  
  # Rename columns 'importance' to 'loading'
  colnames(feature_component_contributions)[colnames(feature_component_contributions) == 'importance'] = 'loading'
  # feature_component_contributions$freq = freq
  if (! is.null(feature_map)) {
    # Put phylogenetic names in table. 
    feature_component_contributions$human_readable = feature_map[rownames(feature_component_contributions), 'name']
  }
  
  # Sort by frequency in cross validation loops, and then by absolute magnitude of loading. 
  feature_component_contributions = feature_component_contributions[
    order(feature_component_contributions$Freq, # Re-order rows
          abs(feature_component_contributions$loading),
          decreasing=TRUE)
    , ]  # All columns
  
  write.csv(feature_component_contributions, file=filename)
  return(feature_component_contributions)
}

#### Find the taxa (in terms of taxonomic names) that possess a given enzyme. Operates on PICRUSt2 data.
extract_taxa_possessing_enzyme = function(
  enzyme,  # E.g. "EC:3.2.1.26"
  enzyme_function='',  # E.g. "beta-fructofuranosidase"
  # Dataframe. PICRUSt2-generated file (e.g. EC_predicted.tsv). 
  # ASV sequences as rows, features (e.g. EC:3.2.1.26) as cols.
  ec_asv_tab, 
  # Dataframe. DADA2-generated file. First column is ASV sequences. Subsequent columns are:
  # "Kingdom,Phylum,Class,Order,Family,Genus,Species". Maps sequences to human readable taxonomy. 
  asv_taxonomy,
  filename=NULL,  # Where to save CSV
  return_table=FALSE
  )
{
  # Create a mask, extract rownames. Returns the ID, which is rRNA sequence. 
  asv_ids_possessing_ec = rownames(ec_asv_tab[ec_asv_tab[, enzyme] != 0, ])  
  asv_possessing_ec = asv_taxonomy[asv_ids_possessing_ec, ]  # Retrieve corresponding taxonomy. 
  asv_possessing_ec$Sequence_16S = rownames(asv_possessing_ec)  # Sequence as last column. 
  if (! is.null(filename))
    write.csv(asv_possessing_ec, file=filename, quote=FALSE, row.names=FALSE)
  cat('Found', length(asv_ids_possessing_ec), 'capable of producing', enzyme, '-', enzyme_function, 
      'out of', dim(ec_asv_tab)[1],'\n')
  if (return_table)
    return(asv_possessing_ec)
}

#### Use sPLSDA as a feature selection algorithm, but not as a classifier. 
# Allows use for other (e.g. nonlinear) classifiers operating over this data. 
# Requires that a sPLS-DA be tuned before use: select_ncomp and select_keepX need to be provided in advance. 
# In standard sPLS-DA usage through mixOmics, the splsda-transformed data for each fold cannot be easily extracted. 
# Hence, this function will do it.
# Implements its own leave one out cross validation (LOO)-partitioning of the data, trains a sPLS-DA on each fold 
# and writes the transformed data to the file system. 
#
# Data written are:
# 1. The training X (features) and Y (response)
# 2. The test X and Y. Files are named with the sample name being withheld in this fold.
# 3. sPLS-DA transformed training X and test X data. 
# (4. sPLS-DA loadings are also provided, though I am not sure they are useful)
# 5. A 2D (or 1D, problem dependent) ordination of the TEST data in sPLS-DA transformed space. This may not be all the 
#    dimensions, but will help with thinking about the data. 
#
# Requires the splsda_transform_newdata.R file be loaded.
#
# NOTE. I have not seriously attempted any classifiers on the sPLS-DA-transformed data. 
# The distance-based classification that ensues may already be optimal. 
extract_loo_folds_splsda_feature_select = function(
  select_ncomp,  # Integer, # of sPLS-DA components to use. Hence, sPLS-DA training on this data need be done beforehand. 
  select_keepX,  # Vector of integers, one value per component. 
  feature_table,  # Samples as rows, features as columns.
  class_labels,  # Vector of factors, one item per sample.  
  ordinate_folds = FALSE,  # 2D (or 1D) ordination of the sPLS-DA trained on the training data. 
  col_per_group = NULL,  # Vector of named variables: c('group' = '#FAFFFC', ...)
  select_logratio = 'CLR',
  problem_label = '.'
  )
{
  
  num_samples = length(class_labels)
  sample_indices = 1:num_samples
  for (fold in 1:num_samples)  # Leave one out
  {
    train_indices = sample_indices[-fold]  # Leaves 'fold' out of the sequence
    test_indices = sample_indices[fold]
    
    train_X = feature_table[train_indices, ]  # The features / taxa
    train_Y = class_labels[train_indices]  # The response
    
    # As usual, R being awesome and all, the return type of this varies with the number of indices given.
    # If multiple, then a matrix is returned and all works well. However, for a single index, a vector is returned and
    # that breaks downstream code.
    # Hence, need to ensure correctly oriented matrix.
    test_X = feature_table[test_indices, ]  # leave one out, there is only a single sample in each test dataset.
    if (length(test_indices) == 1) {
      test_X = t(as.matrix(test_X))
    }
    
    test_Y = class_labels[test_indices]
    test_sample_name = rownames(feature_table)[fold]
    
    # Make directory structure to store the folds in.
    fold_dir = paste(problem_label, '/loo_folds/', fold, sep='')
    dir.create(fold_dir, recursive=TRUE, showWarnings=FALSE)
    
    # Write the test and training X (predictor, taxa) and Y (response, outcome) data to file system.
    # Requires some messing around with data frames to correctly label samples and columns, as otus_relab is a matrix.
    write.csv(as.data.frame(train_X), file=paste(fold_dir, '/train_X_.csv', sep=''),
              quote=FALSE)
    train_Y_df = data.frame(train_Y, row.names=rownames(train_X), check.names=FALSE)
    colnames(train_Y_df) = c('Response')  # Give name to column. 
    write.csv(train_Y_df,
              file=paste(fold_dir, '/train_Y_.csv', sep=''),
              quote=FALSE)
    
    test_X_df = data.frame(test_X, row.names=c(test_sample_name), check.names=FALSE)
    write.csv(test_X_df,
              file=paste(fold_dir, '/test_X_', test_sample_name, '.csv', sep=''),  quote=FALSE)
    test_Y_df = data.frame(test_Y, row.names=c(test_sample_name))
    colnames(test_Y_df) = c('Response')
    write.csv(test_Y_df, file=paste(fold_dir, '/test_Y_', test_sample_name, '.csv', sep=''),  quote=FALSE)
    
    fold_splsda = splsda(train_X, Y=train_Y, ncomp=select_ncomp, logratio=select_logratio, keepX=select_keepX)
    write.csv(fold_splsda$variates$X, file=paste(fold_dir, '/splsda_transformed-train_X.csv', sep=''))
    if (ordinate_folds){
      p = ordinate(fold_splsda, background = TRUE,
                   group = class_labels,
                   plot_title = paste('sPLS-DA: ', problem_label_human, sep=''), name_samples = FALSE,
                   col_per_group_map = col_per_group,
                   show_legend = FALSE,
                   dimensions_mm = c(50, 40), fontsize = 7, point_size = 1, spartan_plot = TRUE,
                   filename = paste(fold_dir, '/', problem_label, 'splsda-training_data.pdf', sep='')
      )
    }
    
    # Write loadings to file.
    write.csv(fold_splsda$loadings$X, file=paste(fold_dir, '/splsda_transformed-train_loadings.csv', sep=''))
    
    # This handles centering, scaling and performig a CLR, all based on data stored within the splsda object.
    res = splsda_transform(object=fold_splsda, newdata=test_X)
    test_transformed = res$newdata.transformed$X
    # Set names consistent with training data.
    colnames(test_transformed) = colnames(fold_splsda$variates$X)
    if (length(test_sample_name) == 1)
      rownames(test_transformed) = c(test_sample_name)
    write.csv(test_transformed, file=paste(fold_dir, '/splsda_transformed-test_X-', test_sample_name, '.csv', sep=''))
  }
}


standard_full_splsda_pipeline = function(
  feature_table,  # DataFrame or matrix. Features as columns, instances as rows. Typically relative abundance. 
  class_labels,  # Vector of factors. One item per instance, in the same order as feature_table. 
  col_per_group,  # Named vector, e.g. c('TRF' = '#F1234FF', ...)
  problem_label_human,  # Human-readable title. 
  problem_label,  # Machine-readable title. (avoid spaces, and capitalisation)
  feature_map = NULL,
  perform_permutation_analysis = FALSE,  # Generate p-values? Extremely computationally expensive. 
  select_max_ncomp,  # How many components to attempt when tuning the sPLS-DA. 
  select_validation = 'loo',
  select_test_keepX = c(1, 2, 3, 4, seq(5, dim(feature_table)[2], 5)),
  select_distance = c('mahalanobis'),
  select_error_mode = 'BER',  # {'BER', 'overall'}
  select_scale = TRUE,
  select_near_zero_var = TRUE,
  logratio_transform = 'CLR',  # {'CLR', 'none'}
  cpus = 8,
  extract_loo_folds = FALSE,
  cex_1D = 5.0,  # Spacing between groups on a 1D swarm plot. 
  ordination_dim_2d = c(50, 40),  # Physical dimensions of 2D ordinations pots. 
  sample_plot_characters = 16  # Single char, else vector of n=len(samples).
  )
{
  dir.create(problem_label, showWarnings = FALSE)
  if (max(select_test_keepX) > dim(feature_table)[2]){
    cat('Warning! Selected to test more features per component than there are features. Curtailing select_test_keepX.')
    select_test_keepX = select_test_keepX[select_test_keepX < dim(feature_table)[2]]
  }
  
  tune_splsda = tune.splsda(X = feature_table, 
                            Y = class_labels,  # Labels
                            ncomp = select_max_ncomp, 
                            logratio = logratio_transform, 
                            test.keepX = select_test_keepX,  
                            validation = select_validation,
                            dist = c(select_distance),  
                            measure = select_error_mode,
                            scale = select_scale,
                            near.zero.var = select_near_zero_var,
                            progressBar = FALSE,
                            cpus = cpus
                            )
  cat('Tuned parameters yield the following performance for each additional component:\n')
  plot(tune_splsda)
  print(tune_splsda$error.rate.class) 
  
  select_ncomp = which.min(colMins(tune_splsda$error.rate))
  cat('Based on optimal performance, selecting', select_ncomp, 'principal component(s) in building the sPLS-DA.\n')
  select_keepX = tune_splsda$choice.keepX[1:select_ncomp]  # Number of features per component. 
  cat('Features per component:\n')
  print(select_keepX)
  tuned_splsda = splsda(X = feature_table, Y = class_labels, 
                        ncomp = select_ncomp, keepX = select_keepX,  # From tuning, above. 
                        logratio = logratio_transform, scale = select_scale,  # Microbiome data is compositional. 
                        near.zero.var = select_near_zero_var)  # Set to true because microbiome data contains many zeros. 
  
  p = ordinate(tuned_splsda, background = TRUE,
               group = class_labels,
               plot_title = paste('sPLS-DA: ', problem_label_human, sep=''), name_samples = FALSE,
               col_per_group_map = col_per_group,
               show_legend = FALSE,
               dimensions_mm = ordination_dim_2d, 
               sample_plot_characters = sample_plot_characters,
               fontsize = 7, point_size = 1, spartan_plot = TRUE,
               filename = paste(problem_label, '/', problem_label, '-splsda.pdf', sep=''),
               cex_1D = cex_1D)
  p  # Display graph.
  tuned_splsda_perf = perf(tuned_splsda, validation="loo", progressBar=FALSE, auc=FALSE)
  cat('sPLS-DA using tuned parameters gives the following performance:\n')
  print(tuned_splsda_perf$error.rate)
  print(tuned_splsda_perf$error.rate.class)
  # Write same results to file system. 
  sink(paste(problem_label, '/', problem_label, '-splsda-performance.txt', sep=''))
  print(tuned_splsda_perf$error.rate)
  print(tuned_splsda_perf$error.rate.class)
  sink()
  
  # Analyse the features informing each component. 
  extract_important_features_all_components(
    model = tuned_splsda, model_perf = tuned_splsda_perf, max_components = select_ncomp, plot_loadings = FALSE, 
    feature_map = feature_map,  # Human reaadble feature names. 
    csv_file_prefix = paste(problem_label, '/', problem_label, '-splsda', sep = '')
  )
  
  # Permutation based statistics. 
  if (perform_permutation_analysis) {
    permutation_analysis_results = statisitcal_significance_permutation_test(
      feature_table = feature_table, select_max_ncomp = select_max_ncomp, select_distance = select_distance,
      select_test_keepX = select_test_keepX, select_error_mode = select_error_mode, select_validation = select_validation,
      select_logratio = logratio_transform,
      data_write_path = paste(problem_label, '/', problem_label, '-splsda_RANDOMISED', sep = ''), 
      graph_title = paste(problem_label_human, ' (', length(unique(class_labels)),' groups)', sep = ''),
      real_model_performance = tuned_splsda_perf, real_model_classes = class_labels, repetitions = 50)
  }
  
  # Extract sPLS-DA-transformed training and test data generated from within a LOO loop. 
  # Enables other external classifiers to be used. 
  if (extract_loo_folds) {
    extract_loo_folds_splsda_feature_select(
      select_ncomp = select_ncomp,
      select_keepX = select_keepX,  # Vector of integers, one value per component. 
      feature_table = feature_table,  # Samples as rows, features as columns.
      class_labels = class_labels,  # Vector of factors, one item per sample.  
      select_logratio = logratio_transform,
      ordinate_folds = TRUE,  # 2D (or 1D) ordination of the sPLS-DA trained on the training data. 
      col_per_group = col_per_group,
      problem_label = problem_label)
  }
}
