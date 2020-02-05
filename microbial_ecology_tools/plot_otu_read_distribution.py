# Analysis reveals the relationship between OTU read number and reads belonging to OTUs by their
# read number. It is intended to help in choosing thresholds for filtering OTU tables. 
#
# A great many OTUs contain very few reads, and a small minority of reads belong to OTUs
# with a low read number. By graphing cumulative distribtions of each you can observe how 
# selecting a given reads/OTU threshold will reduce the total number of reads and OTUs. 
#
# A plot showing boxplots of reads/OTU for each sample is drawn. This highlights any samples 
# that will be unduely affected by a given threshold. Ideally, all samples should contain a 
# similar distribution of reads/OTU. 
#
# $> python plot_otu_abundances.py -i otu_table.biom -o output/directory
#
# Mark Read, 2019. 

# from biom import load_table
import sys
import os
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

def ecdf(observations):
    """
    Empirical Cumulative Distribution Function.
    arguments:
      :observations : list of observations, in any order.

    Returns a list of tuples: [ (<value>,<proportion of observations less than or equal to value>) ]
    Values will be given in ascending order.
    """
    s = sorted(observations)    # non-destructive, preserves original order of observations.
    n = len(s)                    # number of observations.
    cdf = []
    for i, sample in enumerate(s):
        proportion = float(i+1)/n
        tup = (sample, proportion)
        cdf.append(tup)
    return cdf


def plotCDF(dist1, xlabel, filename, name1=None, dist2=None, name2=None, title=None, yrange=None, ylog=False):
    """
    Calculates cumulative distribution functions for the two supplied empirical distributions, and plots the two CDFs
    atop one another.
    """
    plt.clf()
    cdf1 = ecdf(dist1)
    D1X = [record[0] for record in cdf1]
    D1Y = [record[1] for record in cdf1]
    p1, = plt.plot(D1X, D1Y, 'b')
    if dist2:
        cdf2 = ecdf(dist2)
        D2X = [record[0] for record in cdf2]
        D2Y = [record[1] for record in cdf2]
        p2, = plt.plot(D2X, D2Y, 'r')
        #plt.legend([,loc='upper center', bbox_to_anchor=(0.5, 1.1),ncol=5,prop={'size':14})
        plt.legend([p1, p2], [name1, name2], loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.1), prop={'size':16})
    plt.xlabel(xlabel)
    plt.ylabel('Proportion')
    if ylog:
        plt.yscale('log')
    plt.xscale('log')
    if yrange:
        plt.ylim(yrange)
    if title:
        plt.title(title)
    plt.gca().grid(True)   # turn on grid lines.
    font = {'size'   : 18}
    plt.rc('font', **font)
    plt.savefig(filename, dpi=300)


def plotCDFs(dist1, dist2, d1name, d2name, xlabel, filename, title=None, dist3=None, d3name=None):
    """
    Calculates cumulative distribution functions for the two supplied empirical distributions, and plots the two CDFs
    atop one another.
    """
    plt.clf()
    cdf1 = stats.ecdf(dist1)
    cdf2 = stats.ecdf(dist2)
    D1X = [record[0] for record in cdf1]
    D1Y = [record[1] for record in cdf1]
    D2X = [record[0] for record in cdf2]
    D2Y = [record[1] for record in cdf2]
    p1, = plt.plot(D1X, D1Y)
    p2, = plt.plot(D2X, D2Y)
    plt.legend([p1, p2], [d1name, d2name], loc='lower right')
    if dist3:
        cdf3 = stats.ecdf(dist3)
        D3X = [record[0] for record in cdf3]
        D3Y = [record[1] for record in cdf3]
        p3, = plt.plot(D3X, D3Y)
        plt.legend([p1, p2, p3], [d1name, d2name, d3name], loc='lower right')
    plt.xlabel(xlabel)
    plt.ylabel('Proportion')
    if title:
        plt.title(title)
    plt.gca().grid(True)   # turn on grid lines.
    font = {'size'   : 18}
    plt.rc('font', **font)
    plt.savefig(filename, dpi=300)


def main():
    # Read command line arguments. 
    inputFile = sys.argv[sys.argv.index('-i') + 1]  # Input file, an otu BIOM file. 
    outputDir = sys.argv[sys.argv.index('-o') + 1]  # Output directory

    # Create output directory if it does not already exist.
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    
    table = load_table(inputFile)
    observations = table.ids(axis='observation')  # equates to OTUs
    samples_raw = table.ids(axis='sample') 
    
    otu_read_count = dict()  # dict[otu_id] = number of reads otu_id has.
    for otu_id in observations:  # Cycle through each OTU's id. 
        # Returns an np.ndarray (1d=list-like), number of reads for each sample against this OTU.
        reads_in_sample = table.data(otu_id, axis='observation', dense=True)
        # Sum across all the samples to get a count for this OTU.
        otu_read_count[otu_id] = reads_in_sample.sum()

    # read_numbers = list(otu_read_count.values())
    # plot CDF of OTUs by read number. 
    otu_graph_path = outputDir + '/otus_by_reads_cdf.png'
    plotCDF(dist1=otu_read_count.values(), yrange=[0.00, 1.0], xlabel='Reads', filename=otu_graph_path, title='CDF of OTUs by reads')

    # Create list, with an item per read. Register the size of the OTU that each read belongs to. 
    samples_by_otu_read_number = []
    for otu_reads in otu_read_count.values():  # Number of reads each OTU registered.
        for read in list(range(int(otu_reads))):  # Iterate through each of the OTU's reads
            # Log the size of the OTU (in terms of reads) that this current read belongs to. 
            samples_by_otu_read_number.append(otu_reads)

    # plot CDF of number of reads for the OTU that each sample belongs to. 
    sample_graph_path = outputDir + '/read_by_otu_reads_cdf.png'
    plotCDF(dist1=samples_by_otu_read_number, yrange=[0.00, 1.0], xlabel='OTU Reads', filename=sample_graph_path, title='CDF of reads by constituent OTU read number')    

    # combine the above two CDFs into 1 plot. 
    both_graph_path = outputDir + '/otus_reads_cdf.png'
    plotCDF(dist1=samples_by_otu_read_number, name1='reads by otu',
        dist2=otu_read_count.values(), name2='otus',
        yrange=[0.00, 1.0], xlabel='OTU Reads', filename=both_graph_path)    


    # # calculate distribution of OTU read numbers for each sample.
    # samples = dict()
    # for reads in otu_reads.values():
    #     n = len(reads)
    #     # don't log this OTU's read number for the same sample twice. 
    #     # this stores the samples already encountered for this OTU. 
    #     samples_seen = []    
    #     for read in reads:                        
    #         sampleID = read[:read.find('_')]   # extract text to the left of the first underscore. 
    #         # only log if this a read for this sample has not already been encountered in this OTU. 
    #         if sampleID not in samples_seen:
    #             samples_seen.append(sampleID)
    #             # initialize empty list for first encounter of given sample. 
    #             if sampleID not in samples.keys():
    #                 samples[sampleID] = []   
    #             samples[sampleID].append(n)
        
    # keys = samples.keys()
    # keys.sort()
    # data = [samples[key] for key in keys]
    # print keys
    # plt.clf()
    # plt.boxplot(data)
    # plt.xticks(range(1,len(keys)+1), keys)
    # plt.yscale('log')
    # plt.ylabel('Reads/OTU')    
    # fig = plt.gcf()
    # fig.set_size_inches(18.5, 10.5)
    # plt.savefig(outputDir + '/otus_reads_by_sample', dpi=300)

if __name__ == "__main__": main()






