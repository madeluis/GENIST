FILES:
- MAIN FILE:
	- GENIST.m:  GENIST is an algorithm to infer gene regulatory networks from spatial and temporal datasets. The spatial dataset (or any data that can provide information about coexpression) is used by the first step of the algorithm to perform clustering and separate the genes in the network in smaller coexpressed groups. The temporal dataset is used by the second step of the algorithm to infer regulations among the genes, based on Bayesian networks.

- SECONDARY FILES (Called by GENIST):
	- GRN_inference.m
	- GRN_preprocessing.m
	- GRN_preprocessing_TFs_to_genes.m
	- matrix_to_cytoscape_table.m
- OPTIONAL FILES (Arabidopsis specific): 
	- Arabidopsis_TFs.xlsx
	- Locus_Primary_Gene_Symbol_2013.xlsx. 


Running GENIST:

 Inference of Gene Regulatory Networks (GRNs) from spatio temporal data
%
% GENIST(genes_file,time_course_file,clustering_data_file)
% infers regulations among the genes provided in genes_file. The
% regulations are calculated by applying a first step of clustering that
% groups the genes based on their coexpression in the clustering_file data
% (which should be preferently a spatial dataset) followed by a second step
% of a Bayesian Network inference, which uses the data in time_course_file.
    % genes_file must be an excel file containing the list of genes to be
    % included in the network as a column.
    %
    % time_course_file must be an excel file containing expression data of
    % the genes across different time or developmental time points. Each
    % row must correspond to 1 gene. The file can contain any number of
    % genes, and GENIST selects the genes specified in genes_file.
    %
    % clustering_data_file must be an excel file containing expression data
    % of the genes across different spacial confinements. Each
    % row must correspond to 1 gene. The file can contain any number of
    % genes, and GENIST selects the genes specified in genes_file.
    % The names of the genes in all .xlsx files must be consistent. 
    % 
%
%    
% GENIST(genes_file,time_course_file,clustering_data_file,time_lapse,is_load_new_data)
% sets the values of time_lapse and is_load_new_data.
    % time_lapse indicates the lapse between the activation of genes that
    % the algorithm will consider to select the regulators of a gene.
    % If time_lapse is 1, GENIST assumes that genes that are active in one
    % time point can regulate other genes 1 time point later.
    % If time_lapse is 0, GENIST assumes that genes that are active in one
    % time point can regulate other genes in the same time point. Default
    % is 1.
    %
    % If is_load_new_data is true, all data is loaded and saved into a Matlab.
    % If is_load_new_data is false, the data last saved is used to run GENIST.
    % Default is true.
%
%
% GENIST(genes_file,time_course_file,clustering_data_file,time_lapse,is_load_new_data,TF_file,symbol_file)
% sets a list of Transcription Factors (TFs) and symbols of the genes.
    % TF_file must be an excel file containing a list of TFs as a column
    % vector. If TF_file is provided, only the genes (from genes_gile) that
    % are also included in the TF_file will be considered as potential
    % regulators in the network. 
    %
    % symbol_file must be an excel file containing a list of genes in the
    % first column, and a corresponding symbol for each gene in the second
    % column. When symbol_file is provided, the genes in the network for
    % which a symbol was provided are represented by their symbol.
    % The names of the genes in all .xlsx files must be consistent. 
    %
% GENIST returns a plot of the resulting network in MATLAB, as well as a
% file that is saved as 'cityscape_table' in the working directory. This
% file contains a table with all the regulations of the network. This file
% can  be imported into cytocape to generate a plot of the network.
%
%
% Example 1:
% GENIST('gene_list.xlsx','DevTime','spatial_test',1,true,'Arabidopsis_TFs','Locus_Primary_Gene_Symbol_2013')
% Example 2:
% GENIST('gene_list.xlsx','DevTime','spatial_test',1,false)
% Example 3:
% GENIST('gene_list.xlsx','DevTime','spatial_test')
% Example 4:
% GENIST('gene_list.xlsx','DevTime','spatial_test',1,true,[],'Locus_Primary_Gene_Symbol_2013')
 
%
% M. Angels de Luis Balaguer
% Postdoctoral Research Scholar
% North Carolina State University
% 2016
