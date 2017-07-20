FILES:
- MAIN FILE:
	- GENIST2_0.m:  GENIST is an algorithm to infer gene regulatory networks from spatial and temporal datasets. The spatial dataset (or any data that can provide information about coexpression) is used by the first step of the algorithm to perform clustering and separate the genes in the network in smaller coexpressed groups. The temporal dataset is used by the second step of the algorithm to infer regulations among the genes, based on Bayesian networks.

- SECONDARY FILES (Called by GENIST):
	- GRN_inference2_0.m
	- GRN_preprocessing2_0.m
	- matrix_to_cytoscape_table2_0.m
- OPTIONAL FILES (Arabidopsis specific): 
	- Arabidopsis_TFs.xlsx
	- Locus_Primary_Gene_Symbol_2013.xlsx. 


Running GENIST:

% Inference of Gene Regulatory Networks (GRNs) from spatio temporal data
%
%
% GENIST2_0(genes_file,time_course_file,clustering_data_file)
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
% GENIST2_0(genes_file,time_course_file,clustering_data_file,time_lapse,is_load_new_data)
% sets the values of time_lapse and is_load_new_data.
% time_lapse indicates the lapse between the activation of genes that
% the algorithm will consider to select the regulators of a gene.
% If time_lapse is 1, GENIST assumes that genes that are active in one
% time point can regulate other genes in the following time point.
% If time_lapse is 0, GENIST assumes that genes that are active in one
% time point can regulate other genes in the same time point. 
% If time_lapse is [0,1], GENIST assumes that genes that are active in one
% time point can regulate other genes in the same time point or in the
% following time point. Default is 1.
%
% If is_load_new_data is true, all data is loaded and saved into a Matlab.
% If is_load_new_data is false, the data last saved is used to run GENIST.
% Default is true.
%
%
% GENIST2_0(genes_file,time_course_file,clustering_data_file,time_lapse,is_load_new_data,TF_file,symbol_file)
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
%
% GENIST2_0(genes_file,time_course_file,clustering_data_file,time_lapse,is_load_new_data,TF_file,symbol_file,is_reg_fc_th,is_reg_time_percent,n_levels,is_low_conn)
% sets the values of is_reg_fc_th, is_reg_time_percent, n_levels, and 
% is_low_conn.
% is_reg_fc_th is the minimum (threshold) fold-change expression that a 
% gene must experience between two consecutive time points for it to be 
% considered a potential regulator of another gene that is changing its
% expression (also by at least a is_reg_fc_th fold-change) in the same
% (time_lapse == 0) or the next (time_lapse == 1) time point. Default is
% 1.3.
%
% is_reg_time_percent is the minimun percentage of the total time points
% that a potential regulator and a target show simultaneous  
% (if time_lapse == 0) or consecutive (if time_lapse == 1) changes of
% expression for a gene to be considered a potential regulator of the
% target gene. Default is 0.3.  
%
% n_levels is the number of levels in which the time course expression data
% will be discretized for calculating the probabilities that each potential
% regulator regulates each target. Must be an integer. Default is 2.
% Minimum allowed value is 2. Recommended values 2 or 3. 
%
% is_low_conn sets the bottom percentage of regulations of the final
% network that will be condirered too low to be displayed (will be
% discarded). Default is 0.2.
%
%
% GENIST2_0 returns a plot of the resulting network in MATLAB, as well as a
% file that is saved as 'cityscape_table' in the working directory. This
% file contains a table with all the regulations of the network. This file
% can be imported into cytocape to generate a plot of the network.
%
%
% Example 1:
% GENIST2_0('gene_list.xlsx','DevTime','spatial_test',1,true,'Arabidopsis_TFs','Locus_Primary_Gene_Symbol_2013')
% Example 2:
% GENIST2_0('gene_list.xlsx','DevTime','spatial_test',1,false)
% Example 3:
% GENIST2_0('gene_list.xlsx','DevTime','spatial_test')
% Example 4:
% GENIST2_0('gene_list.xlsx','DevTime','spatial_test',1,true,[],'Locus_Primary_Gene_Symbol_2013')
% Example 5:
% GENIST2_0('gene_list.xlsx','DevTime','spatial_test',[0,1],true,[],'Locus_Primary_Gene_Symbol_2013',[],.1,[],0.4)
% Example 6:
% genes = 'gene_list.xlsx';
% time_course = 'DevTime.xlsx';
% clustering_data = '../spatial_test.xlsx';
% TF_list = 'Arabidopsis_TFs.xlsx';
% symbols = 'Locus_Primary_Gene_Symbol_2013.xlsx';
% GENIST2_0(genes,time_course,clustering_data,[0,1],true,TF_list,symbols,1.5,.1,3,0.4)

%
% M. Angels de Luis Balaguer
% Postdoctoral Research Scholar
% North Carolina State University
% 2016
