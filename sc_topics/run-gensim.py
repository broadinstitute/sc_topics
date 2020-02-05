#!/bin/usr/python3

import argparse
import sys
import datetime

import scanpy as sc
import scipy.io
import numpy as np
import pandas as pd
import gensim
from gensim import corpora, models, similarities


## define main function for running analysis
def main(adata_file, subset_col, subset_val, topics, project_directory, output):
    
    ## read in AnnData
    adata = sc.read(adata_file)
    
    if subset_col != 'all':
        adata_use = adata[adata.obs[subset_col] == subset_val].copy()
    else:
        adata_use = adata.copy()
        
        
        
    ## Extract Sparse gbm 
    mat = adata_use.X.transpose(copy=True)
    geneids = adata_use.var_names


    ## Create Vocab list of genes
    id_list = geneids.tolist()
    out = [[]]
    for i in id_list: out.append([i])
    ## Turn into dictionary for use in model
    dictionary = corpora.Dictionary(out)
    ## Convert gbm to a corpus format for model
    corpus = gensim.matutils.Sparse2Corpus(mat)
    corpora.MmCorpus.serialize(project_directory + '/corpus/' + output + '_corpus.mm', corpus)


    ## Latent Dirichlet Allocation ####
    print('Running LDA with ', str(topics) + ' topics')
    lda = models.LdaModel(corpus, id2word=dictionary, num_topics=topics)
    lda.save(project_directory + '/models/' + output + '_lda_model')

    
    ## Cell Topics
    cell_scores_lda = lda.get_document_topics(corpus)
    sp_mat_lda = gensim.matutils.corpus2csc(cell_scores_lda)
    ## topic by cell/documents
    cell_topics_lda = pd.DataFrame(sp_mat_lda.todense()).T
    cell_topics_lda['index'] = adata_use.obs.index.tolist()
    cell_topics_lda.set_index('index', inplace=True)
    cell_topics_lda.to_csv(project_directory + '/outputs/' + output + '_cell_scores_lda.csv')

    
    ## every topic by every gene
    topic_scores_lda = pd.DataFrame(lda.get_topics()).T
    topic_scores_lda.to_csv(project_directory + '/outputs/' + output + '_gene_topics_lda.csv')
    
    
    ## Hierarchical Dirchlet Process ####
    print('Running HDP')
    hdp = models.HdpModel(corpus, id2word=dictionary)
    hdp.save(project_directory + '/models/' + output + '_hdp_model')
    
    
    ## Cell Topics
    cell_scores_hdp = hdp.get_document_topics(corpus)
    sp_mat_hdp = gensim.matutils.corpus2csc(cell_scores_hdp)
    ## topic by cell/documents
    cell_topics_hdp = pd.DataFrame(sp_mat_hdp.todense()).T
    cell_topics_hdp.to_csv(project_directory + '/outputs/'+ output + '_cell_scores_hdp.csv')
    cell_topics_hdp['index'] = adata_use.obs.index.tolist()
    cell_topics_hdp.set_index('index', inplace=True)
    
    ## every topic by every gene
    topic_scores_hdp = pd.DataFrame(hdp.get_topics()).T
    topic_scores_hdp.to_csv(project_directory + '/outputs/' + output + '_gene_topics_hdp.csv')

    

    ## adata_file, subset_col, subset_val, topics, project_directory, output
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', type=str, help='Anndata file to use?')
    parser.add_argument('--subset_col', type=str, default='all', help='What subset column to use?')
    parser.add_argument('--subset_val', type=str, default='all', help='What value to subet to?')
    parser.add_argument('--topics', type=int, default=20, help='How many topics for LDA?')
    parser.add_argument('--project_directory', type=str, help='Project directory?')
    parser.add_argument('--output', type=str, help='What output to prefix to use?')
    args = parser.parse_args()

    main(adata_file = args.file,
         subset_col = args.subset_col,
         subset_val = args.subset_val,
         topics = args.topics,
         project_directory = args.project_directory,
         output = args.output)
