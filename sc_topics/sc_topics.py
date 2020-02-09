from typing import Optional, Union

import numpy as np
import pandas as pd

import scanpy as sc
from anndata import AnnData

from numpy.random.mtrand import RandomState
from sklearn.utils import check_random_state, check_array
from sklearn.decomposition import NMF

# from ._utils import get_init_pos_from_paga, _choose_representation
# from .. import logging as logg
# from .._settings import settings
# from .._compat import Literal

import gensim
from gensim import corpora, models, similarities



def topics(
    adata: AnnData,
    model: str,
    num_topics: int,
    random_state: int,
    update_every: int, # only for ldamodel
    chunksize: int,
    passes: int,
    alpha: str, # 'auto',
    per_word_topics: bool,
    copy: bool    
) -> Optional[AnnData]:
    """\
    Topic modeling description and interface description will go here!


    Parameters
    ----------
    adata
        Annotated data matrix.
    model
        Topic model to use.
    distributed 
        (bool, optional) – Whether distributed computing should be used to accelerate training.
    num_topics
        (int, optional) – The number of requested latent topics to be extracted from the training corpus.
    random_state
        Random seed to use for reproducibility. 
        If `int`, `random_state` is the seed used by the random number generator;
        If `RandomState`, `random_state` is the random number generator;
        If `None`, the random number generator is the `RandomState` instance used
        by `np.random`.
    update_every
        (int, optional) – Number of documents to be iterated through for each update. Set to 0 for batch learning, > 1 for online iterative learning.
    chunksize
        (int, optional) – Number of documents to be used in each training chunk.
    passes
        (int, optional) – Number of passes through the corpus during training.
    alpha 
        ({numpy.ndarray, str}, optional) – Can be set to an 1D array of length equal to the number of expected topics that expresses our a-priori belief for the each topics’ probability. Alternatively default prior selecting strategies can be employed by supplying a string:
            - ’asymmetric’: Uses a fixed normalized asymmetric prior of 1.0 / topicno.
            - ’auto’: Learns an asymmetric prior from the corpus (not available if distributed==True).
    eta 
        ({float, np.array, str}, optional) – A-priori belief on word probability, this can be:
            - scalar for a symmetric prior over topic/word probability,
            - vector of length num_words to denote an asymmetric user defined probability for each word,
            - matrix of shape (num_topics, num_words) to assign a probability for each word-topic combination,
            - the string ‘auto’ to learn the asymmetric prior from the data.
    decay 
        (float, optional) – A number between (0.5, 1] to weight what percentage of the previous lambda value is forgotten when each new document is examined. Corresponds to Kappa from Matthew D. Hoffman, David M. Blei, Francis Bach: “Online Learning for Latent Dirichlet Allocation NIPS’10”.
    per_word_topics
        (bool) – If True, the model also computes a list of topics, sorted in descending order of most likely topics for each word, along with their phi values multiplied by the feature length (i.e. word count).
    copy
        Return a copy instead of writing to adata.



## https://radimrehurek.com/gensim/models/ldamodel.html#gensim.models.ldamodel.LdaModel  #######
offset (float, optional) –

Hyper-parameter that controls how much we will slow down the first steps the first few iterations. Corresponds to Tau_0 from Matthew D. Hoffman, David M. Blei, Francis Bach: “Online Learning for Latent Dirichlet Allocation NIPS’10”.

eval_every (int, optional) – Log perplexity is estimated every that many updates. Setting this to one slows down training by ~2x.

iterations (int, optional) – Maximum number of iterations through the corpus when inferring the topic distribution of a corpus.

gamma_threshold (float, optional) – Minimum change in the value of the gamma parameters to continue iterating.

minimum_probability (float, optional) – Topics with a probability lower than this threshold will be filtered out.

random_state ({np.random.RandomState, int}, optional) – Either a randomState object or a seed to generate one. Useful for reproducibility.

ns_conf (dict of (str, object), optional) – Key word parameters propagated to gensim.utils.getNS() to get a Pyro4 Nameserved. Only used if distributed is set to True.

minimum_phi_value (float, optional) – if per_word_topics is True, this represents a lower bound on the term probabilities.

per_word_topics 

callbacks (list of Callback) – Metric callbacks to log and visualize evaluation metrics of the model during training.

dtype ({numpy.float16, numpy.float32, numpy.float64}, optional) – Data-type to use during calculations inside model. All inputs are also converted.
########################################################################################################

    
    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.
    **X_topics** : `adata.obsm` field
        cell topic scores of data.
    """

    adata_use = adata.copy() if copy else adata

    ## Extract Sparse gbm 
    mat = adata_use.X.transpose()  # copy=True)
    geneids = adata_use.var_names


    ## Create Vocab list of genes
    id_list = geneids.tolist()
    out = [[]]
    for i in id_list: out.append([i])
    ## Turn into dictionary for use in model
    dictionary = corpora.Dictionary(out)
    ## Convert gbm to a corpus format for model
    corpus = gensim.matutils.Sparse2Corpus(mat)
    # corpus = gensim.matutils.Dense2Corpus(mat)
    ## corpora.MmCorpus.serialize(project_directory + '/corpus/' + output + '_corpus.mm', corpus)


    if model == 'nmf':
        nmf = NMF(n_components=num_topics, init='random', random_state=random_state, alpha=decay)
        W = nmf.fit_transform(mat)
        H = nmf.components_

        ## topic by cell/documents
        cell_topics = pd.DataFrame(X)
        cell_topics['index'] = adata.obs.index.tolist()
        cell_topics.set_index('index', inplace=True)

        ## every topic by every gene
        topic_scores = pd.DataFrame(H)


    if model == 'lda':
        ## Latent Dirichlet Allocation ####
        print('Running LDA with ', str(num_topics) + ' topics')
        lda = models.LdaModel(corpus=corpus,
                            id2word=dictionary,
                            num_topics=num_topics,
                            random_state=random_state,
                            update_every=update_every,
                            chunksize=chunksize,
                            passes=passes,
                            alpha='auto',
                            per_word_topics=True)

        ## Cell Topics
        cell_scores = lda.get_document_topics(corpus)
        cell_scores_mat = gensim.matutils.corpus2dense(cell_scores, num_terms=num_topics)

        ## topic by cell/documents
        cell_topics = pd.DataFrame(cell_scores_mat.T)
        cell_topics['index'] = adata.obs.index.tolist()
        cell_topics.set_index('index', inplace=True)

        ## every topic by every gene
        topic_scores = pd.DataFrame(lda.get_topics()).T
        topic_scores['index'] = adata.var_names
        topic_scores.set_index('index', inplace=True)
    

    # if model == "hdp":
    #     ## Hierarchical Dirchlet Process ####
    #     print('Running HDP')
    #     hdp = models.HdpModel(corpus=corpus, 
    #                         id2word=dictionary,
    #                         ## TODO #######
    #                         # max_chunks=None, 
    #                         # max_time=None, 
    #                         # chunksize=256, 
    #                         # kappa=1.0, 
    #                         # tau=64.0, 
    #                         # K=15, 
    #                         # T=150, 
    #                         # alpha=1, 
    #                         # gamma=1, 
    #                         # eta=0.01, 
    #                         # scale=1.0, 
    #                         # var_converge=0.0001, 
    #                         # outputdir=None, 
    #                         # random_state=None)

    #     ## Cell Topics
    #     cell_scores = hdp.get_document_topics(corpus)
    #     sp_mat = gensim.matutils.corpus2csc(cell_scores)

    #     ## topic by cell/documents
    #     cell_topics = pd.DataFrame(sp_mat.todense()).T
    #     cell_topics['index'] = adata.obs.index.tolist()
    #     cell_topics.set_index('index', inplace=True)

    #     ## every topic by every gene
    #     topic_scores = pd.DataFrame(hdp.get_topics()).T


    adata_use.obs = adata_use.obs.join(cell_topics)
    adata_use.var = adata_use.var.join(topic_scores)


    return adata_use # if copy else None

#### UMAP Code ##################################
    # if 'neighbors' not in adata.uns:
    #     raise ValueError(
    #         'Did not find \'neighbors/connectivities\'. Run `sc.pp.neighbors` first.')
    # start = logg.info('computing UMAP')
    # if ('params' not in adata.uns['neighbors']
    #     or adata.uns['neighbors']['params']['method'] != 'umap'):
    #     logg.warning('neighbors/connectivities have not been computed using umap')
    # from umap.umap_ import find_ab_params, simplicial_set_embedding
    # if a is None or b is None:
    #     a, b = find_ab_params(spread, min_dist)
    # else:
    #     a = a
    #     b = b
    # adata.uns['umap'] = {'params':{'a': a, 'b': b}}
    # if isinstance(init_pos, str) and init_pos in adata.obsm.keys():
    #     init_coords = adata.obsm[init_pos]
    # elif isinstance(init_pos, str) and init_pos == 'paga':
    #     init_coords = get_init_pos_from_paga(adata, random_state=random_state)
    # else:
    #     init_coords = init_pos  # Let umap handle it
    # if hasattr(init_coords, "dtype"):
    #     init_coords = check_array(init_coords, dtype=np.float32, accept_sparse=False)

    # random_state = check_random_state(random_state)
    # neigh_params = adata.uns['neighbors']['params']
    # X = _choose_representation(
    #     adata, neigh_params.get('use_rep', None), neigh_params.get('n_pcs', None), silent=True)
    # if method == 'umap':
    #     # the data matrix X is really only used for determining the number of connected components
    #     # for the init condition in the UMAP embedding
    #     n_epochs = 0 if maxiter is None else maxiter
    #     X_umap = simplicial_set_embedding(
    #         X,
    #         adata.uns['neighbors']['connectivities'].tocoo(),
    #         n_components,
    #         alpha,
    #         a,
    #         b,
    #         gamma,
    #         negative_sample_rate,
    #         n_epochs,
    #         init_coords,
    #         random_state,
    #         neigh_params.get('metric', 'euclidean'),
    #         neigh_params.get('metric_kwds', {}),
    #         verbose=settings.verbosity > 3,
    #     )
    # elif method == 'rapids':
    #     metric = neigh_params.get('metric', 'euclidean')
    #     if metric != 'euclidean':
    #         raise ValueError(
    #             f'`sc.pp.neighbors` was called with `metric` {metric!r}, '
    #             "but umap `method` 'rapids' only supports the 'euclidean' metric."
    #         )
    #     from cuml import UMAP
    #     n_neighbors = adata.uns['neighbors']['params']['n_neighbors']
    #     n_epochs = 500 if maxiter is None else maxiter # 0 is not a valid value for rapids, unlike original umap
    #     X_contiguous = np.ascontiguousarray(X, dtype=np.float32)
    #     umap = UMAP(
    #         n_neighbors=n_neighbors,
    #         n_components=n_components,
    #         n_epochs=n_epochs,
    #         learning_rate=alpha,
    #         init=init_pos,
    #         min_dist=min_dist,
    #         spread=spread,
    #         negative_sample_rate=negative_sample_rate,
    #         a=a,
    #         b=b,
    #         verbose=settings.verbosity > 3,
    #     )
    #     X_umap = umap.fit_transform(X_contiguous)
    # adata.obsm['X_umap'] = X_umap  # annotate samples with UMAP coordinates
    # logg.info(
    #     '    finished',
    #     time=start,
    #     deep=(
    #         'added\n'
    #         "    'X_umap', UMAP coordinates (adata.obsm)"
    #     ),
    # )
#########################################################################