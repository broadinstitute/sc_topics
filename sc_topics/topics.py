from typing import Optional, Union

import numpy as np
from anndata import AnnData
from numpy.random.mtrand import RandomState
from sklearn.utils import check_random_state, check_array

from ._utils import get_init_pos_from_paga, _choose_representation
from .. import logging as logg
from .._settings import settings
from .._compat import Literal

import gensim


_InitPos = Literal['paga', 'spectral', 'random']

def topics(
    adata: AnnData,
    model: str,
    num_topics: int,
    random_state: int,
    update_every: int, # only for ldamodel
    chunksize: int,
    passes: int,
    alpha: , # 'auto',
    per_word_topics: bool
) -> Optional[AnnData]:
    """\
    Topic modeling description and interface description will go here!


    Parameters
    ----------
    adata
        Annotated data matrix.
    model
        Topic model to use.
    num_topics
        Th number of topics, if needed to be pre-specified.
    random_state
        Random seed to use for reproducibility. 
        If `int`, `random_state` is the seed used by the random number generator;
        If `RandomState`, `random_state` is the random number generator;
        If `None`, the random number generator is the `RandomState` instance used
        by `np.random`.
    update_every
        Batch size
    alpha
        The initial learning rate for the embedding optimization.
    per_word_topics
        Output per word topics.
    copy
        Return a copy instead of writing to adata.
    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.
    **X_topics** : `adata.obsm` field
        cell topic scores of data.
    """

    adata = adata.copy() if copy else adata

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
    lda = models.LdaModel(corpus=corpus,
                      id2word=dictionary,
                      num_topics=t,
                      random_state=100,
                      update_every=0, # only for ldamodel
                      chunksize=10000,
                      passes=1,
                      alpha='auto',
                      per_word_topics=True
                     )


    return adata if copy else None

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