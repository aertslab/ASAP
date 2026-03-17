from __future__ import annotations

import logging
import os
import pickle
import random
import subprocess
import sys
import tempfile
import time
import warnings
from itertools import chain
from typing import TYPE_CHECKING

import lda
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
import ray
import tmtoolkit
from gensim import matutils, utils
from gensim.models import basemodel
from pycisTopic.utils import loglikelihood, subset_list
from scipy import sparse
import argparse
import polars as pl
import gzip
import numpy as np
from io import StringIO

if TYPE_CHECKING:
    from pycisTopic.cistopic_class import CistopicObject

## NOTE: this code assumes that you ran the mallet function 'train-topics', so that you have the following files in your tmp dir:
## - {number_of_topics}.topics__topickeys.txt
## - {number_of_topics}.topics_doctopics.txt
## - {number_of_topics}.topics_inferencer.mallet
## - {number_of_topics}.topics_state.mallet.gz
## to get these files, see 'train.sh'

print('importing done')
print('defining functions and classes....')

## Arguments
def make_argument_parser():
    """
    Creates an ArgumentParser to read the options for this script from
    sys.argv
    """
    parser = argparse.ArgumentParser(
        description="Load assigned topics and save topic model",)

    parser.add_argument('--tmp', '-tmp', type=str, required=True,
                        help='Tmp dir')
    
    parser.add_argument('--counts_path', '-counts', type=str, required=True,
                        help='Path to mallet counts fi;e.')
    
    parser.add_argument('--topics', '-t', type=int, required=True,
                        help='Number of topics.')

    parser.add_argument('--mallet', '-m', type=str, required=True,
                        help='Path to mallet.')

    parser.add_argument('--ncpu', '-ncpu', type=int, required=True,
                        help='Number of cpus')

    parser.add_argument('--save','-save', type=str, required=True, 
                        help ='path to save model')

    parser.add_argument('--doctopics', '-dt', type=str, required=True,
                        help='Path to doctopics file.')

    parser.add_argument('--cistopic_object_path', '-cto', type=str, required=True,
                        help='Path to cistopic_object.')
                        
    return parser


## CistopicLDAModel class: Nothing has changed here
class CistopicLDAModel:
    """
    cisTopic LDA model class

    :class:`cistopicLdaModel` contains model quality metrics (model coherence (adaptation from Mimno et al., 2011), log-likelihood (Griffiths and Steyvers, 2004), density-based (Cao Juan et al., 2009) and divergence-based (Arun et al., 2010)), topic quality metrics (coherence, marginal distribution and total number of assignments), cell-topic and topic-region distribution, model parameters and model dimensions.

    Parameters
    ----------
    metrics: pd.DataFrame
        :class:`pd.DataFrame` containing model quality metrics, including model coherence (adaptation from Mimno et al., 2011), log-likelihood and density and divergence-based methods (Cao Juan et al., 2009; Arun et al., 2010).
    coherence: pd.DataFrame
        :class:`pd.DataFrame` containing the coherence of each topic (Mimno et al., 2011).
    marginal_distribution: pd.DataFrame
        :class:`pd.DataFrame` containing the marginal distribution for each topic. It can be interpreted as the importance of each topic for the whole corpus.
    topic_ass: pd.DataFrame
        :class:`pd.DataFrame` containing the total number of assignments per topic.
    cell_topic: pd.DataFrame
        :class:`pd.DataFrame` containing the topic cell distributions, with cells as columns, topics as rows and the probability of each topic in each cell as values.
    topic_region: pd.DataFrame
        :class:`pd.DataFrame` containing the topic cell distributions, with topics as columns, regions as rows and the probability of each region in each topic as values.
    parameters: pd.DataFrame
        :class:`pd.DataFrame` containing parameters used for the model.
    n_cells: int
        Number of cells in the model.
    n_regions: int
        Number of regions in the model.
    n_topic: int
        Number of topics in the model.

    References
    ----------
    Mimno, D., Wallach, H., Talley, E., Leenders, M., & McCallum, A. (2011). Optimizing semantic coherence in topic models. In Proceedings of the 2011 Conference on Empirical Methods in Natural Language Processing (pp. 262-272).

    Griffiths, T. L., & Steyvers, M. (2004). Finding scientific topics. Proceedings of the National academy of Sciences, 101(suppl 1), 5228-5235.

    Cao, J., Xia, T., Li, J., Zhang, Y., & Tang, S. (2009). A density-based method for adaptive LDA model selection. Neurocomputing, 72(7-9), 1775-1781.

    Arun, R., Suresh, V., Madhavan, C. V., & Murthy, M. N. (2010). On finding the natural number of topics with latent dirichlet allocation: Some observations. In Pacific-Asia conference on knowledge discovery and data mining (pp. 391-402). Springer, Berlin, Heidelberg.

    """

    def __init__(
        self,
        metrics: pd.DataFrame,
        coherence: pd.DataFrame,
        marg_topic: pd.DataFrame,
        topic_ass: pd.DataFrame,
        cell_topic: pd.DataFrame,
        topic_region: pd.DataFrame,
        parameters: pd.DataFrame,
    ):
        self.metrics = metrics
        self.coherence = coherence
        self.marg_topic = marg_topic
        self.topic_ass = topic_ass
        self.cell_topic = cell_topic
        self.cell_topic_harmony = []
        self.topic_region = topic_region
        self.parameters = parameters
        self.n_cells = cell_topic.shape[1]
        self.n_regions = topic_region.shape[0]
        self.n_topic = cell_topic.shape[0]

    def __str__(self):
        descr = f"CistopicLDAModel with {self.n_topic} topics and n_cells × n_regions = {self.n_cells} × {self.n_regions}"
        return descr

class LDAMallet(utils.SaveLoad, basemodel.BaseTopicModel):
    """
    Wrapper class to run LDA models with Mallet. This class has been adapted from gensim (https://github.com/RaRe-Technologies/gensim/blob/27bbb7015dc6bbe02e00bb1853e7952ac13e7fe0/gensim/models/wrappers/ldamallet.py).

    Parameters
    ----------
    num_topics: int
        The number of topics to use in the model.
    corpus: iterable of iterable of (int, int), optional
        Collection of texts in BoW format. Default: None.
    alpha: float, optional
        Scalar value indicating the symmetric Dirichlet hyperparameter for topic proportions. Default: 50.
    id2word : :class:`gensim.utils.FakeDict`, optional
        Mapping between tokens ids and words from corpus, if not specified - will be inferred from `corpus`. Default: None.
    n_cpu : int, optional
        Number of threads that will be used for training. Default: 1.
    tmp_dir : str, optional
        tmp_dir for produced temporary files. Default: None.
    optimize_interval : int, optional
        Optimize hyperparameters every `optimize_interval` iterations (sometimes leads to Java exception 0 to switch off hyperparameter optimization). Default: 0.
    iterations : int, optional
        Number of training iterations. Default: 150.
    topic_threshold : float, optional
        Threshold of the probability above which we consider a topic. Default: 0.0.
    random_seed: int, optional
        Random seed to ensure consistent results, if 0 - use system clock. Default: 555.
    mallet_path: str
        Path to the mallet binary (e.g. /xxx/Mallet/bin/mallet). Default: "mallet".

    """

    def __init__(
        self,
        num_topics: int,
        corpus: list | None = None,
        alpha: float = 50,
        eta: float = 0.1,
        id2word: utils.FakeDict = None,
        n_cpu: int = 1,
        tmp_dir: str = None,
        optimize_interval: int = 0,
        iterations: int = 150,
        topic_threshold: float = 0.0,
        random_seed: int = 555,
        reuse_corpus: bool = False,
        mallet_path: str = "mallet",
        counts_file: str = None,
    ):
        logger = logging.getLogger("LDAMalletWrapper")
        if id2word is None:
            logger.warning(
                "No id2word mapping provided; initializing from corpus, assuming identity"
            )
            self.num_terms = utils.get_max_id(corpus) + 1
        else:
            self.num_terms = id2word.num_terms

        if self.num_terms == 0:
            raise ValueError("Cannot compute LDA over an empty collection (no terms)")

        self.num_topics = num_topics
        self.topic_threshold = topic_threshold
        self.alpha = alpha
        self.eta = eta
        self.tmp_dir = tmp_dir if tmp_dir else tempfile.gettempdir()
        self.random_label = hex(random.randint(0, 0xFFFFFF))[2:]
        self.n_cpu = n_cpu
        self.optimize_interval = optimize_interval
        self.iterations = iterations
        self.random_seed = random_seed
        self.mallet_path = mallet_path
        self.counts_file = counts_file

        self.word_topics = self.create_regions_topics_frequency_matrix(self.counts_file)
        
    def create_regions_topics_frequency_matrix(self, counts_file):
        """
        Create regions vs topics frequency matrix from Mallet region topics counts file.
        mallet_region_topics_counts_filename
            Mallet region topics counts file.
        """
        no_region_ids = -1
        no_topics = -1
        region_id_topic_counts = []

        with open(self.counts_file, "r") as fh:
            # Column 0: order in which region ID idx was seen in the input corpus file.
            # Column 1: region ID idx
            # Column 3-n: "topic:count" pairs
            #
            # Example:
            # --------
            #
            # 0 12 3:94 11:84 1:84 18:75 17:36 0:31 13:25 4:23 6:22 12:16 9:10 10:6 15:3 7:2 8:1
            # 1 28 8:368 15:267 3:267 17:255 0:245 10:227 16:216 19:201 7:92 18:85 1:58 14:52 9:31 6:17 13:6 2:3
            # 2 33 8:431 16:418 10:354 3:257 17:211 12:146 7:145 9:115 4:108 13:106 18:66 1:60 15:45 6:45 19:33 5:19 14:12 0:1
            # 3 35 7:284 18:230 15:199 10:191 16:164 0:114 4:112 19:107 12:104 13:68 3:49 9:35 1:28 11:25 5:20 17:17 6:11 14:2 8:1
            # 4 57 8:192 3:90 19:88 1:69 18:67 2:63 10:62 17:38 15:37 13:10 4:9 12:2 9:1
            for line in fh:
                columns = line.rstrip().split()
                # Get region ID index from second column.
                region_id_idx = int(columns[1])
                # Get topic index and counts from column 3 till the end by splitting: "topic:count" pairs.
                topics_counts = [
                    (int(topic), int(count))
                    for topic, count in [
                        topic_counts.split(":", 1) for topic_counts in columns[2:]
                    ]
                ]
                # Get topic indices.
                topics_idx = np.array([topic for topic, count in topics_counts])
                # Get counts.
                counts = np.array([count for topic, count in topics_counts])
                # Store region ID index, topics indices and counts till we know how many regions and topics we have.
                region_id_topic_counts.append((region_id_idx, topics_idx, counts))

                # Keep track of the highest seen region ID index and topic index (0-based).
                no_region_ids = max(region_id_idx, no_region_ids)
                no_topics = max(topics_idx.max(), no_topics)

        # Add 1 to region IDs and topics counts to account for start at 0.
        no_region_ids += 1
        no_topics += 1

        # Create regions topics counts matrix and populate it.
        regions_topics_counts = np.zeros((no_topics, no_region_ids), dtype=np.float64)
        for region_idx, topics_idx, counts in region_id_topic_counts:
            regions_topics_counts[topics_idx, region_idx] = counts

        # Create regions topics frequency matrix by dividing all count values for topic
        # by total counts for that topic.

        #np.save('/staging/leuven/stg_00090/ASA/analysis/2024_T2T_ATAC_analysis/SN/src/regions_topics_counts.npy', regions_topics_counts)
        return regions_topics_counts

    def get_topics(self):
        """
        Get topics X words matrix.

        Returns
        -------
        np.ndarray
            Topics X words matrix, shape `num_topics` x `vocabulary_size`.

        """

        regions_topics_counts = np.asarray(self.word_topics, np.float64)

        # Create regions topics frequency matrix by dividing all count values for topic
        # by total counts for that topic.
        regions_topics_frequency = (
            regions_topics_counts / regions_topics_counts.sum(axis=1)[:, None]
        ).astype(np.float32)
        #np.save('/staging/leuven/stg_00090/ASA/analysis/2024_T2T_ATAC_analysis/SN/src/regions_topics_frequency.npy', regions_topics_frequency)

        return regions_topics_frequency
        
## 'run_cgs_model_mallet' function: 
## added state_file and doc_topic_path arguments, so that it can read existing files from dir
def run_cgs_model_mallet(
    binary_matrix: sparse.csr_matrix,
    corpus: list,
    id2word: utils.FakeDict,
    cell_names: list[str],
    n_topics: list[int],
    region_names: list[str],
    n_cpu: int = 1,
    n_iter: int = 500,
    random_state: int = 555,
    alpha: float = 50,
    alpha_by_topic: bool = True,
    eta: float = 0.1,
    eta_by_topic: bool = False,
    top_topics_coh: int = 5,
    tmp_path: str = None,
    save_path: str = None,
    reuse_corpus: bool = False,
    counts_file: str = None,
    doc_topic_path: str = None,
    mallet_path: str = "mallet",
):
    """
    Run Latent Dirichlet Allocation in a model as implemented in Mallet (McCallum, 2002).

    Parameters
    ----------
    binary_matrix: sparse.csr_matrix
        Binary sparse matrix containing cells as columns, regions as rows, and 1 if a regions is considered accessible on a cell (otherwise, 0).
    n_topics: list of int
        A list containing the number of topics to use in each model.
    cell_names: list of str
        List containing cell names as ordered in the binary matrix columns.
    region_names: list of str
        List containing region names as ordered in the binary matrix rows.
    n_cpu: int, optional
        Number of cpus to use for modelling. In this function parallelization is done per model, that is, one model will run entirely in a unique cpu. We recommend to set the number of cpus as the number of models that will be inferred, so all models start at the same time.
    n_iter: int, optional
        Number of iterations for which the Gibbs sampler will be run. Default: 150.
    random_state: int, optional
        Random seed to initialize the models. Default: 555.
    alpha: float, optional
        Scalar value indicating the symmetric Dirichlet hyperparameter for topic proportions. Default: 50.
    alpha_by_topic: bool, optional
        Boolean indicating whether the scalar given in alpha has to be divided by the number of topics. Default: True
    eta: float, optional
        Scalar value indicating the symmetric Dirichlet hyperparameter for topic multinomials. Default: 0.1.
    eta_by_topic: bool, optional
        Boolean indicating whether the scalar given in beta has to be divided by the number of topics. Default: False
    top_topics_coh: int, optional
        Number of topics to use to calculate the model coherence. For each model, the coherence will be calculated as the average of the top coherence values. Default: 5.
    tmp_path: str, optional
        Path to a temporary folder for Mallet. Default: None.
    save_path: str, optional
        Path to save models as independent files as they are completed. This is recommended for large data sets. Default: None.
    reuse_corpus: bool, optional
        Whether to reuse the mallet corpus in the tmp directory. Default: False
    mallet_path: str
        Path to Mallet binary (e.g. "/xxx/Mallet/bin/mallet"). Default: "mallet".

    Return
    ------
    CistopicLDAModel
        A cisTopic LDA model.

    References
    ----------
    McCallum, A. K. (2002). Mallet: A machine learning for language toolkit. http://mallet.cs.umass.edu.

    """
    # Create cisTopic logger
    level = logging.INFO
    log_format = "%(asctime)s %(name)-12s %(levelname)-8s %(message)s"
    handlers = [logging.StreamHandler(stream=sys.stdout)]
    logging.basicConfig(level=level, format=log_format, handlers=handlers)
    log = logging.getLogger("cisTopic")

    # Set models
    if not alpha_by_topic:
        alpha = alpha * n_topics
    if eta_by_topic:
        eta = eta / n_topics

    # Running model
    ## removed start time
    #start = time.time()
    log.info(f"Running model with {n_topics} topics")
    
    ## added the extra state_file argument
    model = LDAMallet(
        corpus=corpus,
        id2word=id2word,
        num_topics=n_topics[0],
        iterations=n_iter,
        alpha=alpha,
        eta=eta,
        n_cpu=n_cpu,
        tmp_dir=tmp_path,
        random_seed=random_state,
        reuse_corpus=reuse_corpus,
        mallet_path=mallet_path,
        counts_file=counts_file,
    )
    
    # Get distributions
    topic_word = model.get_topics()
    doc_topic = (
        pd.read_csv(doc_topic_path, header=None, sep="\t").iloc[:, 2:].to_numpy()
    )

    print('model evaluation')
    # Model evaluation
    cell_cov = np.asarray(binary_matrix.sum(axis=0)).astype(float)
    arun_2010 = tmtoolkit.topicmod.evaluate.metric_arun_2010(
        topic_word, doc_topic, cell_cov
    )
    cao_juan_2009 = tmtoolkit.topicmod.evaluate.metric_cao_juan_2009(topic_word)
    mimno_2011 = tmtoolkit.topicmod.evaluate.metric_coherence_mimno_2011(
        topic_word,
        dtm=binary_matrix.transpose(),
        top_n=20,
        eps=1e-12,
        normalize=True,
        return_mean=False,
    )
    topic_word_assig = model.word_topics
    doc_topic_assig = (doc_topic.T * (cell_cov)).T
    ll = loglikelihood(topic_word_assig, doc_topic_assig, alpha, eta)

    # Organize data
    if len(mimno_2011) <= top_topics_coh:
        metrics = pd.DataFrame(
            [arun_2010, cao_juan_2009, np.mean(mimno_2011), ll],
            index=["Arun_2010", "Cao_Juan_2009", "Mimno_2011", "loglikelihood"],
            columns=["Metric"],
        ).transpose()
    else:
        metrics = pd.DataFrame(
            [
                arun_2010,
                cao_juan_2009,
                np.mean(
                    mimno_2011[
                        np.argpartition(mimno_2011, -top_topics_coh)[-top_topics_coh:]
                    ]
                ),
                ll,
            ],
            index=["Arun_2010", "Cao_Juan_2009", "Mimno_2011", "loglikelihood"],
            columns=["Metric"],
        ).transpose()
    
    ## n_topics is here a list, so we want the first value of the list (= number of topics for which the topic modelling is running)
    
    coherence = pd.DataFrame(
        [range(1, n_topics[0] + 1), mimno_2011], index=["Topic", "Mimno_2011"]
    ).transpose() ## Should iterate over all topics of the model, not just use the total number of topics
    
    marg_topic = pd.DataFrame(
        [
            range(1, n_topics[0] + 1),
            tmtoolkit.topicmod.model_stats.marginal_topic_distrib(doc_topic, cell_cov),
        ],
        index=["Topic", "Marg_Topic"],
    ).transpose() ## Should iterate over all topics of the model, not just use the total number of topics
    

    topic_ass = pd.DataFrame.from_records(
        [
            range(1, n_topics[0] + 1),
            list(chain.from_iterable(model.word_topics.sum(axis=1)[:, None])),
        ],
        index=["Topic", "Assignments"],
    ).transpose() ## Should iterate over all topics of the model, not just use the total number of topics

    cell_topic = pd.DataFrame.from_records(
        doc_topic,
        index=cell_names,
        columns=["Topic" + str(i) for i in range(1, n_topics[0] + 1)],
    ).transpose() ## Should iterate over all topics of the model, not just use the total number of topics

    topic_region = pd.DataFrame.from_records(
        topic_word,
        columns=region_names,
        index=["Topic" + str(i) for i in range(1, n_topics[0] + 1)],
    ).transpose() ## Should iterate over all topics of the model, not just use the total number of topics

    parameters = pd.DataFrame(
        [
            "Mallet",
            n_topics[0],
            n_iter,
            random_state,
            alpha,
            alpha_by_topic,
            eta,
            top_topics_coh,
           
        ],
        index=[
            "package",
            "n_topics",
            "n_iter",
            "random_state",
            "alpha",
            "alpha_by_topic",
            "eta",
            "top_topics_coh",
           
        ],
        columns=["Parameter"],
    )
    # Create object
    model = CistopicLDAModel(
        metrics, coherence, marg_topic, topic_ass, cell_topic, topic_region, parameters
    )
    log.info(f"Model with {n_topics} topics done!")
    ## removed making the directory and save directly at the save_path instead of concatenating dir and file name
    if isinstance(save_path, str):
        log.info(f"Saving model with {n_topics} topics at {save_path}")
        with gzip.open(save_path, "wb") as f:
            pickle.dump(model, f)
    return model
print('defining functions and classes done')
print()

print('executing main starts....')

def main():
    """
    The main executable function
    """
    print('make argument parser')
    parser = make_argument_parser()
    args = parser.parse_args()
    
    ## Print the arguments
    number_of_topics = args.topics
    print('Number of topics: ', str(number_of_topics))
    
    mallet_path = args.mallet
    print('mallet path: ', mallet_path)

    tmp_dir = args.tmp
    print('tmp dir: ', tmp_dir)

    counts_path = args.counts_path
    print('counts_path: ', counts_path)
    
    doctopics = args.doctopics
    print('doctopics: ', doctopics)
    
    save_path = args.save
    print('save_path: ', save_path)

    n_cpus = args.ncpu
    print('number of cpus: ', str(n_cpus))

    cistopic_obj_path = args.cistopic_object_path
    print('cistopic_obj: ', cistopic_obj_path)
    print('Opening cto....')
    infile = open(cistopic_obj_path, 'rb')
    cistopic_obj = pickle.load(infile)
    infile.close()
    print('Opening cto done')
    print('Initializing CTO variables.....')
    binary_matrix = cistopic_obj.binary_matrix
    region_names = cistopic_obj.region_names
    cell_names = cistopic_obj.cell_names
    corpus = matutils.Sparse2Corpus(binary_matrix)
    id2word = utils.FakeDict(len(region_names))
    print('Initializing variables done')
    print()
    print('start running mallet....')
    run_cgs_model_mallet(  
        binary_matrix = binary_matrix,
        corpus = corpus,
        id2word = id2word,
        n_topics = [number_of_topics],
        cell_names = cell_names,
        region_names = region_names,
        n_cpu = n_cpus,
        n_iter = 400,
        random_state = 555,
        alpha = 50,
        alpha_by_topic = True,
        eta = 0.1,
        eta_by_topic = False,
        top_topics_coh = 5,
        tmp_path = tmp_dir,
        save_path = save_path,
        reuse_corpus = True,
        counts_file = counts_path,
        doc_topic_path =  doctopics,
        mallet_path = mallet_path,)
    print('mallet finished')
    print('-------------------------')
    print('JOB DONE!')
    print('-------------------------')
  

if __name__ == "__main__":
    main()