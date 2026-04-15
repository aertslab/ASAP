#!/usr/bin/env python3

"""Variant utils script."""


import math
import keras
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy import stats
from pysam import FastaFile
from functools import lru_cache
from keras.utils import Sequence

import data
# import shap
# import shap_utils
from crested.utils import one_hot_encode_sequence


def get_variant_schema(schema):
    var_schema = {
        'snp': ['chr', 'pos', 'ref', 'alt', 'variant_id'],
        'peak': ['chr', 'start', 'end', 'pos', 'ref', 'alt', 'variant_id'],
        'indel': ['chr', 'start', 'end', 'pos', 'ref', 'alt', 'variant_id']
    }

    return var_schema[schema]


def load_variant_table(variants_table_path, schema):
    variants_table = pd.read_csv(variants_table_path, header=None, sep='\t', names=get_variant_schema(schema))
    variants_table['chr'] = variants_table['chr'].astype(str)

    has_chr_prefix = any('chr' in x.lower() for x in variants_table['chr'].tolist())
    if not has_chr_prefix:
        variants_table['chr'] = 'chr' + variants_table['chr']

    return variants_table


def create_shuffle_table(variants_table, random_seed=None, num_shuffles=None):
    if num_shuffles:
        total_shuffle = len(variants_table) * num_shuffles
        shuf_variants_table = variants_table.sample(total_shuffle,
                                                    random_state=random_seed,
                                                    ignore_index=True,
                                                    replace=True)
        shuf_variants_table['random_seed'] = np.random.permutation(len(shuf_variants_table))
    else:
        shuf_variants_table = pd.DataFrame()  # empty dataframe

    return shuf_variants_table


def get_variants(chrom, pos, input_length, chrom_sizes_dict):
    valid_chrom = chrom in chrom_sizes_dict

    if valid_chrom:
        flank = input_length // 2
        lower_check = (pos - flank > 0)
        upper_check = (pos + flank <= chrom_sizes_dict[chrom])
        in_bounds = lower_check and upper_check
        valid_variant = valid_chrom and in_bounds
        return valid_variant
    else:
        return False


def get_peaks(chrom, start, end, input_length, chrom_sizes_dict):
    valid_chrom = chrom in chrom_sizes_dict

    if valid_chrom:
        flank = (input_length - (end - start)) // 2
        lower_check = (start - flank > 0)
        upper_check = (end + flank <= chrom_sizes_dict[chrom])
        in_bounds = lower_check and upper_check
        valid_peak = valid_chrom and in_bounds
        return valid_peak
    else:
        return False


class PeakGenerator(Sequence):
    def __init__(self,
                 peaks,
                 input_length,
                 genome_fasta,
                 batch_size=512
                 ):

        self.peaks = peaks
        self.num_peaks = self.peaks.shape[0]
        self.input_length = input_length
        self.batch_size = batch_size
        self.genome = FastaFile(genome_fasta)

    @lru_cache(maxsize=50000)
    def cached_fetch(self, chrom, start, end):
        return self.genome.fetch(chrom, start, end).upper()

    def __len__(self):
        return math.ceil(self.num_peaks / self.batch_size)

    def __getitem__(self, index):
        start_index = index * self.batch_size
        end_index = min(self.num_peaks, (index + 1) * self.batch_size)
        peaks_batch = self.peaks.iloc[start_index:end_index].reset_index(drop=True)

        # Generate sequence data
        seqs = self.__seq_generation__(peaks_batch)

        return seqs

    def __seq_generation__(self, peaks_df):
        seqs = np.empty((peaks_df.shape[0], self.input_length, 4), dtype=np.uint8)

        peaks_df['start_ext'], peaks_df['end_ext'] = data.extend_sequence(
            peaks_df['start'].values,
            peaks_df['end'].values,
            self.input_length
        )

        for index, row in peaks_df.iterrows():
            seqs[index, ] = one_hot_encode_sequence(
                self.cached_fetch(row['chr'], row['start_ext'], row['end_ext'])
            )

        return seqs


class VariantGenerator(Sequence):
    def __init__(self,
                 variants,
                 input_length,
                 genome_fasta,
                 batch_size=512,
                 shuffle=False,
                 schema='peak'
                 ):

        self.variants = variants
        self.num_variants = self.variants.shape[0]
        self.input_length = input_length
        self.batch_size = batch_size
        self.shuffle = shuffle
        self.schema = schema
        self.genome = FastaFile(genome_fasta)

    @lru_cache(maxsize=100000)
    def cached_fetch(self, chrom, start, end):
        return self.genome.fetch(chrom, start, end).upper()

    def __len__(self):
        return math.ceil(self.num_variants / self.batch_size)

    def __getitem__(self, index):
        start_index = index * self.batch_size
        end_index = min(self.num_variants, (index + 1) * self.batch_size)
        variants_batch = self.variants.iloc[start_index:end_index].reset_index(drop=True)
        variant_ids = variants_batch['variant_id'].tolist()

        # Generate sequence data
        if self.shuffle:
            ref_seqs, alt_seqs, snp_pos = self.__allele_seq_generation__(variants_batch, variants_batch.random_seed)
        else:
            ref_seqs, alt_seqs, snp_pos = self.__allele_seq_generation__(variants_batch)

        return variant_ids, ref_seqs, alt_seqs, snp_pos

    def __allele_seq_generation__(self, variants_df):

        ref_seqs = np.empty((variants_df.shape[0], self.input_length, 4), dtype=np.uint8)
        alt_seqs = np.empty((variants_df.shape[0], self.input_length, 4), dtype=np.uint8)

        if self.schema == 'snp':
            snp_pos = []
            variants_df['start_ext'], variants_df['end_ext'] = data.extend_sequence(
                (variants_df['pos'] - 1).values,
                variants_df['pos'].values,
                self.input_length
            )
            variants_df['rel_pos'] = variants_df['pos'] - variants_df['start_ext'] - 1

            for index, row in variants_df.iterrows():
                snp = row['rel_pos']
                snp_pos.append(snp)
                ref_seqs[index, ] = alt_seqs[index, ] = one_hot_encode_sequence(
                    self.cached_fetch(row['chr'], row['start_ext'], row['end_ext'])
                )

                # assert ((ref_seqs[index, self.input_length // 2, :] == one_hot_encode_sequence(row['ref'])).sum() == 4), \
                #    f"Got {ref_seqs[index, self.input_length // 2, :]} and {one_hot_encode_sequence(row['ref'])} at index {index}"

                ref_seqs[index, self.input_length // 2, :] = one_hot_encode_sequence(row['ref'])
                alt_seqs[index, self.input_length // 2, :] = one_hot_encode_sequence(row['alt'])

        elif self.schema == 'peak_old':
            snp_pos = []
            variants_df['start_ext'], variants_df['end_ext'] = data.extend_sequence(
                variants_df['start'].values,
                variants_df['end'].values,
                self.input_length
            )
            variants_df['rel_pos'] = variants_df['pos'] - variants_df['start_ext'] - 1

            for index, row in variants_df.iterrows():
                snp = row['rel_pos']
                snp_pos.append(snp)
                ref_seqs[index, ] = alt_seqs[index, ] = one_hot_encode_sequence(
                    self.cached_fetch(row['chr'], row['start_ext'], row['end_ext'])
                )

                assert ((ref_seqs[index, snp, :] == one_hot_encode_sequence(row['ref'])).sum() == 4), \
                    f"Got {ref_seqs[index, snp, :]} and {one_hot_encode_sequence(row['ref'])} at index {index}"

                alt_seqs[index, snp, :] = one_hot_encode_sequence(row['alt'])

        elif self.schema == 'peak':
            snp_pos = []

            variants_df['ref'] = variants_df['ref'].str.upper()
            variants_df['alt'] = variants_df['alt'].str.upper()
            variants_df['ref_len'], variants_df['alt_len'] = variants_df['ref'].str.len(), variants_df['alt'].str.len()
            variants_df['indel_len'] = variants_df[['ref_len', 'alt_len']].max(axis=1)
            max_indel_len = variants_df['indel_len'].max()

            variants_df['start_ext'], variants_df['end_ext'] = data.extend_sequence(
                variants_df['start'].values,
                variants_df['end'].values,
                self.input_length
            )
            variants_df['rel_pos'] = variants_df['pos'] - variants_df['start_ext'] - 1

            # handle indels outside of the peak
            variants_df = data.correct_edge_indels(variants_df)

            for index, row in variants_df.iterrows():
                # indels
                if row['indel_len'] > 1:
                    ref_seq, alt_seq, snp = data.get_indel_sequence(
                        row['chr'],
                        row['start_ext'],
                        row['end_ext'],
                        row['pos'],
                        row['rel_pos'],
                        row['ref'],
                        row['alt'],
                        row['ref_len'],
                        row['alt_len'],
                        genome=self.genome,
                        check_ref_sequence=True
                    )

                    snp_pos.append(snp)
                    ref_seqs[index, ] = one_hot_encode_sequence(ref_seq.upper())
                    alt_seqs[index, ] = one_hot_encode_sequence(alt_seq.upper())

                    if len(row['ref']) == 1:
                        assert ((ref_seqs[index, snp, :] == one_hot_encode_sequence(row['ref'])).sum() == 4), \
                            f"Got {ref_seqs[index, snp, :]} and {one_hot_encode_sequence(row['ref'])} at index {index}"

                # snps
                else:
                    snp = row['rel_pos']
                    snp_pos.append(snp)

                    ref_seqs[index, ] = alt_seqs[index, ] = one_hot_encode_sequence(
                        self.cached_fetch(
                            row['chr'], row['start_ext'], row['end_ext']
                        )
                    )

                    assert ((ref_seqs[index, snp, :] == one_hot_encode_sequence(row['ref'])).sum() == 4), \
                        f"Got {ref_seqs[index, snp, :]} and {one_hot_encode_sequence(row['ref'])} at index {index}"

                    alt_seqs[index, snp, :] = one_hot_encode_sequence(row['alt'])

        if self.shuffle:  # ToDo: incorporate combinatorial variants
            assert seed != -1

        return ref_seqs, alt_seqs, snp_pos


def fetch_peak_predictions(
    model,
    peaks,
    input_length,
    genome_fasta,
    batch_size,
    forward_only=False
):
    preds = np.zeros(shape=(peaks.shape[0], model.output_shape[1]))
    if not forward_only:
        revcomp_preds = np.zeros(shape=(peaks.shape[0], model.output_shape[1]))

    # peak sequence generator
    peak_gen = PeakGenerator(
        peaks=peaks,
        input_length=input_length,
        genome_fasta=genome_fasta,
        batch_size=batch_size
    )

    for i in tqdm(range(len(peak_gen))):
        seqs = peak_gen[i]
        revcomp_seqs = seqs[:, ::-1, ::-1]

        batch_preds = model.predict(seqs, verbose=False)
        if not forward_only:
            revcomp_batch_preds = model.predict(revcomp_seqs, verbose=False)

        preds[i * batch_size:(i + 1) * batch_size] = batch_preds
        if not forward_only:
            revcomp_preds[i * batch_size:(i + 1) * batch_size] = revcomp_batch_preds

    if not forward_only:
        average_preds = np.average([preds, revcomp_preds], axis=0)
        return average_preds
    else:
        return preds


def fetch_variant_predictions(
    model,
    variants_table,
    input_length,
    genome_fasta,
    batch_size,
    shuffle=False,
    forward_only=False,
    schema='peak'
):
    variant_ids = []
    ref_pred = np.zeros(shape=(variants_table.shape[0], model.output_shape[1]))
    alt_pred = np.zeros(shape=(variants_table.shape[0], model.output_shape[1]))
    if not forward_only:
        revcomp_ref_pred = np.zeros(shape=(variants_table.shape[0], model.output_shape[1]))
        revcomp_alt_pred = np.zeros(shape=(variants_table.shape[0], model.output_shape[1]))

    # variant sequence generator
    var_gen = VariantGenerator(
        variants=variants_table,
        input_length=input_length,
        genome_fasta=genome_fasta,
        batch_size=batch_size,
        shuffle=shuffle,
        schema=schema
    )

    for i in tqdm(range(len(var_gen))):

        batch_variant_ids, ref_seqs, alt_seqs, _ = var_gen[i]

        revcomp_ref_seqs = ref_seqs[:, ::-1, ::-1]
        revcomp_alt_seqs = alt_seqs[:, ::-1, ::-1]

        ref_batch_preds = model.predict(ref_seqs, verbose=False)
        alt_batch_preds = model.predict(alt_seqs, verbose=False)
        if not forward_only:
            revcomp_ref_batch_preds = model.predict(revcomp_ref_seqs, verbose=False)
            revcomp_alt_batch_preds = model.predict(revcomp_alt_seqs, verbose=False)

        ref_pred[i * batch_size:(i + 1) * batch_size] = ref_batch_preds
        alt_pred[i * batch_size:(i + 1) * batch_size] = alt_batch_preds
        if not forward_only:
            revcomp_ref_pred[i * batch_size:(i + 1) * batch_size] = revcomp_ref_batch_preds
            revcomp_alt_pred[i * batch_size:(i + 1) * batch_size] = revcomp_alt_batch_preds

        variant_ids.extend(batch_variant_ids)

    if not forward_only:
        average_ref_pred = np.average([ref_pred, revcomp_ref_pred], axis=0)
        average_alt_pred = np.average([alt_pred, revcomp_alt_pred], axis=0)
        return np.array(variant_ids), average_ref_pred, average_alt_pred
    else:
        return np.array(variant_ids), ref_pred, alt_pred


def fetch_peak_shap(
    model,
    peaks,
    input_length,
    genome_fasta,
    batch_size,
    normalize=True
):
    shaps = np.zeros(shape=(peaks.shape[0], model.output_shape[1]))

    # peak sequence generator
    peak_gen = PeakGenerator(
        peaks=peaks,
        input_length=input_length,
        genome_fasta=genome_fasta,
        batch_size=batch_size
    )

    # Get the contribution scores from the correct layer
    if model.layers[-1].activation is keras.activations.sigmoid:
        output = model.layers[-1].output
    elif model.layers[-1].activation is keras.activations.softmax:
        output = shap_utils.get_weighted_meannormed_logits(model)
    else:
        output = model.layers[-1].output

    for i in tqdm(range(len(peak_gen))):
        seqs = peak_gen[i]

        model_explainer = shap.explainers.deep.TFDeepExplainer(
            (model.input, output),
            shap_utils.shuffle_several_times,
            combine_mult_and_diffref=shap_utils.combine_mult_and_diffref
        )

        print("Generating shap scores")  # ToDo: correct shape and normalization
        shap_batch = model_explainer.shap_values(seqs, progress_message=100)
        shaps[i * batch_size:(i + 1) * batch_size] = shap_batch

    return shaps


def fetch_variant_shap(
    model,
    variants_table,
    input_length,
    genome_fasta,
    batch_size,
    shuffle=False,
    schema='peak',
    attribution_length=0,
    normalize=True
):
    variant_ids = []
    ref_shap = np.zeros(shape=(variants_table.shape[0], model.output_shape[1]))
    alt_shap = np.zeros(shape=(variants_table.shape[0], model.output_shape[1]))

    ref_peak_shap = np.zeros(shape=(variants_table.shape[0], model.output_shape[1]))
    alt_peak_shap = np.zeros(shape=(variants_table.shape[0], model.output_shape[1]))

    all_ref_shap = np.zeros(shape=(variants_table.shape[0], input_length, 4, model.output_shape[1]))
    all_ref_seqs = np.zeros(shape=(variants_table.shape[0], input_length, 4))
    all_alt_shap = np.zeros(shape=(variants_table.shape[0], input_length, 4, model.output_shape[1]))
    all_alt_seqs = np.zeros(shape=(variants_table.shape[0], input_length, 4))

    # variant sequence generator
    var_gen = VariantGenerator(
        variants=variants_table,
        input_length=input_length,
        genome_fasta=genome_fasta,
        batch_size=batch_size,
        shuffle=shuffle,
        schema=schema
    )

    # Get the contribution scores from the correct layer
    if model.layers[-1].activation is keras.activations.sigmoid:
        output = model.layers[-1].output
    elif model.layers[-1].activation is keras.activations.softmax:
        output = shap_utils.get_weighted_meannormed_logits(model)
    else:
        output = model.layers[-1].output

    for i in tqdm(range(len(var_gen))):

        batch_variant_ids, ref_seqs, alt_seqs, snp_pos = var_gen[i]

        model_explainer = shap.explainers.deep.TFDeepExplainer(
            (model.input, output),
            shap_utils.shuffle_several_times,
            combine_mult_and_diffref=shap_utils.combine_mult_and_diffref
        )

        print("Generating shap scores")
        ref_shap_batch = model_explainer.shap_values(ref_seqs, progress_message=100)
        alt_shap_batch = model_explainer.shap_values(alt_seqs, progress_message=100)
        ref_shap_batch = np.array(ref_shap_batch)
        alt_shap_batch = np.array(alt_shap_batch)

        if normalize:
            print('Normalizing gradients')
            ref_shap_batch = ref_shap_batch / np.sqrt(np.sum(np.sum(np.square(ref_shap_batch), axis=-1, keepdims=True), axis=-2, keepdims=True))
            alt_shap_batch = alt_shap_batch / np.sqrt(np.sum(np.sum(np.square(alt_shap_batch), axis=-1, keepdims=True), axis=-2, keepdims=True))

        all_ref_shap[i * batch_size:(i + 1) * batch_size] = np.moveaxis(ref_shap_batch, 0, -1)
        all_alt_shap[i * batch_size:(i + 1) * batch_size] = np.moveaxis(alt_shap_batch, 0, -1)
        all_ref_seqs[i * batch_size:(i + 1) * batch_size] = ref_seqs
        all_alt_seqs[i * batch_size:(i + 1) * batch_size] = alt_seqs

        # Get actual (not hypothetical) contribution scores
        ref_shap_batch = np.moveaxis(ref_shap_batch * ref_seqs, 0, -1)
        alt_shap_batch = np.moveaxis(alt_shap_batch * alt_seqs, 0, -1)
        indexes = np.arange(ref_shap_batch.shape[0])

        # ToTest
        peak_length = np.arange(input_length)
        all_pos, _ = np.meshgrid(peak_length, indexes)
        ref_peak_idx = np.argmax(ref_shap_batch[indexes.reshape(-1, 1), all_pos, :, :] != 0, axis=2)[:, :, 0]
        alt_peak_idx = np.argmax(alt_shap_batch[indexes.reshape(-1, 1), all_pos, :, :] != 0, axis=2)[:, :, 0]

        ref_peak_shap[i * batch_size:(i + 1) * batch_size] = ref_shap_batch[indexes.reshape(-1, 1), all_pos, ref_peak_idx, :].mean(axis=1)
        alt_peak_shap[i * batch_size:(i + 1) * batch_size] = alt_shap_batch[indexes.reshape(-1, 1), all_pos, alt_peak_idx, :].mean(axis=1)

        if attribution_length:
            print('Extending context of attribution scores')
            # Get the extended context
            indexes = indexes.reshape(-1, 1)
            range_slices = np.zeros(shape=(ref_shap_batch.shape[0], 2))
            range_slices[:, 0] = np.array(snp_pos) - attribution_length // 2
            range_slices[:, 1] = np.array(snp_pos) + attribution_length // 2
            range_slices = (range_slices[:, 0].reshape(-1, 1) + np.arange(attribution_length)).astype(int)

            # Get the actual (non-zero) attribution scores
            ref_idx = np.argmax(ref_shap_batch[indexes, range_slices, :, :] != 0, axis=2)[:, :, 0]
            alt_idx = np.argmax(alt_shap_batch[indexes, range_slices, :, :] != 0, axis=2)[:, :, 0]

            ref_shap[i * batch_size:(i + 1) * batch_size] = ref_shap_batch[indexes, range_slices, ref_idx, :].mean(axis=1)
            alt_shap[i * batch_size:(i + 1) * batch_size] = alt_shap_batch[indexes, range_slices, alt_idx, :].mean(axis=1)
        else:
            print('Get attribution scores for SNP only')
            # Get the actual (non-zero) attribution scores
            ref_idx = np.argmax(ref_shap_batch[indexes, np.array(snp_pos), :, :] != 0, axis=1)[:, 0]
            alt_idx = np.argmax(alt_shap_batch[indexes, np.array(snp_pos), :, :] != 0, axis=1)[:, 0]

            ref_shap[i * batch_size:(i + 1) * batch_size] = ref_shap_batch[indexes, np.array(snp_pos), ref_idx, :]
            alt_shap[i * batch_size:(i + 1) * batch_size] = alt_shap_batch[indexes, np.array(snp_pos), alt_idx, :]

        variant_ids.extend(batch_variant_ids)

    return np.array(variant_ids), (ref_shap, alt_shap), (ref_peak_shap, alt_peak_shap), ((all_ref_shap, all_ref_seqs), (all_alt_shap, all_alt_seqs))


def get_variant_scores(ref_pred, alt_pred, ref_peak_pred=None, alt_peak_pred=None):
    diff = alt_pred - ref_pred

    if ref_peak_pred is not None and alt_peak_pred is not None:  # Only happens for peak shap scores

        # (ToDo) Check if this makes sense. Would we get negative values here? Maybe subtract only the one term?
        logfc = np.log2(alt_pred - alt_pred.min(axis=0) - ref_pred.min(axis=0)) \
            - np.log2(ref_pred - ref_pred.min(axis=0) - alt_pred.min(axis=0))

        logfc_peak = np.log2(alt_peak_pred - alt_peak_pred.min(axis=0) - ref_peak_pred.min(axis=0)) \
            - np.log2(ref_peak_pred - ref_peak_pred.min(axis=0) - alt_peak_pred.min(axis=0))

        diff_peak = alt_peak_pred - ref_peak_pred
        return (logfc, diff), (logfc_peak, diff_peak)

    else:
        logfc = np.log2(alt_pred) - np.log2(ref_pred)
        return (logfc, diff), (None, None)


def get_variant_scores_with_peaks(ref_pred, alt_pred, preds, ref_peak_pred=None, alt_peak_pred=None):
    N, M = preds.shape
    ref_percentile = np.zeros(shape=ref_pred.shape)
    alt_percentile = np.zeros(shape=alt_pred.shape)

    for i in range(M):
        sorted_preds = np.sort(preds[:, i])
        ref_ranks = np.searchsorted(sorted_preds, ref_pred[:, i], side='right')
        alt_ranks = np.searchsorted(sorted_preds, alt_pred[:, i], side='right')

        ref_percentile[:, i] = (ref_ranks / N) * 100.0
        alt_percentile[:, i] = (alt_ranks / N) * 100.0

    if ref_peak_pred is not None and alt_peak_pred is not None:  # Only happens for peak shap scores
        (logfc, diff), (logfc_peak, diff_peak) = get_variant_scores(ref_pred, alt_pred, ref_peak_pred, alt_peak_pred)

        N, M = preds.shape
        ref_peak_percentile = np.zeros(shape=ref_peak_pred.shape)
        alt_peak_percentile = np.zeros(shape=alt_peak_pred.shape)

        for i in range(M):
            sorted_preds = np.sort(preds[:, i])
            ref_peak_ranks = np.searchsorted(sorted_preds, ref_peak_pred[:, i], side='right')
            alt_peak_ranks = np.searchsorted(sorted_preds, alt_peak_pred[:, i], side='right')

            ref_peak_percentile[:, i] = (ref_peak_ranks / N) * 100.0
            alt_peak_percentile[:, i] = (alt_peak_ranks / N) * 100.0

        return (logfc, diff, ref_percentile, alt_percentile), (logfc_peak, diff_peak, ref_peak_percentile, alt_peak_percentile)

    else:
        (logfc, diff), _ = get_variant_scores(ref_pred, alt_pred)
        return (logfc, diff, ref_percentile, alt_percentile), (None, None, None, None)
