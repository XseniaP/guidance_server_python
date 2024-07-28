from typing import List, Dict
from collections import OrderedDict
import logging
import pandas as pd
import warnings
import networkx as nx
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
from itertools import groupby
import glob

logging.basicConfig(level=logging.WARNING)

ORDERED_MSA_FIELD_NAME = "_Ordered_MSA"
LIST_OF_POINTS_FIELD_NAME = "_List_of_Points"
CERTAINTY_SER_FIELD_NAME = "_Certainty_Ser"
MSA_CERTAINTY_FIELD_NAME = "_Overall_Certainty"
BEST_TAILORED_FIELD_NAME = "_Best_Tailroed_MSA"
GRAPH_WEIGHTS_FIELD_NAME = "_Graph_Weights"
PATH_ON_GRAPH_FIELD_NAME = "_Path_Graph"

DEFAULT_ALIGNMENT_NAME_FIELD = "name"
DEFAULT_ALIGNMENT_DATA_FIELD = "MSA"

NAME_OF_SOURCE_NODE = "_source"
NAME_OF_TARGET_NODE = "_target"
GRAPH_WEIGHT_FIELD_NAME = "weight"


def fasta_reader(fasta_path: str):
    """
    based on: https://www.biostars.org/p/710/
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_path)
    file_name = os.path.basename(fasta_path)
    msa_as_dict = {
        DEFAULT_ALIGNMENT_NAME_FIELD: file_name,
        DEFAULT_ALIGNMENT_DATA_FIELD: {}
    }
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = next(header)[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in next(faiter))
        if header in msa_as_dict[DEFAULT_ALIGNMENT_DATA_FIELD].keys():
            raise ValueError(f"file {fasta_path} has multiple rows with the same value. Specifically, "
                             f"header = {header} appeared twice -> while trying to add it to previous rows: "
                             f"{list(msa_as_dict.keys())}")
        msa_as_dict[DEFAULT_ALIGNMENT_DATA_FIELD][header] = seq
    return msa_as_dict


def get_MSAs_list_from_paths(folder_path: str):
    if os.path.isdir(folder_path):
        MSAs = []
        for file_path in glob.glob(f"{folder_path}/*"):
            try:
                MSAs.append(fasta_reader(file_path))
            except Exception as e:
                logging.error(f'cannot read fasta content from: {file_path}')
                logging.error(e)
        if len(MSAs) == 0:
            raise ValueError(f'no valid fasta files in the folder {folder_path}')
        return MSAs
    raise ValueError(f'the given path is not valid: {folder_path}')


def calc_char_position(row: str, index: int, is_number_gaps: bool):
    if is_number_gaps:
        gaps_amount = row.count('-')
        num_of_chars = index - gaps_amount
        if row[index] == '-':
            return f'-{num_of_chars}'
        return num_of_chars
    else:
        if row[index] == '-':
            return '-'
        gaps_amount = row.count('-')
        return index - gaps_amount


def convert_msa_2_list_of_points(msa: OrderedDict, is_number_gaps: bool = False) -> List:
    columns = []
    MSA_len = len(list(msa.items())[0][1])
    MSA_num_of_rows = len(msa.items())
    MSA_list_of_items = list(msa.items())
    for i in range(MSA_len):
        col = []
        for j in range(MSA_num_of_rows):
            col.append(calc_char_position(MSA_list_of_items[j][1][:i + 1], i, is_number_gaps))
        current_tuple = tuple(col)
        columns.append(current_tuple)

    # print MSA
    for (key, val) in MSA_list_of_items:
        logging.debug(val)

    # print points
    for i in range(MSA_num_of_rows):
        str2print = ''
        for col in columns:
            str2print += f'{str(col[i])} '
        logging.debug(str2print[:-1])
    return columns


def calc_alignment_certainty(MSAs: List[Dict], name_field: str = "name", data_field: str = "MSA") -> pd.Series:
    if not verify_alignment_list(MSAs=MSAs, name_field=name_field, data_field=data_field):
        raise ValueError(f'failed to valid MSAs list')

    for msa in MSAs:
        msa[ORDERED_MSA_FIELD_NAME] = OrderedDict(sorted(msa[data_field].items()))
        msa[LIST_OF_POINTS_FIELD_NAME] = convert_msa_2_list_of_points(msa[ORDERED_MSA_FIELD_NAME], False)

    logging.debug(MSAs)
    overall_certainty = {}
    for msa in MSAs:
        other_MSAs = MSAs[:]  # fastest way to copy
        other_MSAs.remove(msa)
        list_of_columns = msa[LIST_OF_POINTS_FIELD_NAME]
        result_dict = {}
        for idx_col, col in enumerate(list_of_columns):
            col_score = 0
            for other_msa in other_MSAs:
                if col in other_msa[LIST_OF_POINTS_FIELD_NAME]:
                    col_score += 1 / len(other_MSAs)
            result_dict[f"row_{idx_col}"] = col_score
        ser = pd.Series(data=result_dict, index=list(result_dict.keys()))
        msa[CERTAINTY_SER_FIELD_NAME] = ser
        msa[MSA_CERTAINTY_FIELD_NAME] = ser.mean()
        overall_certainty[msa[name_field]] = ser.mean()
    overall_certainty = pd.Series(data=overall_certainty, index=list(overall_certainty.keys()))
    overall_certainty = overall_certainty.sort_values(ascending=False)
    return overall_certainty


def get_alignment_with_max_certainty(
        MSAs: List[Dict] = None,
        path_to_folder_with_fasta_files: str = None,
        name_field: str = "name",
        data_field: str = "MSA",
        is_visualize_best_msa: bool = False,
        path_to_save_png: str = None,
    ) -> Dict:
    if not MSAs and not path_to_folder_with_fasta_files:
        raise ValueError(f'please provide a list of MSAs or a path to folder with fasta files')

    if not MSAs:
        MSAs = get_MSAs_list_from_paths(path_to_folder_with_fasta_files)

    if not verify_alignment_list(MSAs=MSAs, name_field=name_field, data_field=data_field):
        raise ValueError(f'failed to valid MSAs list')

    overall_certainty = calc_alignment_certainty(MSAs, name_field, data_field)
    top_alignment_name = overall_certainty.index[0]

    for msa in MSAs:
        if msa[name_field] == top_alignment_name:
            if is_visualize_best_msa:
                visualize_msa(target_msa=msa, MSAs=MSAs, path_to_save_png=path_to_save_png, data_field=data_field, name_field=name_field)

            return msa
    warnings.warn(
        f"top alignment name ({top_alignment_name}) isn't in list of MSAs names: {[m[name_field] for m in MSAs]}")


def add_edge_or_update_weight(graph, node1, node2, num_of_alignments):
    if graph.has_edge(node1, node2):
        graph[node1][node2][GRAPH_WEIGHT_FIELD_NAME] -= 1 / num_of_alignments
    else:
        graph.add_edge(node1, node2, weight=1)


def add_nodes_and_edges(graph, list_of_points, num_of_alignments):
    last_node = NAME_OF_SOURCE_NODE
    for col in list_of_points:
        node_name = '_'.join([str(i) for i in col])  # first column of the alignment
        add_edge_or_update_weight(graph, last_node, node_name, num_of_alignments)
        last_node = node_name
    add_edge_or_update_weight(graph, last_node, NAME_OF_TARGET_NODE, num_of_alignments)


def build_graph_from_alignments(
        MSAs: List[Dict],
        name_field: str = "name",
        data_field: str = "MSA"
    ) -> nx.DiGraph:
    graph = nx.DiGraph()
    graph.add_node(NAME_OF_SOURCE_NODE)
    graph.add_node(NAME_OF_TARGET_NODE)

    for msa in MSAs:
        msa[ORDERED_MSA_FIELD_NAME] = OrderedDict(sorted(msa[data_field].items()))
        msa[LIST_OF_POINTS_FIELD_NAME] = convert_msa_2_list_of_points(msa[ORDERED_MSA_FIELD_NAME], True)

    logging.info(MSAs)

    for i, msa in enumerate(MSAs):
        add_nodes_and_edges(graph, msa[LIST_OF_POINTS_FIELD_NAME], len(MSAs))

    return graph


def get_tailored_alignment(
        MSAs: List[Dict] = None,
        path_to_folder_with_fasta_files: str = None,
        name_field: str = "name",
        data_field: str = "MSA",
        is_calc_certainty_for_all_paths: bool = False,
        is_visualize_best_msa: bool = False,
        path_to_save_png: str = None,
    ) -> dict:
    if not MSAs and not path_to_folder_with_fasta_files:
        raise ValueError(f'please provide a list of MSAs or a path to folder with fasta files')

    if not MSAs:
        MSAs = get_MSAs_list_from_paths(path_to_folder_with_fasta_files)

    if not verify_alignment_list(MSAs=MSAs, name_field=name_field, data_field=data_field):
        raise ValueError(f'failed to valid MSAs list')

    unaligned_sequences = OrderedDict(sorted(MSAs[0][data_field].items()))
    unaligned_sequences = OrderedDict({key: value.replace("-", "") for key, value in unaligned_sequences.items()})

    graph = build_graph_from_alignments(MSAs=MSAs, name_field=name_field, data_field=data_field)

    best_path = nx.shortest_path(graph, NAME_OF_SOURCE_NODE, NAME_OF_TARGET_NODE, weight=GRAPH_WEIGHT_FIELD_NAME)
    best_alignment = calc_alignment_from_graph_path(unaligned_sequences=unaligned_sequences, path_on_graph=best_path, name_field=name_field, data_field=data_field)
    best_alignment[GRAPH_WEIGHTS_FIELD_NAME] = calc_edges_weight_on_path(graph, best_path, len(MSAs))

    if is_visualize_best_msa:
        visualize_msa(target_msa=best_alignment, MSAs=MSAs, path_to_save_png=path_to_save_png, data_field=data_field, name_field=name_field)

    # if not return certainty for all alignments, return the best
    if not is_calc_certainty_for_all_paths:
        return best_alignment

    for msa in MSAs:
        msa[PATH_ON_GRAPH_FIELD_NAME] = ['_'.join([str(c) for c in point]) for point in msa[LIST_OF_POINTS_FIELD_NAME]]
        msa[GRAPH_WEIGHTS_FIELD_NAME] = calc_edges_weight_on_path(graph=graph, path=msa[PATH_ON_GRAPH_FIELD_NAME], num_of_alignments=len(MSAs))

    return MSAs


def calc_edges_weight_on_path(graph: nx.Graph, path: List, num_of_alignments) -> float:
    last_node = ''
    path_certainty = []
    path_weights = []
    for i, node in enumerate(path):
        if i != 0:
            weight = graph[last_node][node][GRAPH_WEIGHT_FIELD_NAME]
            column_certainty = 1 - weight + 1 /num_of_alignments
            path_certainty.append(column_certainty)
            path_weights.append(weight)
        last_node = node
    return sum(path_weights) / len(path_weights)


def calc_alignment_from_graph_path(unaligned_sequences: OrderedDict, path_on_graph: List, name_field: str, data_field: str) -> Dict:
    unaligned_sequences_list = list(unaligned_sequences.values())
    aligned_result = [''] * len(unaligned_sequences_list)
    if not NAME_OF_SOURCE_NODE == path_on_graph[0]:
        warnings.warn(
            f'first node in path should be {NAME_OF_SOURCE_NODE}, received value of {path_on_graph[0]} -> invalid path')
    if not NAME_OF_TARGET_NODE == path_on_graph[-1]:
        warnings.warn(
            f'last node in path should be {NAME_OF_TARGET_NODE}, received value of {path_on_graph[-1]} -> invalid path')

    # don't check the first and last nodes as they are static
    for node in path_on_graph[1:-1]:
        columns_values = node.split('_')
        for col_idx, col_value in enumerate(columns_values):
            if col_value[0] == '-':
                aligned_result[col_idx] += '-'
            else:
                aligned_result[col_idx] += unaligned_sequences_list[col_idx][int(col_value)]

    aligned_dict = {}
    aligned_dict[data_field] = {key: alignment_str for key, alignment_str in zip(list(unaligned_sequences.keys()), aligned_result)}
    aligned_dict[name_field] = BEST_TAILORED_FIELD_NAME

    return aligned_dict


def verify_alignment_list(MSAs: List[Dict], name_field: str = "name", data_field: str = "MSA") -> bool:
    general_unaligned_sequences = OrderedDict(sorted(MSAs[0][data_field].items()))
    general_unaligned_sequences = OrderedDict(
        {key: value.replace("-", "") for key, value in general_unaligned_sequences.items()})

    for msa in MSAs:
        specific_unaligned_sequences = OrderedDict(sorted(msa[data_field].items()))
        specific_unaligned_sequences = OrderedDict(
            {key: value.replace("-", "") for key, value in specific_unaligned_sequences.items()})
        for general_item, specific_item in zip(general_unaligned_sequences.items(),
                                               specific_unaligned_sequences.items()):
            if general_item != specific_item:
                raise ValueError(f"some of your MSAs doesn't have the same rows names or don't have the same "
                                 f"unaligned sequenes. Specifically, MSA = {MSAs[0][name_field]} has a row named "
                                 f"{general_item[0]} and unaligned sequence of {general_item[1]} while the other MSA "
                                 f"{msa[name_field]} have a row named {specific_item[0]} with unaligned sequence of "
                                 f"{specific_item[1]}")

        len_first_row = len(list(msa[data_field].values())[0])
        for row in list(msa[data_field].values()):
            if len(row) != len_first_row:
                raise ValueError(f"rows for MSA {msa[name_field]} doesn't have the same length."
                                 f"There is a row with a length of {len_first_row} and a row with a length of {len(row)}")

    return True


def visualize_msa(
        target_msa: Dict,
        MSAs: List[Dict],
        path_to_save_png: str,
        data_field: str = "MSA",
        name_field: str = "name"
    ):
    if CERTAINTY_SER_FIELD_NAME not in target_msa:
        other_MSAs = MSAs[:]  # fastest way to copy
        if target_msa in MSAs:
            other_MSAs.remove(target_msa)

        for msa in MSAs:
            if ORDERED_MSA_FIELD_NAME not in msa:
                msa[ORDERED_MSA_FIELD_NAME] = OrderedDict(sorted(msa[data_field].items()))
            if LIST_OF_POINTS_FIELD_NAME not in msa:
                msa[LIST_OF_POINTS_FIELD_NAME] = convert_msa_2_list_of_points(msa[ORDERED_MSA_FIELD_NAME], True)

        if ORDERED_MSA_FIELD_NAME not in target_msa:
            target_msa[ORDERED_MSA_FIELD_NAME] = OrderedDict(sorted(target_msa[data_field].items()))
        if LIST_OF_POINTS_FIELD_NAME not in target_msa:
            target_msa[LIST_OF_POINTS_FIELD_NAME] = convert_msa_2_list_of_points(target_msa[ORDERED_MSA_FIELD_NAME], True)

        list_of_columns = target_msa[LIST_OF_POINTS_FIELD_NAME]
        result_dict = {}
        for idx_col, col in enumerate(list_of_columns):
            col_score = 0
            for other_msa in other_MSAs:
                if col in other_msa[LIST_OF_POINTS_FIELD_NAME]:
                    col_score += 1 / len(other_MSAs)
            result_dict[f"row_{idx_col}"] = col_score
        ser = pd.Series(data=result_dict, index=list(result_dict.keys()))
        target_msa[CERTAINTY_SER_FIELD_NAME] = ser
        target_msa[MSA_CERTAINTY_FIELD_NAME] = ser.mean()

    msa_np = np.array([])
    y_labels = []
    for row_idx, (title, row) in enumerate(target_msa[ORDERED_MSA_FIELD_NAME].items()):
        if row_idx == 0:
            msa_np = np.array([list(row)])
        else:
            msa_np = np.append(msa_np, [list(row)], axis=0)
        y_labels.append(title)
    labels_np = target_msa[CERTAINTY_SER_FIELD_NAME].to_numpy()
    x_labels = labels_np.copy()
    x_labels = np.around(x_labels, decimals=2)
    labels_np = labels_np.reshape((1, -1))
    labels_np = np.repeat(labels_np, [len(target_msa[ORDERED_MSA_FIELD_NAME].items())], axis=0)

    s = sns.heatmap(1 - labels_np, annot=msa_np, cmap='RdYlGn_r',
                vmax=1, vmin=0, xticklabels=x_labels, yticklabels=y_labels, cbar=False, fmt='', linewidths=2, square=True)
    s.set(xlabel='Certainty')
    s.set_title(target_msa[name_field])
    fig = s.get_figure()
    number_of_rows = labels_np.shape[0]
    number_of_columns = labels_np.shape[1]
    fig.set_size_inches(number_of_columns * 0.3 + 2, number_of_rows * 0.3 + 2)
    plt.tight_layout()
    if os.path.isdir(path_to_save_png):
        logging.info(f'saving visualization of the alignment to the following folder: {path_to_save_png}')
        fig.savefig(os.path.join(path_to_save_png, target_msa[name_field] + '.png'), bbox_inches='tight')
    elif path_to_save_png is not None:
        logging.info(f'saving visualization of the alignment to: {path_to_save_png}')
        fig.savefig(os.path.join(path_to_save_png, target_msa[name_field] + '.png'), bbox_inches='tight')
    else:
        logging.info(f'path to save fig is none, showing the alignment')
        plt.show()
