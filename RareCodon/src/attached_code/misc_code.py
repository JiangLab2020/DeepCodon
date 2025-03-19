from collections import Counter

import Levenshtein
import numpy as np
import pandas as pd


def codon_to_vector(codon_string, rare_table):
    vectors = ""

    for i in range(0, len(codon_string), 3):
        codon = codon_string[i : i + 3]
        vectors += str(rare_table.get(codon))

    return vectors


def vector_with_slidingWindow(vector, sliding_window=5):
    length = len(vector)
    result = []

    for i in range(length):
        start = max(0, i - sliding_window // 2)
        end = min(length, i + sliding_window // 2 + 1)
        window = vector[start:end]
        majority = 0
        for w in window:
            if w == "1":
                majority += 1
        majority = str(majority)
        result.append(majority)
    return "".join(result)


def levenshtein_similarity(str1, str2):
    max_len = max(len(str1), len(str2))
    if max_len == 0:
        return 1.0
    distance = Levenshtein.distance(str1, str2)
    return 1 - (distance / max_len)


def filter_by_levenshtein_similarity(df, col_name, keep_ratio=0.8):
    data = df[col_name].values

    n = len(data)
    similarities = np.zeros((n, n))

    for i in range(n):
        for j in range(i, n):
            sim = levenshtein_similarity(data[i], data[j])
            similarities[i, j] = sim
            similarities[j, i] = sim

    avg_similarities = similarities.mean(axis=1)
    num_to_keep = int(np.ceil(keep_ratio * n))
    sorted_indices = np.argsort(avg_similarities)[::-1]
    keep_indices = sorted_indices[:num_to_keep]
    return df.iloc[keep_indices].copy()


def rare_codon_AND(df, col_name):
    strings = df[col_name].tolist()

    if len(strings) == 0:
        return ""
    result = strings[0]
    for s in strings[1:]:
        result = "".join(
            [
                "1" if result[i] == "1" and s[i] == "1" else "0"
                for i in range(len(result))
            ]
        )

    return result


def rare_codon_addAND(df, col_name):
    strings = df[col_name].tolist()
    if len(strings) == 0:
        return ""
    result = np.zeros(len(strings[0]), dtype=int)
    for s in strings:
        result += np.array([int(c) for c in s])
    result = np.floor(result / len(strings)).astype(int)
    return "".join(map(str, result))
