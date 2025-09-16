#!/usr/bin/env codon
import sys
import os
from typing import Dict, List, Set, Optional


def reverse_complement(key: str) -> str:
    """
    Return the reverse complement of a DNA sequence.
    """
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

    key_list = list(key[::-1])
    for i in range(len(key_list)):
        key_list[i] = complement[key_list[i]]
    return ''.join(key_list)


class Node:
    def __init__(self, kmer: str):
        self._children: Set[int] = set()
        self._count: int = 0
        self.kmer: str = kmer
        self.visited: bool = False
        self.depth: int = 0
        self.max_depth_child: Optional[int] = None

    def add_child(self, child_idx: int) -> None:
        self._children.add(child_idx)

    def increase(self) -> None:
        self._count += 1

    def reset(self) -> None:
        self.visited = False
        self.depth = 0
        self.max_depth_child = None

    def get_count(self) -> int:
        return self._count

    def get_children(self) -> List[int]:
        return list(self._children)

    def remove_children(self, target: Set[int]) -> None:
        self._children = self._children - target


class DBG:
    def __init__(self, k: int, data_list: List[List[str]]):
        self.k: int = k
        self.nodes: Dict[int, Node] = {}
        self.kmer2idx: Dict[str, int] = {}
        self.kmer_count: int = 0

        self._check(data_list)
        self._build(data_list)

    def _check(self, data_list: List[List[str]]) -> None:
        assert len(data_list) > 0
        assert self.k <= len(data_list[0][0])

    def _build(self, data_list: List[List[str]]) -> None:
        for data in data_list:
            for original in data:
                rc = reverse_complement(original)
                for i in range(len(original) - self.k - 1):
                    self._add_arc(original[i:i + self.k], original[i + 1:i + 1 + self.k])
                    self._add_arc(rc[i:i + self.k], rc[i + 1:i + 1 + self.k])

    def show_count_distribution(self) -> None:
        count = [0] * 30
        for idx in self.nodes:
            node_count = self.nodes[idx].get_count()
            if node_count < 30:
                count[node_count] += 1
        print(count[0:10])

    def _add_node(self, kmer: str) -> int:
        if kmer not in self.kmer2idx:
            self.kmer2idx[kmer] = self.kmer_count
            self.nodes[self.kmer_count] = Node(kmer)
            self.kmer_count += 1

        idx = self.kmer2idx[kmer]
        self.nodes[idx].increase()
        return idx

    def _add_arc(self, kmer1: str, kmer2: str) -> None:
        idx1 = self._add_node(kmer1)
        idx2 = self._add_node(kmer2)
        self.nodes[idx1].add_child(idx2)

    def _get_count(self, child_idx: int) -> int:
        return self.nodes[child_idx].get_count()

    def _get_sorted_children(self, idx: int) -> List[int]:
        children = self.nodes[idx].get_children()
        children.sort(key=self._get_count, reverse=True)
        return children

    def _get_depth(self, idx: int) -> int:
        if not self.nodes[idx].visited:
            self.nodes[idx].visited = True
            children = self._get_sorted_children(idx)
            max_depth = 0
            max_child = None

            for child in children:
                depth = self._get_depth(child)
                if depth > max_depth:
                    max_depth = depth
                    max_child = child

            self.nodes[idx].depth = max_depth + 1
            self.nodes[idx].max_depth_child = max_child

        return self.nodes[idx].depth

    def _reset(self) -> None:
        for idx in self.nodes:
            self.nodes[idx].reset()

    def _get_longest_path(self) -> List[int]:
        max_depth = 0
        max_idx = None

        for idx in self.nodes:
            depth = self._get_depth(idx)
            if depth > max_depth:
                max_depth = depth
                max_idx = idx

        path = []
        while max_idx is not None:
            path.append(max_idx)
            max_idx = self.nodes[max_idx].max_depth_child

        return path

    def _delete_path(self, path: List[int]) -> None:
        for idx in path:
            del self.nodes[idx]

        path_set = set(path)
        for idx in self.nodes:
            self.nodes[idx].remove_children(path_set)

    def _concat_path(self, path: List[int]) -> Optional[str]:
        if len(path) < 1:
            return None

        concat = self.nodes[path[0]].kmer
        for i in range(1, len(path)):
            concat += self.nodes[path[i]].kmer[-1]

        return concat

    def get_longest_contig(self) -> Optional[str]:
        self._reset()
        path = self._get_longest_path()
        contig = self._concat_path(path)
        self._delete_path(path)
        return contig