# SurvivalClusteringTree

Outcome-guided clustering by recursive partitioning methods.

## Installation

```
remotes::install_github("luyouepiusf/SurvivalClusteringTree")
```

## Introduction

Suppose that we are given a dataset of samples with survival outcomes and sample features. The main objective of the method is to identify several clusters of samples with similar survival outcomes and similar features. The package implements a novel statistical machine learning method to perform the task.

## Algorithm

1. Perform bootstrap aggregation. For each bootstrap dataset, build a survival tree by recursively sweeping through all variables to find the best split in the sample space, and repeat until no further splits can be found.

2. Define a distance between each pair of samples in the original dataset with respect to the survival tree. The distances between samples in the same terminal node are 0. For samples that are not in the same terminal node, there is a unique path on the tree connecting the two terminal nodes, the sum of the absolute log-rank test z-scores on the path is the distance between them.

3. Aggregate the distances by averaging pairwise distances as the final distance.

4. Apply hierarchical clustering analysis to the distance matrix to identify clusters.
