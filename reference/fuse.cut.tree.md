# Cut Hierarchical Clustering Tree into Clusters

Divides the clustering tree into a specified number of clusters.

## Usage

``` r
fuse.cut.tree(tree, k)
```

## Arguments

- tree:

  Clustering tree

- k:

  Number of clusters

## Value

A vector indicating which cluster each element in the original data
frame belonged to

## Examples

``` r
# Example: Cutting small tree in 2 segments
tree <- matrix(c(
-1, -2,  49.53106,  49.53106,  1.14473,
-3, -4,  78.49604,  78.49604,  1.14473,
-5, -6, 147.07154, 147.07154,  1.14473,
1,  2,  72.98287, 201.00997,  1.14473,
4,  3, 106.38879, 454.47029,  1.14473
), ncol = 5, byrow = TRUE)

segments <- fuse.cut.tree(tree, 2)
segments
#> [1] 1 1 1 1 2 2
```
