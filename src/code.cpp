#include <Rcpp.h>
using namespace Rcpp;

// Find the Jaccard Similarity Coefficient between all cells (rows) of an input
// matrix of k-nearest neighbor IDs for each cell
//
// @param knn_ids, a num_cells by num_neighbors matrix of nearest-neighbor indices
// @return jaccards, a matrix in which each row is a tuple (from, to, jac) of indices
//                   and jaccard coefficients between each pair of cells in the dataset.
//
// Citation: This code is heavily inspired by that of Chen Hao in the GitHub
//           package "Rphenograph" at github.com/JinmiaoChenLab/Rphenograph/ as well as
//           the Python package for PhenoGraph, written by Jacob Levine and hosted at
//           github.com/dpeerlab/PhenoGraph/


// [[Rcpp::export]]
NumericMatrix find_jaccard_coefficients(NumericMatrix knn_ids) {
  int num_cells = knn_ids.nrow();
  int num_neighbors = knn_ids.ncol();
  NumericMatrix jaccards(num_cells * num_neighbors, 3);
  int row_index = 0;
  for (int cell = 0; cell < num_cells; cell++) {
    for (int neighbor = 0; neighbor < num_neighbors; neighbor++) {
      int neighbor_index = knn_ids(cell, neighbor) - 1;
      NumericVector node_cell = knn_ids(cell, _);
      NumericVector node_neighbor = knn_ids(neighbor_index , _);
      // count number of mutual neighbors between node_cell and node_neighbor
      int node_intersection = intersect(node_cell, node_neighbor).size();
      // count unique cells in the neighborhood of node_cell and node_neighbor
      int node_union = union_(node_cell, node_neighbor).size();
      // find jaccard coefficient and report tuples
      jaccards(row_index, 0) = cell + 1;
      jaccards(row_index, 1) = neighbor_index + 1;
      jaccards(row_index, 2) = (1.0 * node_intersection) / node_union;
      row_index++;
    }
  }

  return jaccards;
}
