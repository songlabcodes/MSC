#' MSC results from 8k PBMC data from 10x benchmark data set
#'
#' @format ## `pbmc_8k_msc_results`
#' 
#' 
#' A list containing the MSC results
#' \describe{
#'   \item{modules}{A list of cell clusters}
#'   \item{module.table}{A data.frame containing cell cluster statistics (sizes, compactness, intra-cluster connectivity etc).}
#'   \item{cell.network}{An igraph object containing the cell network}
#'   \item{alpha.value}{The alpha parameter value used to calculate the cluster compactness}
#'   \item{pruned.table}{Refined list of cell clusters after applying compactness and intra-cluster connectivity filters.}
#' }
#' @source <https://www.10xgenomics.com/resources/datasets/8-k-pbm-cs-from-a-healthy-donor-2-standard-2-0-1>
#' @name pbmc_8k_msc_results
"pbmc_8k_msc_results"
#' simMix1 data set from pipeComp study
#'
#' @format ## `simMix1`
#' 
#' 
#' A list containing the MSC results
#' \describe{
#'   \item{simMix}{SingleCellExperiment object}
#' }
#' @source <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02136-7>
#' @name simMix1
"simMix1"