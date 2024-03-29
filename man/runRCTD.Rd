% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Deconvolution_family.R
\name{runRCTD}
\alias{runRCTD}
\title{runRCTD cell-type deconvolution}
\usage{
runRCTD(object, ref, cell_type_var, overwrite = T, return.RCTD = F)
}
\arguments{
\item{object}{The variale in the Seurat meta.data file containing the cell type group}

\item{ref}{The reference seurat object}

\item{cell_type_var}{The variale in the Seurat meta.data file containing the cell type group}
}
\value{

}
\description{
Run the a cell type deconvolution using the SpaceXR/RCTD algorithm from SPATA objects
}
\author{
Dieter Henrik Heiland
}
