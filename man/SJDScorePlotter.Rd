% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SJDScorePlotter.R
\name{SJDScorePlotter}
\alias{SJDScorePlotter}
\title{Plot SJD score}
\usage{
SJDScorePlotter(
  SJDalg,
  scores,
  lbb,
  info,
  SampleMetaNamesTable,
  clrs2end = c("plum", "purple", "blue", "blue4", "black", "darkred", "red", "orange",
    "yellow"),
  clrs1end = c("black", "black", "black", "darkred", "red", "orange", "yellow")
)
}
\arguments{
\item{SJDalg}{SJD algorithm to plot i.e 'twoStageLCA'}

\item{scores}{score list of the SJD algorithm i.e twoStageLCA$score_list}

\item{lbb}{dataset label i.e 'NeuroGenesis4'}

\item{info}{list of sample meta data matrices}

\item{SampleMetaNamesTable}{dataframe containing column information of each sample meta data matrices}

\item{clrs2end}{color scale for result scores from other algorithms. Default: c("plum","purple","blue","blue4","black","darkred","red","orange","yellow")}

\item{clrs1end}{color scale for result scores from sepNMF, concatNMF and jointNMF algorithms. Default: c("black","black","black","darkred","red","orange","yellow")}
}
\value{
A list containing ggplot object
}
\description{
plot dimensionality reduction scores for each SJD algorithm for dataset analyzed by SJD
}
\examples{

library(ggplot2)

data(NeuroGenesis4.afterWrap)
data(NeuroGenesis4.info)

SampleMetaNamesTable = data.frame(
   row.names = names(NeuroGenesis4),
   Type = c('Yaxis','Yaxis','2Dscatter','2Dscatter'),
   XaxisColumn = c("X","DAYx","tSNE_1","tsne1:ch1"),
   YaxisColumn = c("PJDscores","PJDscores","tSNE_2","tsne2:ch1"),
   COLaxisColumn = c("color","colorBYlabelsX","PJDscores","PJDscores"),
   PCHColumn = c("","","",""),
   inset = c(TRUE, TRUE, TRUE, TRUE),
   insetLOC = c("topright", "topright", "topright", "topright"),
   insetZoom = c(0.3, 0.3, 0.3, 0.3),
   ordDECREASE=c(FALSE, FALSE, FALSE, FALSE),
   CLRfoldPRB=c(0.5, 0.5, 0.5, 0.5)
)

grp = list(
Shared.All.4 = c(1 : 4),
Shared.bulk.2 = c(1, 2),
Shared.sc.2 = c(3, 4),
Hs.Meisnr.1 = c(1),
Hs.AZ.1 = c(2),
Gesch.1 = c(3),
Telley.1 = c(4)
)

dims = c(2, 2, 2, 2, 2, 2, 2)

twoStageLCA.out = twoStageLCA(dataset = NeuroGenesis4.afterWrap, group = grp, comp_num = dims)

SJDScorePlotter.obj = SJDScorePlotter(
    SJDalg = "twoStageLCA",
    scores = twoStageLCA.out$score_list,
    lbb = "NeuroGenesis4.p2",
    info = NeuroGenesis4.info,
    SampleMetaNamesTable = SampleMetaNamesTable
)

}
\keyword{images}
