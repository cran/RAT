% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RAT.R
\name{i.map}
\alias{i.map}
\title{Map of international collaboration.}
\usage{
i.map(
  biblio,
  homeCountry = NULL,
  ext = c(-180, 180, -55, 90),
  sea.col = "white",
  country.col = "grey",
  country.border.col = "black",
  country.border.tick = 0.3,
  line.curvature = 0.1,
  line.size = 0.8,
  line.alpha = 0.4,
  line.color = "black",
  country.point.color = "white",
  country.point.line = "black",
  country.point.alpha = 0.8,
  country.size.proportional = FALSE,
  country.point.size = 1,
  homeCountry.point.color = "darkgrey",
  homeCountry.point.line = "black",
  homeCountry.point.alpha = 0.8,
  homeCountry.point.size = 5
)
}
\arguments{
\item{biblio}{A data.frame exported from Web of Science as tab delimited text, full record OR a vector with country frequencies where names are the country names.}

\item{homeCountry}{A character string specifying the country of origin of the researcher. Look at map$country for the complete list. If NULL, the country with most hits in Web of Science is used.}

\item{ext}{extent of the bounding box of the map in decimal degrees (minX, maxX, minY, maxY).}

\item{sea.col}{A character indicating the color of the sea.}

\item{country.col}{A character indicating the color of the countries in the world.}

\item{country.border.col}{A character indicating the color of the border among countries.}

\item{country.border.tick}{An integer value defining the size of the border line among countries.}

\item{line.curvature}{An integer value defining the curvature of the line connecting the home country with the countries of collaborators.}

\item{line.size}{An integer value defining the size of the line connecting the home country with the countries of collaborators.}

\item{line.alpha}{An integer value defining the transparency of the line connecting the home country with the countries of collaborators.}

\item{line.color}{A character indicating the color of the line connecting the home country with the countries of collaborators.}

\item{country.point.color}{A character indicating the color of the vertex representing each country.}

\item{country.point.line}{A character indicating the color of line of the vertex representing each country.}

\item{country.point.alpha}{An integer value defining the transparency of the vertex representing each country.}

\item{country.size.proportional}{Logical. If TRUE, the size of each country is proportional to the number of collaborations.}

\item{country.point.size}{An integer value defining the size of vertex representing each country. Ignored if country.size.proportional = TRUE.}

\item{homeCountry.point.color}{A character indicating the color of the vertex representing the home country.}

\item{homeCountry.point.line}{A character indicating the color of the line of the vertex representing the home country.}

\item{homeCountry.point.alpha}{An integer value defining the transparency of the vertex representing the home country.}

\item{homeCountry.point.size}{An integer value defining the size of vertex representing the home country.}
}
\value{
A map with the network of collaborations.
}
\description{
Generates a network of international collaboration.
}
\details{
The network connects the researcher with all their collaborators.
}
\examples{
data(biblio)
i.map(biblio, country.size.proportional = TRUE)

biblio = c(5, 3, 2, 1)
names(biblio) = c("Finland", "Portugal", "Brazil", "Italy")
i.map(biblio)
}
