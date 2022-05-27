#####RAT - Research Assessment Tools
#####Version 0.3.0 (2022-05-27)
#####By Pedro Cardoso & Stefano Mammola
#####Maintainer: pedro.cardoso@helsinki.fi
#####Reference: Cardoso, P., Fukushima, C.S. & Mammola, S. (subm.) Quantifying the internationalization and representativeness of research.
#####Changed from v0.2.0:
#####Deprecated function wos and built own code from WoS export

#####required packages
library("ggplot2")
library("graphics")
library("mapproj")
library("stats")
library("stringr")
library("utils")
#' @import ggplot2
#' @import graphics
#' @import mapproj
#' @import stats
#' @import stringr
#' @import utils

globalVariables(c("map", "world", "region", "x", "y"))

#Function to standardize country names
stdCountries <- function(countries){
  for(i in 1:length(countries)){

    #try to find country in map$country
    if(countries[i] %in% map$country){
      countries[i] = map[which(map$country == countries[i]), ]$stdCountry

    #if not found identify best match using fuzzy matching (Levenshtein edit distance)
    } else {
      d = adist(map$country, countries[i])
      countries[i] = map$country[which(d == min(d))][1]
    }
  }

  return(countries)
}

#Function to get country counts
getCountries <- function(biblio){

  if(!("C1" %in% colnames(biblio)))
    stop("biblio should be a tab-delimited file from WoS, exported using the option export as full record. Are you importing a WoS file with a different format?")

  #prepare data
  data(map, envir = environment())
  countries = matrix(NA, nrow = nrow(biblio), ncol = nrow(map))
  colnames(countries) = map$country

  #get country count
  for(i in 1:nrow(countries))
    countries[i, ] = str_detect(biblio$C1[i], map$country)
  countries = ifelse(countries, 1, 0)
  colnames(countries) = map$stdCountry
  countries = t(rowsum(t(countries), colnames(countries)))
  countries[countries > 1] = 1
  countries = colSums(countries)
  countries = countries[order(countries, decreasing = TRUE)]
  countries = countries[which(countries > 0)]

  return(countries)
}

################################################################################
################################MAIN FUNCTIONS##################################
################################################################################

#' H-index.
#' @description Calculates the h-index.
#' @param biblio A data.frame exported from Web of Science as tab delimited text, full record.
#' @param fulldata if TRUE returns publication and citation counts.
#' @details The h-index is a measure of scientific output calculated as the h number of papers with more than h citations (Hirsch, 2005).
#' @return The h-index value. If fulldata = TRUE a list with full data.
#' @references Hirsch, J.E. (2005). An index to quantify an individual's scientific research output. PNAS, 102: 16569–16572. doi:10.1073/pnas.0507655102.
#' @examples data(biblio)
#' h.index(biblio)
#' h.index(biblio, TRUE)
#' @export
h.index <- function(biblio, fulldata = FALSE){

  countries = getCountries(biblio)

  #order biblio
  biblio = biblio[order(biblio$TC, decreasing = TRUE), ]

  #get citation data from each paper
  citations = biblio$TC

  #calculate h as the h papers with more than h citations
  h = length(which(citations >= 1:length(citations)))

  #extra data
  if(fulldata){
    publications = biblio[, c("AU", "PY", "TI", "TC")]
    colnames(publications) = c("Authors", "Year", "Title", "Citations")
    publicationCount = nrow(publications)
    citations = sum(publications$Citations)
    h = list(h_index = h, n_publications = publicationCount, citations = citations, publications = publications)
  }

  return(h)
}

#' I-index.
#' @description Calculates the i-index (internationalization).
#' @param biblio A data.frame exported from Web of Science as tab delimited text, full record.
#' @param r if TRUE the i-index is multiplied by the r-index, i.e., weighted according to the expected distribution of GDP values of collaborating countries.
#' @param h if TRUE the i-index is divided by the h-index to create a measure independent of the latter.
#' @param homeCountry A character string specifying the country of origin of the researcher to calculate the r-index if r = TRUE. Look at map$country for the complete list. If NULL, the country with most hits in Web of Science is used.
#' @param logbase The log base for building the octaves of the r-index if r = TRUE.
#' @param fulldata if TRUE returns publication and citation counts.
#' @details The i-index (internationalization) is a measure of scientific collaborations across countries. Calculated as the i number of co-author countries in more than i papers (Cardoso et al. subm.).
#' The weighted version of the index multiplies its raw value by the square rooted difference between observed and expected distribution of GDP per capita of countries constituting the index (function RAT::represent).
#' The standardized distribution divides the i-index (weighted or not) by the h-index as these two are usually correlated.
#' @return The i-index value. If fulldata = TRUE a list with full data.
#' @references Cardoso, P., Fukushima, C.S. & Mammola, S. (subm.) Quantifying the internationalization and representativeness of research.
#' @examples data(biblio)
#' i.index(biblio)
#' i.index(biblio, r = TRUE, fulldata = TRUE)
#' i.index(biblio, r = TRUE, h = TRUE, logbase = 10, fulldata = TRUE)
#' @export
i.index <- function(biblio, r = FALSE, h = FALSE, homeCountry = NULL, logbase = 2, fulldata = FALSE){

  countries = getCountries(biblio)

  #calculate i as the i countries in more than i papers
  i = length(which(countries >= 1:length(countries)))

  #weight by GDP distribution
  if(r)
    i = i * r.index(biblio, homeCountry, logbase)

  #standardize by h-index
  if(h)
    i = i / h.index(biblio)

  #extra data
  if(fulldata){
    n_countries = length(countries)
    i = list(i_index = i, n_countries = n_countries, countries = countries)
  }

  return(i)
}

#' R-index.
#' @description Calculates the r-index (representativeness).
#' @param biblio A data.frame exported from Web of Science as tab delimited text, full record.
#' @param homeCountry A character string specifying the country of origin of the researcher. Look at map$country for the complete list. If NULL, the country with most hits in Web of Science is used.
#' @param logbase The log base for building the octaves.
#' @param plot plots the expected and observed distribution of collaborations according to GDP.
#' @details The r-index (representativeness) is a measure of the overlap between observed and expected distributions of GDP per capita of collaborating countries (Cardoso et al. subm.).
#' The abundance distribution of log(GDP per capita) of countries in the collaborators list is calculated (using octaves). This is compared with the global distribution of GDPs by using the overlap of both lists.
#' @return The r-index value.
#' @references Cardoso, P., Fukushima, C.S. & Mammola, S. (subm.) Quantifying the internationalization and representativeness of research.
#' @examples data(biblio)
#' r.index(biblio)
#' r.index(biblio, plot = TRUE)
#' @export
r.index <- function(biblio, homeCountry = NULL, logbase = 2, plot = FALSE){

  countries = getCountries(biblio)

  #exclude home country
  if(is.null(homeCountry))
    homeCountry = names(countries)[1] #take the country with most matches
  countries = countries[-which(names(countries) == homeCountry)]

  #return 0 if no countries beyond own country
  if (length(countries) == 0)
    return(0)

  #observed GDP distribution by octaves
  gdp = c()
  for(i in 1:length(countries))
    gdp = c(gdp, rep(as.numeric(map[which(map$stdCountry == names(countries[i]))[1], 3]), countries[i]))
  gdp = as.integer(log(gdp, logbase))
  gdp = data.frame(table(gdp)/sum(table(gdp)))

  #expected GDP distribution by octaves
  global = as.numeric(unique(map[, 2:3])$gdpPerCapita)
  global = as.integer(log(global, logbase))
  global = data.frame(table(global)/sum(table(global)))

  #calculate overlap as the sum of min values
  colnames(gdp) = colnames(global) = c("octave", "frequence")
  overlap = merge(global, gdp, by = "octave", all.x = TRUE, all.y = TRUE)
  overlap[is.na(overlap)] = 0

  if(plot)
    barplot(t(as.matrix(overlap[, 2:3])), legend.text = c("expected", "observed"), beside = TRUE)

  overlap = sum(apply(overlap[, 2:3], 1, min))

  return(overlap)
}

#' Map of international collaboration.
#' @description Generates a network of international collaboration.
#' @param biblio A data.frame exported from Web of Science as tab delimited text, full record.
#' @param homeCountry A character string specifying the country of origin of the researcher. Look at map$country for the complete list. If NULL, the country with most hits in Web of Science is used.
#' @param ext extent of the bounding box of the map in decimal degrees (minX, maxX, minY, maxY).
#' @param sea.col A character indicating the color of the sea.
#' @param country.col A character indicating the color of the countries in the world.
#' @param country.border.col A character indicating the color of the border among countries.
#' @param country.border.tick An integer value defining the size of the border line among countries.
#' @param line.curvature An integer value defining the curvature of the line connecting the home country with the countries of collaborators.
#' @param line.size An integer value defining the size of the line connecting the home country with the countries of collaborators.
#' @param line.alpha  An integer value defining the transparency of the line connecting the home country with the countries of collaborators.
#' @param line.color A character indicating the color of the line connecting the home country with the countries of collaborators.
#' @param country.point.color A character indicating the color of the vertex representing each country.
#' @param country.point.line A character indicating the color of line of the vertex representing each country.
#' @param country.point.alpha An integer value defining the transparency of the vertex representing each country.
#' @param country.size.proportional Logical. If TRUE, the size of each country is proportional to the number of collaborations.
#' @param country.point.size An integer value defining the size of vertex representing each country. Ignored if country.size.proportional = TRUE.
#' @param homeCountry.point.color A character indicating the color of the vertex representing the home country.
#' @param homeCountry.point.line A character indicating the color of the line of the vertex representing the home country.
#' @param homeCountry.point.alpha An integer value defining the transparency of the vertex representing the home country.
#' @param homeCountry.point.size An integer value defining the size of vertex representing the home country.
#' @details The network connects the researcher with all their collaborators.
#' @return A map with the network of collaborations.
#' @examples data(biblio)
#' i.map(biblio, country.size.proportional = TRUE)
#' @export
i.map <- function(biblio, homeCountry = NULL,
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
                  country.point.line  = "black",
                  country.point.alpha  = 0.8,
                  country.size.proportional = FALSE,
                  country.point.size = 1,
                  homeCountry.point.color = "darkgrey",
                  homeCountry.point.line  = "black",
                  homeCountry.point.alpha  = 0.8,
                  homeCountry.point.size = 5
){

  countries = getCountries(biblio)

  #get home country if needed
  if(is.null(homeCountry))
    homeCountry = names(countries)[1] #take the country with most matches

  #Get world data
  data(map, envir = environment())
  world <- ggplot2::map_data("world")

  #add coordinates
  map <- map[map[,1] %in% names(countries), ]
  homeCountry <- map[map$country == homeCountry, ]

  #Set size for countries
  if(country.size.proportional == TRUE)
    country.point.size = countries[order(names(countries))]
  if(is.null(homeCountry.point.size))
    homeCountry.point.size = sqrt(max(country.point.size, na.rm = TRUE))

  #convert as.numeric
  homeCountry$x <- as.numeric(homeCountry$x)
  homeCountry$y <- as.numeric(homeCountry$y)
  map$x   <- as.numeric(map$x)
  map$y   <- as.numeric(map$y)

  #Make the plot
  net <- ggplot() +
    #plot map
    geom_map(map = world, data = world,
             aes(map_id = region),
             color = country.border.col,
             fill = country.col,
             size = country.border.tick) +

    #map range
    xlim(ext[1], ext[2]) +
    ylim(ext[3], ext[4]) +

    #plot line
    geom_curve(data = map, aes(x = homeCountry$x,
                               y = homeCountry$y,
                               xend = jitter(x, 0.000001), #jitter to avoid problems with identical points (if any)
                               yend = jitter(y, 0.000001)), #jitter to avoid problems with identical points (if any)
                               curvature = line.curvature,
                               size = line.size,
                               alpha = line.alpha,
                               color = line.color) +

    #plot countries
    geom_point(data = map,
               aes(x = x, y = y),
               size = sqrt(country.point.size),
               colour = country.point.line,
               fill = country.point.color,
               alpha = country.point.alpha,
               shape = 21,
               stroke = 0.8) +

    #plot homeCountry
    geom_point(data = homeCountry,
               aes(x = x, y = y),
               colour = homeCountry.point.line,
               fill = homeCountry.point.color,
               size = homeCountry.point.size,
               alpha = homeCountry.point.alpha,
               shape = 21,
               stroke = 0.8)+

    theme_bw() +
    theme(
      axis.line = element_blank(), axis.text.x = element_blank(),
      axis.text.y = element_blank(), axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),legend.position="none",
      panel.background = element_rect(fill = sea.col, colour = sea.col),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(1,1,1,1), 'cm'),
      plot.background = element_rect(fill = sea.col, colour = sea.col))

  return(net)
}

#' Matrix matching country names, coordinates and GDP.
#'
#' A dataset that links author countries with the map using the coordinates and with GDP per capita.
#' Current GDP values are for 2020 (World Bank data: https://data.worldbank.org/indicator/NY.GDP.PCAP.PP.CD)
#'
#' @docType data
#' @keywords datasets
#' @name map
#' @usage data(map)
#' @format A data.frame with countries and corresponding coordinates.
NULL

#' biblio file for testing.
#'
#' A dataset from Web of Science, exported as tab delimited text, full record.
#'
#' @docType data
#' @keywords datasets
#' @name biblio
#' @usage data(biblio)
#' @format A data.frame with bibliographical data.
NULL

