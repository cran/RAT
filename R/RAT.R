#####RAT - Research Assessment Tools
#####Version 0.2.0 (2022-04-05)
#####By Pedro Cardoso & Stefano Mammola
#####Maintainer: pedro.cardoso@helsinki.fi
#####Reference: Cardoso, P., Fukushima, C.S. & Mammola, S. (subm.) Quantifying the international collaboration of researchers and research institutions.
#####Changed from v0.1.1:
#####Added function r.index
#####Added parameters r, h, homeCountry, and logbase to i.index
#####Added parameter ext to i.map
#####Added algorithm to standardize names
#####Added algorithm to identify matches when country name is not recognized

#####required packages
library("ggplot2")
library("graphics")
library("mapproj")
library("stats")
library("utils")
library("wosr")
#' @import ggplot2
#' @import graphics
#' @import mapproj
#' @import stats
#' @import utils
#' @importFrom wosr auth pull_wos

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

################################################################################
################################MAIN FUNCTIONS##################################
################################################################################

#' Web of Science.
#' @description Downloads data from Web of Science.
#' @param id ResearcherID or name of the researcher as 'Surname, First name'. Beware that in case of multiple authors with the same name all publications will be returned. In alternative, use the search notation of Web of Science for regular searches beyond an author's name.
#' @param user username for Web of Science. Not needed if you have access through your institution.
#' @param pass password or Web of Science. Not needed if you have access through your institution.
#' @param coll Web of Science collections to include.
#' @details The data downloaded can be used in h.index and i.index functions to avoid constant downloads.
#' @return A list with publication data.
#' @examples wos("C-2482-2012")
#' @export
wos <- function(id, user = NULL, pass = NULL, coll = c("SCI", "SSCI", "AHCI", "ISTP", "ISSHP","BSCI", "BHCI", "ESCI")){

  if(!grepl("=", id)){
    if(grepl("\\d", id))
      query = paste("AI =", id)
    else
      query = paste("AU =", id)
  } else {
    query = id
  }

  id = tryCatch({
    wosr::pull_wos(query = query, editions = coll, sid = wosr::auth(username = user, password = pass), silent = TRUE)
  },
  error = function(cond) {
    message("Could not connect to server.")
  })
  return(id)
}

#' H-index.
#' @description Calculates the h-index based on Web of Science data.
#' @param id A list obtained with function 'wos'.
#' @param fulldata if TRUE returns publication and citation counts.
#' @details The h-index is a measure of scientific output calculated as the h number of papers with more than h citations (Hirsch, 2005).
#' @return The h-index value. If fulldata = TRUE a list with full data.
#' @references Hirsch, J.E. (2005). An index to quantify an individual's scientific research output. PNAS, 102: 16569–16572. doi:10.1073/pnas.0507655102.
#' @examples id = wos("C-2482-2012")
#' h.index(id)
#' h.index(id, TRUE)
#' @export
h.index <- function(id, fulldata = FALSE){

  if(is.null(id)) return()

  #get citation data from each paper
  citations = id$publication$tot_cites
  citations = citations[order(citations, decreasing = TRUE)]

  #calculate h as the h papers with more than h citations
  h = length(which(citations >= 1:length(citations)))

  #extra data
  if(fulldata){
    publications = id$publication[, c('title', 'doi', 'tot_cites')]
    colnames(publications)[3] = 'citations'
    publications = publications[order(publications$citations, decreasing = TRUE),]
    publicationCount = nrow(publications)
    citations = sum(publications$citations)
    h = list(h_index = h, n_publications = publicationCount, citations = citations, publications = publications)
  }

  return(h)
}

#' I-index.
#' @description Calculates the i-index based on Web of Science data.
#' @param id A list obtained with function 'wos'.
#' @param r if TRUE the i-index is multiplied by the r-index, i.e., weighted according to the expected distribution of GDP values of collaborating countries.
#' @param h if TRUE the i-index is divided by the h-index to create a measure independent of the latter.
#' @param homeCountry A character string specifying the country of origin of the researcher to calculate the r-index if r = TRUE. Look at map$country for the complete list. If NULL, the country with most hits in Web of Science is used.
#' @param logbase The log base for building the octaves of the r-index if r = TRUE.
#' @param fulldata if TRUE returns publication and citation counts.
#' @details The i-index (internationalization) is a measure of scientific collaborations across countries. Calculated as the i number of co-author countries in more than i papers (Cardoso et al. subm.).
#' The weighted version of the index multiplies its raw value by the square rooted difference between observed and expected distribution of GDP per capita of countries constituting the index (function RAT::represent).
#' The standardized distribution divides the i-index (weighted or not) by the h-index as these two are usually correlated.
#' @return The i-index value. If fulldata = TRUE a list with full data.
#' @references Cardoso, P., Fukushima, C.S. & Mammola, S. (subm.) Quantifying the international collaboration attitude of scholars.
#' @examples id = wos("C-2482-2012")
#' i.index(id)
#' i.index(id, r = TRUE)
#' i.index(id, r = TRUE, h = TRUE, logbase = 10, fulldata = TRUE)
#' @export
i.index <- function(id, r = FALSE, h = FALSE, homeCountry = NULL, logbase = 2, fulldata = FALSE){

  if(is.null(id)) return()

  #get country data from each paper
  countries = id$address[, c('ut', 'country')]
  countries$country = stdCountries(countries$country)
  countries = unique(countries)
  countries = table(countries$country)
  countries = countries[order(countries, decreasing = TRUE)]

  #calculate i as the i countries in more than i papers
  i = length(which(countries >= 1:length(countries)))

  #weight by GDP distribution
  if(r)
    i = i * r.index(id, homeCountry, logbase)

  #standardize by h-index
  if(h)
    i = i / h.index(id)

  #extra data
  if(fulldata){
    n_countries = length(countries)
    i = list(i_index = i, n_countries = n_countries, countries = countries)
  }

  return(i)
}

#' R-index.
#' @description Calculates the r-index based on Web of Science data.
#' @param id A list obtained with function 'wos'.
#' @param homeCountry A character string specifying the country of origin of the researcher. Look at map$country for the complete list. If NULL, the country with most hits in Web of Science is used.
#' @param logbase The log base for building the octaves.
#' @param plot plots the expected and observed distribution of collaborations according to GDP.
#' @details The r-index (representativeness) is a measure of the overlap between observed and expected distributions of GDP per capita of collaborating countries (Cardoso et al. subm.).
#' The abundance distribution of log(GDP per capita) of countries in the collaborators list is calculated (using octaves). This is compared with the global distribution of GDPs by using the overlap of both lists.
#' @return The r-index value.
#' @references Cardoso, P., Fukushima, C.S. & Mammola, S. (subm.) Quantifying the international collaboration attitude of scholars.
#' @examples id = wos("C-2482-2012")
#' r.index(id)
#' r.index(id, logbase = 10, plot = TRUE)
#' @export
r.index <- function(id, homeCountry = NULL, logbase = 2, plot = FALSE){

  if(is.null(id)) return()

  #get country data from each paper and exclude home country
  countries = id$address[, c('ut', 'country')]
  countries$country = stdCountries(countries$country)
  countries = unique(countries)$country
  if(is.null(homeCountry)){
    count = table(countries)
    count = count[order(count, decreasing = TRUE)]
    homeCountry = names(count)[1] #take the country with most matches
  }
  countries = countries[countries != homeCountry]

  #return 0 if no countries beyond own country
  if (length(countries) == 0)
    return(0)

  #observed GDP distribution by octaves
  gdp = c()
  for(i in 1:length(countries))
    gdp[i] = as.numeric(map[which(map$stdCountry == countries[i])[1], 3])
  gdp = as.integer(log(gdp, logbase))
  gdp = data.frame(table(gdp)/sum(table(gdp)))

  #global GDP distribution by octaves
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
#' @param id A list obtained with function 'wos'.
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
#' @examples id = wos("C-2482-2012")
#' i.map(id, country.size.proportional = TRUE)
#' @export
i.map <- function(id, homeCountry = NULL,
                  ext = c(-180, 180, -90, 90),
                  sea.col = "black",
                  country.col = "grey",
                  country.border.col = "black",
                  country.border.tick = 0.3,
                  line.curvature = 0.1,
                  line.size = 0.8,
                  line.alpha = 0.4,
                  line.color = "yellow",
                  country.point.color = "yellow",
                  country.point.line  = "yellow",
                  country.point.alpha  = 0.8,
                  country.size.proportional = FALSE,
                  country.point.size = 1,
                  homeCountry.point.color = "purple",
                  homeCountry.point.line  = "black",
                  homeCountry.point.alpha  = 0.8,
                  homeCountry.point.size = 5
){

  if(is.null(id)) return()

  #Get world data
  data(map, envir = environment())
  world <- ggplot2::map_data("world")

  #get country data from each paper
  count = i.index(id, fulldata = TRUE)$countries

  #estimate HomeCountry if needed
  if(is.null(homeCountry))
    homeCountry = names(count)[1] #take the country with most matches

  #add coordinates
  countries <- map[map[,1] %in% names(count), ]
  homeCountry <- countries[countries$country %in% homeCountry, ]

  #Set size for countries
  if(country.size.proportional == TRUE){
    country.point.size = count[sort(names(count))]
    country.point.size = country.point.size[names(country.point.size) %in% countries$country]
  }
  if(is.null(homeCountry.point.size))
    homeCountry.point.size = sqrt(max(country.point.size, na.rm = TRUE))

  #convert as.numeric
  homeCountry$x <- as.numeric(homeCountry$x)
  homeCountry$y <- as.numeric(homeCountry$y)
  countries$x   <- as.numeric(countries$x)
  countries$y   <- as.numeric(countries$y)

  #Make the plot
  net <- ggplot() +
    #plot map
    geom_map(map = world, data = world,
             aes(map_id = region),
             color = country.border.col,
             fill = country.col,
             size = country.border.tick) +

    #map range
    xlim(ext[1],ext[2])+
    ylim(ext[3],ext[4])+

    #plot line
    geom_curve(data = countries, aes(x = homeCountry$x,
                                     y = homeCountry$y,
                                     xend = jitter(x, 0.000001), #jitter to avoid problems with identical points (if any)
                                     yend = jitter(y, 0.000001)), #jitter to avoid problems with identical points (if any)
               curvature = line.curvature,
               size = line.size,
               alpha = line.alpha,
               color = line.color) +

    #plot countries
    geom_point(data = countries,
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
#' Current GDP values are for 2020 (World Bank data: https://data.worldbank.org/indicator/NY.GDP.PCAP.PP.CD )
#'
#' @docType data
#' @keywords datasets
#' @name map
#' @usage data(map)
#' @format A matrix with countries and corresponding coordinates.
NULL
