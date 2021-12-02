#####RAT - Research Assessment Tools
#####Version 0.1.0 (2021-12-01)
#####By Pedro Cardoso & Stefano Mammola
#####Maintainer: pedro.cardoso@helsinki.fi
#####Reference: Cardoso, P., Fukushima, C.S. & Mammola, S. (subm.) Quantifying the international collaboration of researchers and research institutions.
#####Changed from v0.0.0:
#####Everything

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
#' @examples wos("A-8820-2008")
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

  id = try(wosr::pull_wos(query = query, editions = coll, sid = auth(username = user, password = pass)), silent = TRUE)

  if(class(id) == "try-error"){
    warning("Could not connect to server.")
    return()
  } else {
    return(id)
  }
}

#' H-index.
#' @description Calculates the h-index based on Web of Science data.
#' @param id A list obtained with function 'wos'.
#' @param fulldata if TRUE returns publication and citation counts.
#' @details The h-index is a measure of scientific output calculated as the h number of papers with more than h citations (Hirsch, 2005).
#' @return The h-index value. If fulldata = TRUE a list with full data.
#' @references Hirsch, J.E. (2005). An index to quantify an individual's scientific research output. PNAS, 102: 16569–16572. doi:10.1073/pnas.0507655102.
#' @examples id = wos("A-8820-2008")
#' h.index(id)
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
#' @param fulldata if TRUE returns publication and citation counts.
#' @details The i-index is a measure of scientific collaborations across countries. Calculated as the i number of co-author countries in more than i papers (Cardoso et al. subm.).
#' @return The i-index value. If fulldata = TRUE a list with full data.
#' @references Cardoso, P., Fukushima, C.S. & Mammola, S. (subm.) Quantifying the international collaboration attitude of scholars.
#' @examples stefano = wos("I-1518-2019")
#' i.index(stefano)
#' pedro = wos("A-8820-2008")
#' i.index(pedro, fulldata = TRUE)
#' @export
i.index <- function(id, fulldata = FALSE){

  if(is.null(id)) return()

  #get country data from each paper
  countries = unique(id$address[, c('ut', 'country')])
  countries = table(countries$country)
  countries = countries[order(countries, decreasing = TRUE)]

  #calculate i as the i countries in more than i papers
  i = length(which(countries >= 1:length(countries)))

  #extra data
  if(fulldata){
    n_countries = length(countries)
    i = list(i_index = i, n_countries = n_countries, countries = countries)
  }

  return(i)
}

#' Map of international collaboration.
#' @description Generates a network of international collaboration.
#' @param id A list obtained with function 'wos'.
#' @param homeCountry A character string specifying the country of origin of the researcher. Look at map$country for the complete list. If NULL, the country with most hits in Web of Science is used.
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
#' @examples id = wos("A-8820-2008")
#' i.map(id, country.size.proportional = TRUE)
#' @export
i.map <- function(id, homeCountry = NULL,
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
  countries   <- map[map[,1] %in% names(count),]
  homeCountry <- countries[countries$country %in% homeCountry,]

  #Set size for countries
  if(country.size.proportional == TRUE){
    country.point.size = count[sort(names(count))]
    country.point.size = country.point.size[names(country.point.size) %in% countries$country]

  }
  if(is.null(homeCountry.point.size))
    homeCountry.point.size = sqrt(max(country.point.size, na.rm = TRUE))

  #Make the plot
  net <- ggplot() +
    #plot map
    geom_map(map = world, data = world,
             aes(map_id = region),
             color = country.border.col,
             fill = country.col,
             size = country.border.tick) +

    #plot line
    geom_curve(data = countries, aes(x = homeCountry$x, #jitter to avoid problems with identical points (if any)
                                     y = homeCountry$y,
                                     xend = jitter(x,0.000001),
                                     yend = jitter(y,0.000001)),
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
      plot.background = element_rect(fill = sea.col, colour = sea.col))

  return(net)
}

#' Matrix matching country names and coordinates.
#'
#' A dataset that links author countries with the map using the coordinates.
#'
#' @docType data
#' @keywords datasets
#' @name map
#' @usage data(map)
#' @format A matrix with countries and corresponding coordinates.
NULL
