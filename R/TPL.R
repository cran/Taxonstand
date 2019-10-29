TPL <-
function(splist, genus = NULL, species = NULL, infrasp = NULL, infra = TRUE, corr = TRUE, diffchar = 2, max.distance = 1, version = "1.1", encoding = "UTF-8", author = TRUE, drop.lower.level = FALSE, file = "", silent = TRUE, repeats = 6) {

  splist2 <- NULL

  try(splist2 <- splist, silent = TRUE)

  if (!is.null(splist2) && (!is.null(genus) || !is.null(species) || !is.null(infrasp))) {

    stop("Argument 'splist' incompatible with arguments 'genus' and 'species'")

  } else if (is.null(splist2) && ((is.null(genus) && !is.null(species)) || (!is.null(genus) && is.null(species)))) {

    stop("Arguments 'genus' and 'species' must be provided")

  } else if (is.null(splist2) && !is.null(genus) && !is.null(species)) {

    if (infra == TRUE && !is.null(infrasp)) {

      splist <- paste(genus, species, infrasp)

    } else if (infra == FALSE || is.null(infrasp)) {

      splist <- paste(genus, species)

    }

  }



  TPLck2 <- function(d) {

    a <- NULL

    if (silent == FALSE) {

      print(paste("Checking", as.character(d), "in The Plant List"))

    }

    counter <- 0

    a <- try(TPLck(sp = d, infra = infra, corr = corr, diffchar = diffchar, max.distance = max.distance, version = version, encoding = encoding, author = author, drop.lower.level = drop.lower.level), silent = FALSE)

    while(class(a) == "try-error" && counter < repeats) {

      a <- try(TPLck(sp = d, infra = infra, corr = corr, diffchar = diffchar, max.distance = max.distance, version = version, encoding = encoding, author = author, drop.lower.level = drop.lower.level), silent = TRUE)

      counter <- counter + 1

    }

    if (class(a) == "try-error") {

      # a <- NULL

      a <- data.frame(Taxon = d, Genus = NA, Hybrid.marker = NA, Species = NA, Abbrev = NA, Infraspecific.rank = NA, Infraspecific = NA,
                      Authority = NA, ID = NA, Plant.Name.Index = NA, TPL.version = NA, Taxonomic.status = NA,
                      Family = NA, New.Genus = NA, New.Hybrid.marker = NA, New.Species = NA, New.Infraspecific.rank = NA, New.Infraspecific = NA,
                      New.Authority = NA, New.ID = NA, New.Taxonomic.status = NA, Typo = NA, WFormat = NA, Higher.level = NA, Date = NA, Tax_res = NA, 
                      stringsAsFactors = FALSE)

    }

    invisible(a)

  }



  if (length(splist) < 5) {

    results <- do.call("rbind", lapply(splist, TPLck2))

  } else {

    op <- pbapply::pboptions()

    pbapply::pboptions(type = "txt")

    results <- do.call("rbind", pbapply::pblapply(splist, TPLck2))

    pbapply::pboptions(op)

  }



  if (infra == FALSE) {

    results <- results[, -c(4, 13)]

  }

  if (nchar(file) > 0) {

    write.csv(results, file = file, row.names = FALSE)

  }

  return(results)

}
