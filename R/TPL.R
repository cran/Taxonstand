TPL <-
function(splist, genus=NULL, species=NULL, infrasp=NULL, infra=TRUE, abbrev=TRUE, corr=FALSE, diffchar=2, max.distance=1, file="") {
	splist2 <- NULL
	try(splist2 <- splist, silent=TRUE)
	if(!is.null(splist2) && (!is.null(genus) || !is.null(species) || !is.null(infrasp))) {
		stop("argument 'splist' incompatible with arguments 'genus' and 'species'")
		} else
	if(is.null(splist2) && ((is.null(genus) && !is.null(species)) || (!is.null(genus) && is.null(species)))) {
		stop("arguments 'genus' and 'species' must be provided")
		} else
	if(!is.null(splist2) && is.null(genus) && is.null(species) && is.null(infrasp)) {
		if(abbrev==TRUE) {
		abbr <- rep(NA, length(splist))
		splist0 <- splist
		for(i in 1:length(splist)) {abbr[i] <- ifelse(length(grep("cf. ", splist[i], fixed=TRUE))>0, "cf.", abbr[i])}
		splist <- sub("cf. ", "", splist, fixed=TRUE)
		for(i in 1:length(splist)) {abbr[i] <- ifelse(length(grep("cf.", splist[i], fixed=TRUE))>0, "cf.", abbr[i])}
		splist <- sub("cf.", "", splist, fixed=TRUE)
		for(i in 1:length(splist)) {abbr[i] <- ifelse(length(grep("cf ", splist[i], fixed=TRUE))>0, "cf.", abbr[i])}
		splist <- sub("cf ", "", splist, fixed=TRUE)
		for(i in 1:length(splist)) {abbr[i] <- ifelse(length(grep("aff. ", splist[i], fixed=TRUE))>0, "aff.", abbr[i])}
		splist <- sub("aff. ", "", splist, fixed=TRUE)
		for(i in 1:length(splist)) {abbr[i] <- ifelse(length(grep("aff.", splist[i], fixed=TRUE))>0, "aff.", abbr[i])}
		splist <- sub("aff.", "", splist, fixed=TRUE)
		for(i in 1:length(splist)) {abbr[i] <- ifelse(length(grep("aff ", splist[i], fixed=TRUE))>0, "aff.", abbr[i])}
		splist <- sub("aff ", "", splist, fixed=TRUE)
		for(i in 1:length(splist)) {abbr[i] <- ifelse(length(grep("s.l. ", splist[i], fixed=TRUE))>0, "s.l.", abbr[i])}
		splist <- sub("s.l. ", "", splist, fixed=TRUE)
		for(i in 1:length(splist)) {abbr[i] <- ifelse(length(grep("s.l.", splist[i], fixed=TRUE))>0, "s.l.", abbr[i])}
		splist <- sub("s.l.", "", splist, fixed=TRUE)
		for(i in 1:length(splist)) {abbr[i] <- ifelse(length(grep("s.l ", splist[i], fixed=TRUE))>0, "s.l.", abbr[i])}
		splist <- sub("s.l ", "", splist, fixed=TRUE)
		for(i in 1:length(splist)) {abbr[i] <- ifelse(length(grep("s.str. ", splist[i], fixed=TRUE))>0, "s.str.", abbr[i])}
		splist <- sub("s.str. ", "", splist, fixed=TRUE)
		for(i in 1:length(splist)) {abbr[i] <- ifelse(length(grep("s.str.", splist[i], fixed=TRUE))>0, "s.str.", abbr[i])}
		splist <- sub("s.str.", "", splist, fixed=TRUE)
		for(i in 1:length(splist)) {abbr[i] <- ifelse(length(grep("s.str ", splist[i], fixed=TRUE))>0, "s.str.", abbr[i])}
		splist <- sub("s.str ", "", splist, fixed=TRUE)
		ns <- sum(is.na(match(splist0,splist)))
		print(paste(ns, "substitutions of standard annotations in specific epithet"))
		splist <- sub("subsp. ", "", splist, fixed=TRUE)		
		splist <- sub("subsp.", "", splist, fixed=TRUE)
		splist <- sub("subsp ", "", splist, fixed=TRUE)		
		splist <- sub("var. ", "", splist, fixed=TRUE)
		splist <- sub("var.", "", splist, fixed=TRUE)
		splist <- sub("var ", "", splist, fixed=TRUE)
		ns <- sum(is.na(match(splist0,splist)))
		print(paste(ns, "substitutions of subsp. and var. in infraspecific epithet"))
		}
		splist <- gsub("      ", " ", splist, fixed=TRUE)
		splist <- gsub("     ", " ", splist, fixed=TRUE)
		splist <- gsub("    ", " ", splist, fixed=TRUE)
		splist <- gsub("   ", " ", splist, fixed=TRUE)
		splist <- gsub("  ", " ", splist, fixed=TRUE)
		f<-function(x) {unlist(strsplit(x," "))[1]}
		genus<-unlist(lapply(splist,f))
		f<-function(x) {unlist(strsplit(x," "))[2]}
		species<-unlist(lapply(splist,f))
		f<-function(x) {unlist(strsplit(x," "))[3]}
		infrasp<-unlist(lapply(splist,f))
	} else
	if(is.null(splist2) && (!is.null(genus) && !is.null(species))) {
		if(length(genus) != length(species)) {
			stop("'genus' and 'species' are not vectors of the same length")
		}
		if(abbrev==TRUE) {
		species0 <- species
		abbr <- rep(NA, length(species))
		for(i in 1:length(species)) {abbr[i] <- ifelse(length(grep("cf. ", species[i], fixed=TRUE))>0, "cf.", abbr[i])}
		species <- sub("cf. ", "", species, fixed=TRUE)		
		for(i in 1:length(species)) {abbr[i] <- ifelse(length(grep("cf.", species[i], fixed=TRUE))>0, "cf.", abbr[i])}
		species <- sub("cf.", "", species, fixed=TRUE)
		for(i in 1:length(species)) {abbr[i] <- ifelse(length(grep("cf ", species[i], fixed=TRUE))>0, "cf.", abbr[i])}
		species <- sub("cf ", "", species, fixed=TRUE)
		for(i in 1:length(species)) {abbr[i] <- ifelse(length(grep("aff. ", species[i], fixed=TRUE))>0, "aff.", abbr[i])}
		species <- sub("aff. ", "", species, fixed=TRUE)
		for(i in 1:length(species)) {abbr[i] <- ifelse(length(grep("aff.", species[i], fixed=TRUE))>0, "aff.", abbr[i])}
		species <- sub("aff.", "", species, fixed=TRUE)
		for(i in 1:length(species)) {abbr[i] <- ifelse(length(grep("aff ", species[i], fixed=TRUE))>0, "aff.", abbr[i])}
		species <- sub("aff ", "", species, fixed=TRUE)
		for(i in 1:length(species)) {abbr[i] <- ifelse(length(grep("s.l. ", species[i], fixed=TRUE))>0, "s.l.", abbr[i])}
		species <- sub("s.l. ", "", species, fixed=TRUE)
		for(i in 1:length(species)) {abbr[i] <- ifelse(length(grep("s.l.", species[i], fixed=TRUE))>0, "s.l.", abbr[i])}
		species <- sub("s.l.", "", species, fixed=TRUE)
		for(i in 1:length(species)) {abbr[i] <- ifelse(length(grep("s.l ", species[i], fixed=TRUE))>0, "s.l.", abbr[i])}
		species <- sub("s.l ", "", species, fixed=TRUE)
		for(i in 1:length(species)) {abbr[i] <- ifelse(length(grep("s.str. ", species[i], fixed=TRUE))>0, "s.str.", abbr[i])}		
		species <- sub("s.str. ", "", species, fixed=TRUE)
		for(i in 1:length(species)) {abbr[i] <- ifelse(length(grep("s.str.", species[i], fixed=TRUE))>0, "s.str.", abbr[i])}		
		species <- sub("s.str.", "", species, fixed=TRUE)
		for(i in 1:length(species)) {abbr[i] <- ifelse(length(grep("s.str ", species[i], fixed=TRUE))>0, "s.str.", abbr[i])}		
		species <- sub("s.str ", "", species, fixed=TRUE)
		ns <- sum(is.na(match(species0,species)))
		print(paste(ns, "substitutions of standard annotations in specific epithet"))
		}
		genus <- gsub("  ", "", genus, fixed=TRUE)
		genus <- gsub(" ", "", genus, fixed=TRUE)
		species <- gsub("  ", "", species, fixed=TRUE)
		species <- gsub(" ", "", species, fixed=TRUE)
		infrasp <- gsub("  ", "", infrasp, fixed=TRUE)
		infrasp <- gsub(" ", "", infrasp, fixed=TRUE)
	}
d <- paste(genus, species, infrasp)

TPLck2 <- function(d) {
	TPLck(sp=d, corr=corr, diffchar=diffchar, max.distance=max.distance, infra=infra)
	}
results <- do.call("rbind", lapply(d, TPLck2))
results$Infraspecific <- ifelse(results$Infraspecific=="NA", "", as.character(results$Infraspecific))
results <- data.frame(results[,1:2], Abbrev=as.character(abbr), results[,-c(1:2)])
if(infra==FALSE) {
	results <- results[,-c(4,10)]
	}
if(nchar(file)>0) {
	write.csv(results, file=file, row.names=F)
	}
return(results)
}
