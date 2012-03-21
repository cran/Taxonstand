TPLck <-
function(sp, corr=FALSE, diffchar=2, max.distance=1, infra=TRUE) {
	genus <- unlist(strsplit(sp," "))[1]
	species <- unlist(strsplit(sp," "))[2]
	infrasp <- ifelse(length(unlist(strsplit(sp," "))) > 2, unlist(strsplit(sp," "))[3], "")

	searchstring <- paste("http://www.theplantlist.org/tpl/search?q=",genus,"+",species, "&csv=true", sep="")
	table.sp <- NULL	
	try(table.sp <- read.csv(searchstring, header=TRUE, sep=",", fill=TRUE),silent=TRUE)
	Genus <- as.character(genus)
	Species <- as.character(species)
	Infraspecific <- as.character(infrasp)
	marker <- FALSE

	# Loop 1
	if(is.null(table.sp)) {
		Family <- ""
		Taxonomic.status <- ""
		Plant.Name.Index <- FALSE
		New.Genus <- Genus
		New.Species <- Species
		New.Infraspecific <- Infraspecific
		Authority <- ""
		Typo <- FALSE
		WFormat <- TRUE
	} else
	# Loop 2
	if(!is.null(table.sp)) {
		k <- dim(table.sp)[2]
		z <- dim(table.sp)[1]
		# Loop 2.1
		if(k == 1) {
		Family <- ""
		Plant.Name.Index <- FALSE
		Taxonomic.status <- ""
		New.Genus <- Genus
		New.Species <- Species
		New.Infraspecific <- Infraspecific
		Authority <- ""
		Typo <- FALSE
		WFormat <- FALSE
		} else
		# Loop 2.2
		if(k > 1) {
			# Loop 2.2.1
			if(z == 0) {
			Family <- ""
			Plant.Name.Index <- FALSE
			Taxonomic.status <- ""
			New.Genus <- Genus
			New.Species <- Species
			New.Infraspecific <- Infraspecific
			Authority <- ""
			Typo <- FALSE
			WFormat <- FALSE
			} else
			# Loop 2.2.2
			if(z > 1) {
				# Loop 2.2.2.1
				if (length(levels(factor(paste(table.sp$Genus, table.sp$Species)))) > 1 && corr==T) {
				spx <- length(agrep(species, "sp", max.distance=0)) + length(agrep(species, "sp.", max.distance=0))
				mf <- c(as.character(1:1000))
				is.mf <- agrep(species, mf, max.distance=list(deletions=1), value=TRUE)
				cck <- agrep(species, table.sp$Species, value=TRUE, max.distance=max.distance)
				ddf <- abs(nchar(cck) - nchar(species))
					if(length(cck) > 0) {
					cck <- cck[ddf==min(ddf)]
					ddf <- abs(nchar(cck) - nchar(species))
					}
				levs <- length(levels(factor(cck)))
					if(length(is.mf) == 0 && length(cck) > 0 && ddf <= diffchar && levs == 1 && spx == 0) {
					searchstring <- paste("http://www.theplantlist.org/tpl/search?q=",genus,"+",cck[1], "&csv=true", sep="")
					try(table.sp <- read.table(searchstring, header=TRUE, sep=",", fill=TRUE), silent=T)
					k <- dim(table.sp)[2]
					z <- dim(table.sp)[1]
					marker <- TRUE
					}
				}
				# Loop 2.2.2.2
				if (length(levels(factor(paste(table.sp$Genus, table.sp$Species)))) > 1) {
				Family <- ""
				Plant.Name.Index <- FALSE
				Taxonomic.status <- ""
				New.Genus <- Genus
				New.Species <- Species
				New.Infraspecific <- Infraspecific
				Authority <- ""
				Typo <- FALSE
				WFormat <- FALSE
				}
				# Loop 2.2.2.3
				if (length(levels(factor(paste(table.sp$Genus, table.sp$Species)))) == 1) {
				grep1 <- grep(infrasp, table.sp$Infraspecific.epithet, value=TRUE)
				ngrep <- nchar(grep1)
					# Loop 2.2.2.3.A
					if (infra==TRUE && length(grep(infrasp, table.sp$Infraspecific.epithet))>0 && abs(ngrep-(nchar(infrasp)))==0) {
					table.sp <- table.sp[table.sp$Infraspecific.epithet==infrasp,]
					table.sp$Taxonomic.status.in.TPL <- as.factor(as.character(table.sp$Taxonomic.status.in.TPL))
					} else
					# Loop 2.2.2.3.B
					if (infra==FALSE || length(grep(infrasp, table.sp$Infraspecific.epithet))==0 || abs(ngrep-(nchar(infrasp)))>0) {
					table.sp <- table.sp[table.sp$Infraspecific.epithet=="" || is.na(table.sp$Infraspecific.epithet),]
					table.sp$Taxonomic.status.in.TPL <- as.factor(as.character(table.sp$Taxonomic.status.in.TPL))
					}
					k <- dim(table.sp)[2]
					z <- dim(table.sp)[1]
					# Loop 2.2.2.3.1
					if (z == 0) {
					Family <- ""
					Plant.Name.Index <- FALSE
					Taxonomic.status <- ""
					New.Genus <- Genus
					New.Species <- Species
					New.Infraspecific <- Infraspecific
					Authority <- ""
					Typo <- ifelse(marker==TRUE, TRUE, FALSE)
					WFormat <- FALSE
					} else
					# Loop 2.2.2.3.2 (Search for accepted names)
					if (sum(levels(table.sp$Taxonomic.status.in.TPL)=="Accepted")>0) {
					table.sp <- table.sp[table.sp$Taxonomic.status.in.TPL=="Accepted", ]
					Plant.Name.Index <- TRUE
					Taxonomic.status <- as.character(table.sp$Taxonomic.status.in.TPL[1])
					Family <- as.character(table.sp$Family[1])
					New.Genus <- as.character(table.sp$Genus[1])
					New.Species <- as.character(table.sp$Species[1])
						if (infra == T && length(grep(infrasp, table.sp$Infraspecific.epithet))>0) {
						New.Infraspecific <- as.character(table.sp$Infraspecific.epithet[1])
						} else
						if (infra == F || length(grep(infrasp, table.sp$Infraspecific.epithet))==0) {
						New.Infraspecific <- ""
						}
					Authority <- as.character(table.sp$Authorship[1])
					Typo <- ifelse(marker==TRUE, TRUE, FALSE)
					WFormat <- FALSE
					} else
					# Loop 2.2.2.3.3 (Search for synonyms)
					if (sum(levels(table.sp$Taxonomic.status.in.TPL)=="Accepted")==0 && sum(levels(table.sp$Taxonomic.status.in.TPL)=="Synonym")>0) {
					table.sp.id <- as.character(table.sp$ID[[1]])
					at <- readLines(paste("http://www.theplantlist.org/tpl/record/", table.sp.id, sep=""))
					az <- "<p>This name is a <a href=\"/about/#synonym\">synonym</a> of"
					n <- pmatch(az, at)
					nsen <- at[n]
					nsen <- unlist(strsplit(unlist(strsplit(nsen, split=">")), "<"))
					searchstring <- paste("http://www.theplantlist.org/tpl/search?q=",nsen[11],"+",nsen[15], "&csv=true", sep="")
					kup <- length(grep("var.", nsen)) + length(grep("subsp.", nsen))
						if(infra==T && kup > 0) {
						infrasp <- nsen[23]
						} else
						if(kup == 0) {
						infrasp <- ""
						}
					try(table.sp <- read.table(searchstring, header=TRUE, sep=",", fill=TRUE),silent=TRUE)
					grep1 <- grep(infrasp, table.sp$Infraspecific.epithet, value=TRUE)
					ngrep <- nchar(grep1)
					# Loop 2.2.2.3.3.A
						if (infra==TRUE && length(grep(infrasp, table.sp$Infraspecific.epithet))>0 && abs(ngrep-(nchar(infrasp)))==0) {
						table.sp <- table.sp[table.sp$Infraspecific.epithet==infrasp,]
						table.sp$Taxonomic.status.in.TPL <- as.factor(as.character(table.sp$Taxonomic.status.in.TPL))
						} else
					# Loop 2.2.2.3.3.B
						if (infra==FALSE || length(grep(infrasp, table.sp$Infraspecific.epithet))==0 || abs(ngrep-(nchar(infrasp)))>0) {
						table.sp <- table.sp[table.sp$Infraspecific.epithet=="" || is.na(table.sp$Infraspecific.epithet),]
						table.sp$Taxonomic.status.in.TPL <- as.factor(as.character(table.sp$Taxonomic.status.in.TPL))
						}
					Plant.Name.Index <- TRUE
					Taxonomic.status <- "Synonym"
					Family <- as.character(table.sp$Family[1])
					New.Genus <- as.character(table.sp$Genus[1])
					New.Species <- as.character(table.sp$Species[1])
					New.Infraspecific <- as.character(table.sp$Infraspecific.epithet[1])
					Authority <- as.character(table.sp$Authorship[1])
					Typo <- ifelse(marker==TRUE, TRUE, FALSE)
					WFormat <- FALSE
					} else
					# Loop 2.2.2.3.4 (Search for unresolved names)
					if (sum(levels(table.sp$Taxonomic.status.in.TPL)=="Accepted")==0 && sum(levels(table.sp$Taxonomic.status.in.TPL)=="Synonym")==0 && sum(levels(table.sp$Taxonomic.status.in.TPL)=="Unresolved")>0) {
					Plant.Name.Index <- TRUE
					Taxonomic.status <- "Unresolved"
					Family <- as.character(table.sp$Family[1])
					New.Genus <- as.character(table.sp$Genus[1])
					New.Species <-  as.character(table.sp$Species[1])
						New.Infraspecific <- as.character(table.sp$Infraspecific.epithet[1])
					Authority <- as.character(table.sp$Authorship[1])
					Typo <- ifelse(marker==TRUE, TRUE, FALSE)
					WFormat <- FALSE
					}}} else
			# Loop 2.2.3
			if(z == 1) {
				# Loop 2.2.3.1
				if (is.na(table.sp$Taxonomic.status.in.TPL)) {
				Taxonomic.status <- ""
				Plant.Name.Index <- FALSE
				Family <- ""
				New.Genus <- Genus
				New.Species <- Species
				New.Infraspecific <- Infraspecific
				Authority <- ""
				Typo <- ifelse(marker==TRUE, TRUE, FALSE)
				WFormat <- FALSE
				} else
				# Loop 2.2.3.2
				if (table.sp$Taxonomic.status.in.TPL=="Synonym") {
				at <- readLines(paste("http://www.theplantlist.org/tpl/search?q=", genus,"+", species, sep=""))
				az <- "<p>This name is a <a href=\"/about/#synonym\">synonym</a> of"
				n <- pmatch(az, at)
				nsen <- at[n]
				nsen <- unlist(strsplit(unlist(strsplit(nsen, split=">")), "<"))
				searchstring <- paste("http://www.theplantlist.org/tpl/search?q=",nsen[11],"+",nsen[15], "&csv=true", sep="")
				kup <- length(grep("var.", nsen)) + length(grep("subsp.", nsen))
					if(infra==T && kup > 0) {
					infrasp <- nsen[23]
					} else
					if(kup == 0) {
					infrasp <- ""
					}
				try(table.sp <- read.table(searchstring, header=TRUE, sep=",", fill=TRUE),silent=TRUE)
				Plant.Name.Index <- TRUE
				Taxonomic.status <- "Synonym"
				Family <- as.character(table.sp$Family[1])
				New.Genus <- nsen[11]
				New.Species <- nsen[15]
				grep1 <- grep(infrasp, table.sp$Infraspecific.epithet, value=TRUE)
				ngrep <- nchar(grep1)
					# Loop 2.2.2.3.3.A
					if (infra==TRUE && length(grep(infrasp, table.sp$Infraspecific.epithet))>0 && abs(ngrep-(nchar(infrasp)))==0) {
					table.sp <- table.sp[table.sp$Infraspecific.epithet==infrasp,]
					table.sp$Taxonomic.status.in.TPL <- as.factor(as.character(table.sp$Taxonomic.status.in.TPL))
					} else
					# Loop 2.2.2.3.3.B
					if (infra==FALSE || length(grep(infrasp, table.sp$Infraspecific.epithet))==0 || abs(ngrep-(nchar(infrasp)))>0) {
					table.sp <- table.sp[table.sp$Infraspecific.epithet=="" || is.na(table.sp$Infraspecific.epithet),]
					table.sp$Taxonomic.status.in.TPL <- as.factor(as.character(table.sp$Taxonomic.status.in.TPL))
					}
				New.Infraspecific <- as.character(table.sp$Infraspecific.epithet[1])
				Authority <- nsen[19]
				Typo <- ifelse(marker==TRUE, TRUE, FALSE)
				WFormat <- FALSE
				} else
				# Loop 2.2.3.3
				if (table.sp$Taxonomic.status.in.TPL=="Accepted"||table.sp$Taxonomic.status.in.TPL=="Unresolved") {
				Plant.Name.Index <- TRUE
				Taxonomic.status <- as.character(table.sp$Taxonomic.status.in.TPL[1])
				Family <- as.character(table.sp$Family)
				New.Genus <- as.character(table.sp$Genus)
				New.Species <- as.character(table.sp$Species)
				New.Infraspecific <- as.character(table.sp$Infraspecific.epithet)
				Authority <- as.character(table.sp$Authorship)
				Typo <- ifelse(marker==TRUE, TRUE, FALSE)
				WFormat <- FALSE
				}
			} # End of loop 2.2.3
} # End of loop 2.1
} # End of loop 2
results <- data.frame(Genus, Species, Infraspecific, Plant.Name.Index, Taxonomic.status, Family, New.Genus, New.Species, New.Infraspecific, Authority, Typo, WFormat)
}
