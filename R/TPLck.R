TPLck <-
function(sp, corr=FALSE, diffchar=2, max.distance=1, infra=TRUE, version="1.1") {
sp <- as.character(sp)
genus <- unlist(strsplit(sp," "))[1]
species <- unlist(strsplit(sp," "))[2]
infrasp <- ifelse(length(unlist(strsplit(sp," "))) > 2, unlist(strsplit(sp," "))[3], "")
vv <- ifelse(version == "1.0", "", version)

searchstring <- paste("http://www.theplantlist.org/tpl", vv, "/search?q=",genus,"+",species, "&csv=true", sep="")
table.sp <- NULL
try(table.sp <- read.csv(searchstring, header=TRUE, sep=",", fill=TRUE, colClasses="character", as.is = TRUE),silent=TRUE)
Genus <- genus
Species <- species
Infraspecific <- infrasp
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
if (length(unique(paste(table.sp$Genus, table.sp$Species))) > 1 && corr==TRUE) {
spx <- length(agrep(species, "sp", max.distance=0)) + length(agrep(species, "sp.", max.distance=0))
mf <- c(as.character(1:1000))
is.mf <- agrep(species, mf, max.distance=list(deletions=1), value=TRUE)
cck <- agrep(species, table.sp$Species, value=TRUE, max.distance=max.distance)
ddf <- abs(nchar(cck) - nchar(species))
if(length(cck) > 0) {
cck <- cck[ddf==min(ddf)]
ddf <- abs(nchar(cck) - nchar(species))
}
levs <- length(unique(cck))
if(length(is.mf) == 0 && length(cck) > 0 && ddf <= diffchar && levs == 1 && spx == 0) {
searchstring <- paste("http://www.theplantlist.org/tpl", vv, "/search?q=",genus,"+",cck[1], "&csv=true", sep="")
try(table.sp <- read.table(searchstring, header=TRUE, sep=",", fill=TRUE,colClasses = "character", as.is = TRUE), silent=T)
k <- dim(table.sp)[2]
z <- dim(table.sp)[1]
marker <- TRUE
}
}
# Loop 2.2.2.2
if (length(unique(paste(table.sp$Genus, table.sp$Species))) > 1) {
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
if (length(unique(paste(table.sp$Genus, table.sp$Species))) == 1) {
grep1 <- grep(infrasp, table.sp$Infraspecific.epithet, value=TRUE)
ngrep <- nchar(grep1)
# Loop 2.2.2.3.A
if (infra==TRUE && length(grep(infrasp, table.sp$Infraspecific.epithet))>0 && abs(ngrep-(nchar(infrasp)))==0) {
table.sp <- table.sp[table.sp$Infraspecific.epithet==infrasp,]
table.sp$Taxonomic.status.in.TPL <- table.sp$Taxonomic.status.in.TPL
} else
# Loop 2.2.2.3.B
if (infra==FALSE || length(grep(infrasp, table.sp$Infraspecific.epithet))==0 || abs(ngrep-(nchar(infrasp)))>0) {
# Loop 2.2.2.3.B.1
if((table.sp$Infraspecific.epithet=="" || is.na(table.sp$Infraspecific.epithet))==FALSE && length(grep(Species, table.sp$Infraspecific.epithet))>0) {
table.sp <- table.sp[table.sp$Infraspecific.epithet==Species, ]
table.sp$Taxonomic.status.in.TPL <- table.sp$Taxonomic.status.in.TPL
} else
# Loop 2.2.2.3.B.2
if((table.sp$Infraspecific.epithet=="" || is.na(table.sp$Infraspecific.epithet))==FALSE && length(grep(Species, table.sp$Infraspecific.epithet))==0) {
table.sp <- table.sp[table.sp$Infraspecific.epithet=="" || is.na(table.sp$Infraspecific.epithet),]
table.sp$Taxonomic.status.in.TPL <- table.sp$Taxonomic.status.in.TPL
}
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
if (any(table.sp$Taxonomic.status.in.TPL=="Accepted")) {
table.sp <- table.sp[table.sp$Taxonomic.status.in.TPL=="Accepted", ]
Plant.Name.Index <- TRUE
Taxonomic.status <- table.sp$Taxonomic.status.in.TPL[1]
Family <- table.sp$Family[1]
New.Genus <- table.sp$Genus[1]
New.Species <- table.sp$Species[1]
if (infra == T && length(grep(infrasp, table.sp$Infraspecific.epithet))>0) {
New.Infraspecific <- table.sp$Infraspecific.epithet[1]
} else
if (infra == F || length(grep(infrasp, table.sp$Infraspecific.epithet))==0) {
New.Infraspecific <- ""
}
Authority <- table.sp$Authorship[1]
Typo <- ifelse(marker==TRUE, TRUE, FALSE)
WFormat <- FALSE
} else
# Loop 2.2.2.3.3 (Search for synonyms)
if (sum(table.sp$Taxonomic.status.in.TPL=="Accepted")==0 && sum(table.sp$Taxonomic.status.in.TPL=="Synonym")>0) {
table.sp <- table.sp[table.sp$Taxonomic.status.in.TPL=="Synonym", ]
#table.sp.id <- table.sp$ID[[1]]
table.sp.id <- table.sp[1,1]
at <- readLines(paste("http://www.theplantlist.org/tpl", vv, "/record/", table.sp.id, sep=""))
if (version=="1.1") {
az <- "<p>This name is a <a href=\"/1.1/about/#synonym\">synonym</a> of"
} else
if (version=="1.0") {
az <- "<p>This name is a <a href=\"/about/#synonym\">synonym</a> of"
}
n <- pmatch(az, at)
nsen <- at[n]
nsen <- unlist(strsplit(unlist(strsplit(nsen, split=">")), "<"))
if (version=="1.1") {
searchstring <- paste("http://www.theplantlist.org/tpl", vv, "/search?q=",nsen[13],"+",ifelse(nsen[17]=="\u00D7", nsen[21], nsen[17]), "&csv=true", sep="")
} else
if (version=="1.0") {
searchstring <- paste("http://www.theplantlist.org/tpl", vv, "/search?q=",nsen[11],"+",ifelse(nsen[15]=="\u00D7", nsen[19], nsen[15]), "&csv=true", sep="")
}
kup <- length(grep("var.", nsen)) + length(grep("subsp.", nsen))
if(infra==T && kup > 0) {
if (version=="1.1") {
infrasp <- ifelse(nsen[17]=="\u00D7", nsen[29], nsen[25])
} else
if (version=="1.0") {
infrasp <- ifelse(nsen[17]=="\u00D7", nsen[27], nsen[23])
}
} else
if(kup == 0) {
infrasp <- ""
}
try(table.sp <- read.table(searchstring, header=TRUE, sep=",", fill=TRUE,
                                                                   colClasses="character", as.is = TRUE),silent=TRUE)
grep1 <- grep(infrasp, table.sp$Infraspecific.epithet, value=TRUE)
ngrep <- nchar(grep1)
# Loop 2.2.2.3.3.A
if (infra==TRUE && length(grep(infrasp, table.sp$Infraspecific.epithet))>0 && abs(ngrep-(nchar(infrasp)))==0) {
table.sp <- table.sp[table.sp$Infraspecific.epithet==infrasp,]
table.sp$Taxonomic.status.in.TPL <- table.sp$Taxonomic.status.in.TPL
} else
# Loop 2.2.2.3.3.B
if (infra==FALSE || length(grep(infrasp, table.sp$Infraspecific.epithet))==0 || abs(ngrep-(nchar(infrasp)))>0) {
table.sp <- table.sp[table.sp$Infraspecific.epithet=="" || is.na(table.sp$Infraspecific.epithet),]
table.sp$Taxonomic.status.in.TPL <- table.sp$Taxonomic.status.in.TPL
}
Plant.Name.Index <- TRUE
Taxonomic.status <- "Synonym"
Family <- table.sp$Family[1]
New.Genus <- table.sp$Genus[1]
New.Species <- table.sp$Species[1]
New.Infraspecific <- table.sp$Infraspecific.epithet[1]
Authority <- table.sp$Authorship[1]
Typo <- ifelse(marker==TRUE, TRUE, FALSE)
WFormat <- FALSE
} else
# Loop 2.2.2.3.4 (Search for unresolved names)
if (sum(table.sp$Taxonomic.status.in.TPL=="Accepted")==0 && sum(table.sp$Taxonomic.status.in.TPL=="Synonym")==0 && sum(table.sp$Taxonomic.status.in.TPL=="Unresolved")>0) {
Plant.Name.Index <- TRUE
Taxonomic.status <- "Unresolved"
Family <- table.sp$Family[1]
New.Genus <- table.sp$Genus[1]
New.Species <- table.sp$Species[1]
New.Infraspecific <- table.sp$Infraspecific.epithet[1]
Authority <- table.sp$Authorship[1]
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
if (table.sp$Taxonomic.status.in.TPL=="Synonym"||table.sp$Taxonomic.status.in.TPL=="Misapplied") {
at <- readLines(paste("http://www.theplantlist.org/tpl", vv, "/search?q=", genus,"+", species, sep=""))
if (table.sp$Taxonomic.status.in.TPL=="Synonym") {
if (version=="1.1") {
az <- "<p>This name is a <a href=\"/1.1/about/#synonym\">synonym</a> of"
} else
if (version=="1.0") {
az <- "<p>This name is a <a href=\"/about/#synonym\">synonym</a> of"
}
} else 
if (table.sp$Taxonomic.status.in.TPL=="Misapplied") {
if (version=="1.1") {
az <- "<p>In the past this name has been <a href=\"/1.1/about/#misapplied\">erroneously used</a> to refer to"
} else 
if (version=="1.0") {
az <- "<p>In the past this name has been <a href=\"/about/#misapplied\">erroneously used</a> to refer to"
}
}
n <- pmatch(az, at)
nsen <- at[n]
nsen <- unlist(strsplit(unlist(strsplit(nsen, split=">")), "<"))
if (version=="1.1") {
searchstring <- paste("http://www.theplantlist.org/tpl", vv, "/search?q=",nsen[13],"+",ifelse(nsen[17]=="\u00D7", nsen[21], nsen[17]), "&csv=true", sep="")
} else
if (version=="1.0") {
searchstring <- paste("http://www.theplantlist.org/tpl", vv, "/search?q=",nsen[11],"+",ifelse(nsen[15]=="\u00D7", nsen[19], nsen[15]), "&csv=true", sep="")
}
kup <- length(grep("var.", nsen)) + length(grep("subsp.", nsen))
if(infra==T && kup > 0) {
if (version=="1.1") {
infrasp <- ifelse(nsen[17]=="\u00D7", nsen[29], nsen[25])
} else
if (version=="1.0") {
infrasp <- ifelse(nsen[17]=="\u00D7", nsen[27], nsen[23])
}
} else
if(kup == 0) {
infrasp <- ""
}
if (table.sp$Taxonomic.status.in.TPL=="Synonym") {
Taxonomic.status <- "Synonym"
} else 
if (table.sp$Taxonomic.status.in.TPL=="Misapplied") {
Taxonomic.status <- "Misapplied"
}
try(table.sp <- read.table(searchstring, header=TRUE, sep=",", fill=TRUE, colClasses="character"),silent=TRUE)
Plant.Name.Index <- TRUE
Family <- table.sp$Family[1]
if (version=="1.1") {
New.Genus <- nsen[13]
New.Species <- ifelse(nsen[17]=="\u00D7", nsen[21], nsen[17])
Authority <- ifelse(nsen[17]=="\u00D7", nsen[25], nsen[21])
} else
if (version=="1.0") {
New.Genus <- nsen[11]
New.Species <- ifelse(nsen[15]=="\u00D7", nsen[19], nsen[15])
Authority <- ifelse(nsen[15]=="\u00D7", nsen[23], nsen[19])
}
grep1 <- grep(infrasp, table.sp$Infraspecific.epithet, value=TRUE)
ngrep <- nchar(grep1)
# Loop 2.2.2.3.3.A
if (infra==TRUE && length(grep(infrasp, table.sp$Infraspecific.epithet))>0 && abs(ngrep-(nchar(infrasp)))==0) {
table.sp <- table.sp[table.sp$Infraspecific.epithet==infrasp,]
table.sp$Taxonomic.status.in.TPL <- table.sp$Taxonomic.status.in.TPL
} else
# Loop 2.2.2.3.3.B
if (infra==FALSE || length(grep(infrasp, table.sp$Infraspecific.epithet))==0 || abs(ngrep-(nchar(infrasp)))>0) {
table.sp <- table.sp[table.sp$Infraspecific.epithet=="" || is.na(table.sp$Infraspecific.epithet),]
table.sp$Taxonomic.status.in.TPL <- table.sp$Taxonomic.status.in.TPL
}
New.Infraspecific <- table.sp$Infraspecific.epithet[1]
Typo <- ifelse(marker==TRUE, TRUE, FALSE)
WFormat <- FALSE
} else
# Loop 2.2.3.3
if (table.sp$Taxonomic.status.in.TPL=="Accepted"||table.sp$Taxonomic.status.in.TPL=="Unresolved") {
Plant.Name.Index <- TRUE
Taxonomic.status <- table.sp$Taxonomic.status.in.TPL[1]
Family <- table.sp$Family
New.Genus <- table.sp$Genus
New.Species <- table.sp$Species
New.Infraspecific <- table.sp$Infraspecific.epithet
Authority <- table.sp$Authorship
Typo <- ifelse(marker==TRUE, TRUE, FALSE)
WFormat <- FALSE
}
} # End of loop 2.2.3
} # End of loop 2.1
} # End of loop 2
results <- data.frame(Genus, Species, Infraspecific, Plant.Name.Index, TPL_version=version, Taxonomic.status, Family, New.Genus, New.Species, New.Infraspecific, Authority, Typo, WFormat, stringsAsFactors = FALSE)
}
