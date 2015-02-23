TPLck <- function(sp, id=NULL, corr=TRUE, diffchar=2, max.distance=1, infra=TRUE, abbrev=TRUE, version="1.1", encoding="UTF-8") {

  ## Format Species Name ##

  sp <- as.character(sp)
abbr <- NA

#Remove infraspecific abbrev from species name, add stdz abbrev to final output
if(abbrev==TRUE) {
  vec0 <- c("nothossp. ", "nothossp.", " nothossp ", "nothosubsp. ",
            "nothosubsp.", " nothosubsp ", "cultivar. ", "cultivar.",
            " cultivar ", "cv. ", "cv."," cv ", "subfo. ", " subfo ", "subf. ", "subf.",
            " subf ", " subproles ", "cf. ", "cf.", " cf ", "aff. ",
            "aff.", " aff ", "s.l. ", "s.l.", "s.l ", "s.str. ",
            "s.str.", "s.str ", "×", "x. ", "x.", " x ", "X. ",
            "X.", " X ", "X ", "f. ", "f.", " f ", "fo. ", "fo.",
            " fo ", " forma ", "subvar.", " subvar ", "subvar. ", "var. ",
            "var.", " var ", "subsp. ", "subsp.", " subsp ",
            "ssp. ", "ssp.", " ssp ", " gama ", " grex ", "lus. ",
            "lus.", " lus ", " lusus ", "monstr. ", " monstr ",
            "nm. ", "nm.", " nm ", "prol. ", "prol.", " prol ",
            " proles ", " race ")
  vec1 <- c(rep('nothosubsp.', 6), rep('cv.', 6), rep('subf.',5), "subproles",
            rep('cf.', 3), rep("aff.",3), rep('s.l.', 3), rep("s.str.", 3), rep("×",8),
            rep("f.",7), rep("subvar.",3), rep("var.",3), rep('subsp.',6)," gama ",
            " grex ", rep("lus.",4), rep("monstr.",2), rep("nm.",3), rep("prol.",4),
            " race ") #Standardize abbrev

for(j in 1:length(vec0)) {
abbr <- ifelse(length(grep(vec0[j], sp, fixed=TRUE))>0, vec1[j], abbr)
sp <- sub(vec0[j], " ", sp, fixed=TRUE)
}
abbr <- gsub(" ", "", abbr, fixed=TRUE)
sp <- gsub("      ", " ", sp, fixed=TRUE)
sp <- gsub("     ", " ", sp, fixed=TRUE)
sp <- gsub("    ", " ", sp, fixed=TRUE)
sp <- gsub("   ", " ", sp, fixed=TRUE)
sp <- gsub("  ", " ", sp, fixed=TRUE)
sp <- gsub("_", "", sp, fixed=TRUE)
sp <- ifelse(substr(sp, 1, 1)==" ", substr(sp, 2, nchar(sp)), sp)
}
Abbrev <- abbr

#Format name for TPL
genus <- unlist(strsplit(sp," "))[1]
species <- unlist(strsplit(sp," "))[2]
if(corr==TRUE) {species <- tolower(species)}
infrasp <- ifelse(length(unlist(strsplit(sp," "))) > 2, unlist(strsplit(sp," "))[3], "")
if(corr==TRUE) {infrasp <- tolower(infrasp)}

## Search TPL for Species Name ##

vv <- ifelse(version == "1.0", "", version)

searchstring <- paste("http://www.theplantlist.org/tpl", vv, "/search?q=",genus,"+",species, "&csv=true", sep="")
table.sp <- NULL
try(table.sp <- read.csv(searchstring, header=TRUE, sep=",", fill=TRUE, colClasses="character", as.is = TRUE, encoding=encoding),silent=TRUE)

#Set up fields
Genus <- genus
Species <- species
Infraspecific <- infrasp
New.Hybrid.marker <- ""
marker <- FALSE
marker.infra <- FALSE
Multi <- FALSE #to mark if multiple valid names/synonyms are returned


## Function to fill in fields if species not found ##

notfound <- function(id, Genus, Species, Abbrev, Infraspecific, version,
                     New.Hybrid.marker){
  Family <- ""
  Taxonomic.status <- ""
  Plant.Name.Index <- FALSE
  New.Genus <- Genus
  New.Species <- Species
  New.Abbrev <- Abbrev
  New.Infraspecific <- Infraspecific
  Authority <- ""
  Typo <- FALSE
  WFormat <- FALSE
  ID <- ""
  New.ID <- ""
  if(is.null(id)){
    return(data.frame(Genus, Species, Abbrev, Infraspecific,
                      ID, Plant.Name.Index, TPL_version = version, Taxonomic.status,
                      Family, New.Genus, New.Hybrid.marker, New.Species, New.Abbrev,
                      New.Infraspecific, Authority, New.ID, Typo, Multi, WFormat,
                      stringsAsFactors = FALSE))
  } else {
    return(data.frame(UserID = id, Genus, Species, Abbrev, Infraspecific,
                      ID, Plant.Name.Index, TPL_version = version, Taxonomic.status,
                      Family, New.Genus, New.Hybrid.marker, New.Species, New.Abbrev,
                      New.Infraspecific, Authority, New.ID, Typo, Multi, WFormat,
                      stringsAsFactors = FALSE))
  }
}


## If Nothing is Found ##

if(is.null(table.sp)) {
  x <- notfound(id, Genus, Species, Abbrev, Infraspecific, version, New.Hybrid.marker)
  x$WFormat <- T
  return(x)
} else

if(!is.null(table.sp)) {
k <- dim(table.sp)[2]
z <- dim(table.sp)[1]

if (k == 1|z == 0|all(is.na(table.sp$Taxonomic.status.in.TPL))) {
  return(notfound(id, Genus, Species, Abbrev, Infraspecific, version,
                  New.Hybrid.marker))
}


## If Multiple Entries Are Found ##

else if (z > 0) {

  ## If Matching to Genus Only ##

  if(is.na(species)){
    if('Accepted' %in% table.sp$Taxonomic.status.in.TPL){
      table.sp <- table.sp[table.sp$Taxonomic.status.in.TPL=="Accepted",]
      #Consider Genus valid, take first accepted species to assign family, etc.

      Plant.Name.Index <- TRUE
      Taxonomic.status <- table.sp$Taxonomic.status.in.TPL[1]
      Family <- table.sp$Family[1]
      New.Genus <- table.sp$Genus[1]
      New.Hybrid.marker <- ''
      New.Species <- ''
      New.Abbrev <- ""
      New.Infraspecific <- ""
      Authority <- ''
      Typo <- FALSE
      WFormat <- FALSE
      ID <- ''
      New.ID <- ID

    } else if('Synonym' %in% table.sp$Taxonomic.status.in.TPL){

      #Find most common currently accepted genus, take family, etc from that

      if('H' %in% table.sp$Confidence.level){
        table.sp <- table.sp[table.sp$Confidence.level=="H",]
      }

      samp <- sample(1:nrow(table.sp), min(ceiling(0.5*nrow(table.sp)), 15))
      newgen <- vector()

      for(i in samp){
        at <- readLines(paste("http://www.theplantlist.org/tpl",vv, "/record/",
                              table.sp[i,1], sep = ""), encoding = encoding)
        if (version == "1.1") {
          az <- "<p>This name is a <a href=\"/1.1/about/#synonym\">synonym</a> of"
        }
        else if (version == "1.0") {
          warning('Genus-level search not implemented for version 1.0')
          return(notfound(id, Genus, Species, Abbrev, Infraspecific, version,
                          New.Hybrid.marker))
        }

        syn <- at[pmatch(az, at)] #Looks for line that matches to current name
        syn <- unlist(strsplit(unlist(strsplit(syn, split = ">")), "<"))
        newgen[[i]] <- ifelse(syn[13]=='×', paste(syn[13], syn[17], sep=""), syn[13])
      }

      newgen <- newgen[samp]

      if(length(unique(newgen))>1) {
        warning(paste('Multiple possible synonyms for',genus))
        Multi <- TRUE
        freq <- vector()

        for(i in 1:length(unique(newgen))) {
          freq[i] <- length(newgen[newgen==unique(newgen)[i]])
        }
        ng <- unique(newgen)[which(freq==max(freq))] #Use the most frequent genus name
      } else {
        ng <- unique(newgen)
      }

      searchstring <- paste("http://www.theplantlist.org/tpl",vv,
                            "/search?q=", ng,"&csv=true",sep = "")

      #Search using updated genus name
      try(table.sp <- read.table(searchstring, header =TRUE, sep = ",", fill =TRUE,
                                 colClasses = "character", as.is = TRUE,
                                 encoding = encoding), silent = TRUE)

      #Keep only accepted, high-quality names if there's a choice
      if(nrow(table.sp)>1){
        if('Accepted' %in% table.sp$Taxonomic.status.in.TPL) {
          table.sp <- table.sp[table.sp$Taxonomic.status.in.TPL=="Accepted",]
        }
        if ('H' %in% table.sp$Confidence.level) {
          table.sp <- table.sp[table.sp$Confidence.level=='H',]
        } else if ('M' %in% table.sp$Confidence.level) {
          table.sp <- table.sp[table.sp$Confidence.level=='M',]
        }
      }

      Plant.Name.Index <- TRUE
      Taxonomic.status <- 'Synonym'
      Family <- table.sp$Family[1]
      New.Genus <- ng
      New.Hybrid.marker <- ''
      New.Species <- ''
      New.Abbrev <- ""
      New.Infraspecific <- ""
      Authority <- ''
      Typo <- FALSE
      WFormat <- FALSE
      ID <- ''
      New.ID <- ID
    } else { #If all names are listed as unresolved or misapplied, return as not found
      return(notfound(id, Genus, Species, Abbrev, Infraspecific, version,
                      New.Hybrid.marker))
    }

    if(is.null(id)){
      return(data.frame(Genus, Species, Abbrev, Infraspecific,
                        ID, Plant.Name.Index, TPL_version = version, Taxonomic.status,
                        Family, New.Genus, New.Hybrid.marker, New.Species, New.Abbrev,
                        New.Infraspecific, Authority, New.ID, Typo, Multi, WFormat,
                        stringsAsFactors = FALSE))
    } else {
      return(data.frame(UserID = id, Genus, Species, Abbrev, Infraspecific,
                        ID, Plant.Name.Index, TPL_version = version, Taxonomic.status,
                        Family, New.Genus, New.Hybrid.marker, New.Species, New.Abbrev,
                        New.Infraspecific, Authority, New.ID, Typo, Multi, WFormat,
                        stringsAsFactors = FALSE))
    }
  }

  ## If Doing Fuzzy Matching ##

if (length(unique(paste(table.sp$Genus, table.sp$Species))) > 1 && corr==TRUE) {
spx <- length(agrep(species, "sp", max.distance=0)) + length(agrep(species, "sp.", max.distance=0))
mf <- c(as.character(1:1000))
is.mf <- agrep(species, mf, max.distance=list(deletions=1), value=TRUE)
#Run fuzzy match w/ set max distance between given and matched sp epithet
cck <- agrep(species, table.sp$Species, value=TRUE, max.distance=max.distance)
ddf <- abs(nchar(cck) - nchar(species))
if(length(cck) > 0) {
cck <- cck[ddf==min(ddf)]
#If mult matches, select smallest diff in number of char from original
ddf <- abs(nchar(cck) - nchar(species))
}
levs <- length(unique(cck))

## Search for Best Fuzzy Match ##

if(length(is.mf) == 0 && length(cck) > 0 && ddf <= diffchar && levs == 1 && spx == 0) {
searchstring <- paste("http://www.theplantlist.org/tpl", vv, "/search?q=",genus,"+",cck[1], "&csv=true", sep="")
try(table.sp <- read.table(searchstring, header=TRUE, sep=",", fill=TRUE,colClasses = "character", as.is = TRUE, encoding=encoding), silent=T)
k <- dim(table.sp)[2]
z <- dim(table.sp)[1]
marker <- TRUE
}
}

if (length(unique(paste(table.sp$Genus, table.sp$Species))) > 1) {
  return(notfound(id, Genus, Species, Abbrev, Infraspecific, version,
                  New.Hybrid.marker))
  #If still returning entire genus list, or corr=F, match was not found
  #   could potentially alter code here to search for 2nd best match?
}

if (length(unique(paste(table.sp$Genus, table.sp$Species))) == 1) {
grep1 <- grep(infrasp, table.sp$Infraspecific.epithet, value=TRUE)
ngrep <- nchar(grep1)

## If no Appropriate Infrasp Match, Fuzzy Match Infra Names ##

if ((length(ngrep) == 0 || abs(ngrep-nchar(infrasp)) > 0) && corr==TRUE && !is.na(infrasp) && nchar(infrasp)>0) {
cck.infra <- agrep(infrasp, table.sp$Infraspecific.epithet, value=TRUE, max.distance=max.distance)
ddf.infra <- abs(nchar(cck.infra) - nchar(infrasp))
if(length(cck.infra) > 0) {
cck.infra <- cck.infra[ddf.infra==min(ddf.infra)]
}
#Remove spelling variants (to prevent unnecessary "multi" matches)
if(length(cck.infra) > 0 && "Spelling variant" %in% table.sp[,12] &&
     sum(table.sp[,12]=="Spelling variant") < length(cck.infra)){
  cck.infra <-
    cck.infra[!cck.infra %in% table.sp[table.sp[,12] %in% 'Spelling variant',
                                       'Infraspecific.epithet']]
  #table.sp[,12] is Nomenclatural.status.from.original.data.source
}
levs <- length(unique(cck.infra))
if(length(cck.infra) > 0 && levs == 1) {
infrasp <- unique(cck.infra)
grep1 <- grep(infrasp, table.sp$Infraspecific.epithet, value=TRUE)
ngrep <- nchar(grep1)
marker.infra <- TRUE
}
}

  #Match to Infrasp Name w/ same # of char
  if (infra==TRUE && length(grep(infrasp, table.sp$Infraspecific.epithet))>0 && abs(ngrep-(nchar(infrasp)))==0) {
table.sp <- table.sp[table.sp$Infraspecific.epithet==infrasp,]
if (dim(table.sp)[1]>1 && sum(!is.na(grep(Abbrev, table.sp$Infraspecific.rank, ignore.case=TRUE)))>0) {
table.sp <- table.sp[table.sp$Infraspecific.rank==Abbrev,]
}
} else

if (infra==FALSE || length(grep(infrasp, table.sp$Infraspecific.epithet))==0 || abs(ngrep-(nchar(infrasp)))>0) {

#Or match to infrasp name w/ NA or blank
if(table.sp$Infraspecific.epithet %in% c("", NA) && length(grep(Species, table.sp$Infraspecific.epithet))==0) {
table.sp <- table.sp[table.sp$Infraspecific.epithet %in% c("", NA),]
} else

  #Or match to infrasp name = original species epithet
if(!table.sp$Infraspecific.epithet %in% c("", NA) && length(grep(Species, table.sp$Infraspecific.epithet))>0) {
table.sp <- table.sp[table.sp$Infraspecific.epithet==Species, ]
}
}
#Update k and z (now that table.sp has been reduced)
k <- dim(table.sp)[2]
z <- dim(table.sp)[1]


## Evaluate New Species Table to Choose Best Match ##

if (z == 0) {
  #Nothing was a sufficiently good match, return 'no match'
  return(notfound(id, Genus, Species, Abbrev, Infraspecific, version,
                  New.Hybrid.marker))
}

#If any name is listed as "Accepted"  or High/Med Quality, Use that name
#   (uses 1st "Accepted" Name listed)
else if (any(table.sp$Taxonomic.status.in.TPL=="Accepted")) {
table.sp <- table.sp[table.sp$Taxonomic.status.in.TPL=="Accepted", ]
if(nrow(table.sp)>1){
  if (sum(table.sp$Confidence.level == "H") > 0) {
    table.sp <- table.sp[table.sp$Confidence.level == "H", ]
  }
  else if (sum(table.sp$Confidence.level == "M") > 0) {
    table.sp <- table.sp[table.sp$Confidence.level == "M", ]
  }
  if (nrow(table.sp) > 1) {
    warning(paste(sp, "has more than one valid accepted name"))
    Multi <- TRUE
  }
}
Plant.Name.Index <- TRUE
Taxonomic.status <- table.sp$Taxonomic.status.in.TPL[1]
Family <- table.sp$Family[1]
New.Genus <- table.sp$Genus[1]
New.Hybrid.marker <- table.sp$Species.hybrid.marker[1]
New.Species <- table.sp$Species[1]
if (infra == T && length(grep(infrasp, table.sp$Infraspecific.epithet))>0) {
  New.Abbrev <- table.sp$Infraspecific.rank[1]
  New.Infraspecific <- table.sp$Infraspecific.epithet[1]
} else
if (infra == F || length(grep(infrasp, table.sp$Infraspecific.epithet))==0) {
  New.Abbrev <- ""
  New.Infraspecific <- ""
}
Authority <- table.sp$Authorship[1]
Typo <- ifelse(marker==TRUE || marker.infra==TRUE, TRUE, FALSE)
WFormat <- FALSE
ID <- table.sp[1,1]
New.ID <- ID
} else
# (Search for synonyms or missapplied names)
if (sum(table.sp$Taxonomic.status.in.TPL=="Accepted")==0 && (sum(table.sp$Taxonomic.status.in.TPL=="Synonym")>0||sum(table.sp$Taxonomic.status.in.TPL=="Misapplied")>0)) {
if(sum(table.sp$Taxonomic.status.in.TPL=="Synonym")>0) {
table.sp <- table.sp[table.sp$Taxonomic.status.in.TPL=="Synonym", ]
} else
if(sum(table.sp$Taxonomic.status.in.TPL=="Misapplied")>0) {
table.sp <- table.sp[table.sp$Taxonomic.status.in.TPL=="Misapplied", ]
}
if(sum(table.sp$Confidence.level=="H")>0) {
table.sp <- table.sp[table.sp$Confidence.level=="H", ]
} else
if(sum(table.sp$Confidence.level=="M")>0) {
table.sp <- table.sp[table.sp$Confidence.level=="M", ]
}
if(nrow(table.sp)>1) {
warning(paste(sp, "has more than one valid synonym"))
Multi <- TRUE
}

if (sum(table.sp$Taxonomic.status.in.TPL == "Synonym") > 0) {
  Taxonomic.status <- "Synonym"
}
else if (sum(table.sp$Taxonomic.status.in.TPL == "Misapplied") > 0) {
  Taxonomic.status <- "Misapplied"
}

## Search for Updated Name of Synonyms in Text ##

table.sp.id <- table.sp[1,1]
ID <- table.sp.id
at <- readLines(paste("http://www.theplantlist.org/tpl", vv, "/record/", table.sp.id, sep=""), encoding=encoding)
if (sum(table.sp$Taxonomic.status.in.TPL=="Synonym")>0) {
if (version=="1.1") {
az <- "<p>This name is a <a href=\"/1.1/about/#synonym\">synonym</a> of"
} else
if (version=="1.0") {
az <- "<p>This name is a <a href=\"/about/#synonym\">synonym</a> of"
}
} else
if (sum(table.sp$Taxonomic.status.in.TPL=="Misapplied")>0) {
if (version=="1.1") {
az <- "<p>In the past this name has been <a href=\"/1.1/about/#misapplied\">erroneously used</a> to refer to"
} else
if (version=="1.0") {
az <- "<p>In the past this name has been <a href=\"/about/#misapplied\">erroneously used</a> to refer to"
}
}
n <- pmatch(az, at) #Looks for line that matches to current name
nsen <- at[n]
nsen <- unlist(strsplit(unlist(strsplit(nsen, split=">")), "<"))
if (version=="1.1") {
searchstring <- paste("http://www.theplantlist.org/tpl", vv, "/search?q=",nsen[13],"+",ifelse(nsen[17]=="×", nsen[21], nsen[17]), "&csv=true", sep="")
} else
if (version=="1.0") {
searchstring <- paste("http://www.theplantlist.org/tpl", vv, "/search?q=",nsen[11],"+",ifelse(nsen[15]=="×", nsen[19], nsen[15]), "&csv=true", sep="")
}
kup <- length(grep("var.", nsen)) + length(grep("subsp.", nsen))
if(infra==T && kup > 0) {
if (version=="1.1") {
infrasp <- ifelse(nsen[17]=="×", nsen[29], nsen[25])
} else
if (version=="1.0") {
infrasp <- ifelse(nsen[15]=="×", nsen[27], nsen[23])
}
} else
if(kup == 0) {
infrasp <- ""
}


## Find, Search by Updated Name (for Synonyms) ##

try(table.sp <- read.table(searchstring, header=TRUE, sep=",", fill=TRUE, colClasses="character", as.is = TRUE, encoding=encoding),silent=TRUE)
if (version=="1.1" && nsen[17]=="×") {
table.sp <- table.sp[table.sp$Species.hybrid.marker=="×", ]
} else
if (version=="1.0" && nsen[15]=="×") {
table.sp <- table.sp[table.sp$Species.hybrid.marker=="×", ]
}

if (infra == TRUE && any(table.sp$Infraspecific.epithet == infrasp)) {
table.sp <- table.sp[table.sp$Infraspecific.epithet==infrasp,]
} else

if (infra == FALSE || all(table.sp$Infraspecific.epithet != infrasp)) {
  table.sp <- table.sp[table.sp$Infraspecific.epithet %in% c('',NA),]
}

if(nrow(table.sp)>1){
  if('Accepted' %in% table.sp$Taxonomic.status.in.TPL) {
    table.sp <- table.sp[table.sp$Taxonomic.status.in.TPL=="Accepted",]
  } else if ('H' %in% table.sp$Confidence.level) {
    table.sp <- table.sp[table.sp$Confidence.level=='H',]
  } else if ('M' %in% table.sp$Confidence.level) {
    table.sp <- table.sp[table.sp$Confidence.level=='M',]
  }

  if(nrow(table.sp)>1){
    warning(paste(sp, "has more than one valid updated name for synonym"))
    Multi <- TRUE
  }
}

#If still mutl entries, take first one
Plant.Name.Index <- TRUE
Family <- table.sp$Family[1]
New.Genus <- table.sp$Genus[1]
New.Hybrid.marker <- table.sp$Species.hybrid.marker[1]
New.Species <- table.sp$Species[1]
New.Abbrev <- table.sp$Infraspecific.rank[1]
New.Infraspecific <- table.sp$Infraspecific.epithet[1]
Authority <- table.sp$Authorship[1]
Typo <- ifelse(marker==TRUE || marker.infra==TRUE, TRUE, FALSE)
WFormat <- FALSE
New.ID <- table.sp[1,1]
} else
# (Search for unresolved names)
if (sum(table.sp$Taxonomic.status.in.TPL=="Accepted")==0 && sum(table.sp$Taxonomic.status.in.TPL=="Synonym")==0 && sum(table.sp$Taxonomic.status.in.TPL=="Unresolved")>0) {
Plant.Name.Index <- TRUE
Taxonomic.status <- "Unresolved"
Family <- table.sp$Family[1]
New.Genus <- table.sp$Genus[1]
New.Hybrid.marker <- table.sp$Species.hybrid.marker[1]
New.Species <- table.sp$Species[1]
New.Abbrev <- table.sp$Infraspecific.rank[1]
New.Infraspecific <- table.sp$Infraspecific.epithet[1]
Authority <- table.sp$Authorship[1]
Typo <- ifelse(marker==TRUE || marker.infra==TRUE, TRUE, FALSE)
WFormat <- FALSE
ID <- table.sp[1,1]
New.ID <- ID
}
}
if(is.null(id)){
  results <- data.frame(Genus, Species, Abbrev, Infraspecific,
                        ID, Plant.Name.Index, TPL_version = version,
                        Taxonomic.status, Family, New.Genus, New.Hybrid.marker,
                        New.Species, New.Abbrev, New.Infraspecific, Authority,
                        New.ID, Typo, Multi, WFormat, stringsAsFactors = FALSE)

} else {
  UserID <- id
  results <- data.frame(UserID, Genus, Species, Abbrev, Infraspecific,
                        ID, Plant.Name.Index, TPL_version = version,
                        Taxonomic.status, Family, New.Genus, New.Hybrid.marker,
                        New.Species, New.Abbrev, New.Infraspecific, Authority,
                        New.ID, Typo, Multi, WFormat, stringsAsFactors = FALSE)

}
return(results)

}

}
}