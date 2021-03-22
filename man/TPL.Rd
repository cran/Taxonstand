\name{TPL}
\alias{TPL}
\title{
Standardize a list of plant names according to The Plant List.
}
\description{
Connects to The Plant List (TPL) website and validates the names of a vector of
plant taxa, replacing synonyms by accepted names and removing orthographical
errors in plant names.
}
\usage{
TPL(splist, genus = NULL, species = NULL, infrasp = NULL, 
infra = TRUE, corr = TRUE, diffchar = 2, max.distance = 1, 
version = "1.1", encoding = "UTF-8", author = TRUE, 
drop.lower.level = FALSE, file = "", silent = TRUE, repeats = 6)
}
\arguments{
\item{splist}{A character vector specifying the input taxa, each element
including genus and specific epithet and, potentially, author name and
infraspecific abbreviation and epithet.}

\item{genus}{A character vector containing the genera of plant taxon names.}

\item{species}{A character vector containing the specific epithets of plant
taxon names.}

\item{infrasp}{A character vector containing the infraspecific epithets of
plant taxon names.}

\item{infra}{Logical. If TRUE (default), infraspecific epithets are used to
match taxon names in TPL.}

\item{corr}{Logical. If TRUE (default), spelling errors are corrected (only)
in the specific and infraspecific epithets prior to taxonomic standardization.}

\item{diffchar}{A number indicating the maximum difference between the number
of characters in corrected and original taxon names. Not used if corr = FALSE.}

\item{max.distance}{A number indicating the maximum distance allowed for a
match in \code{\link[base]{agrep}} when performing corrections of
spelling errors in specific epithets. Not used if corr = FALSE.}

\item{version}{A character vector indicating whether to connect to the newest
version of TPL (1.1) or to the older one (1.0). Defaults to "1.1".}

\item{encoding}{Encoding to be assumed for input strings from TPL website;
defaults to "UTF-8" (see \code{\link[utils]{read.csv}} and
\code{\link[base]{Encoding}}).}

\item{author}{Logical. If TRUE (default), the function tries to extract
author names from the input taxon (see \strong{Details}).}

\item{drop.lower.level}{Logical. If TRUE, the variety is dropped from the
input taxon if both subspecies and variety are given, and the forma is dropped
from the input taxon if both subspecies or variety and forma are given. If
specific and subspecific epithet are identical, the subspecies [variety] part
is dropped instead, and variety [forma] is kept. Defaults to FALSE.}

\item{file}{Either a character string naming a file or a connection open for
writing. "" (default) indicates output to the console.}

\item{silent}{Logical. If FALSE, the function prints the taxon name that is
currently processed in the output. Defaults to TRUE.}

\item{repeats}{A number indicating how many times \code{TPLck} should be
called if no connection to TPL website can be established (temporarily).}
}
\details{
The procedure used for taxonomic standardization is based on function
\code{\link{TPLck}}. A progress bar indicates the proportion of taxon names processed
so far. In case the TPL website cannot be reached temporarily, the function returns
an error but repeats trying to match the given taxon multiple times (see
\code{repeats}). If standardization is still not successful, the input taxon is
returned in field 'Taxon' with NA in all other fields.
}
\value{
A \code{data.frame} with the following components:

\item{$Taxon}{Original taxon name as provided in input.}

\item{$Genus}{Original genus name as provided in input.}

\item{$Hybrid.marker}{Hybrid marker, if taxon is indicated as hybrid in the input.}

\item{$Species}{Original specific epithet as provided in input.}

\item{$Abbrev}{Original abbreviation other than infraspecific rank included in
input taxon, including "cf.", "aff.", "agg.", "nom. cons.", "nom. cons. prop.",
"nom. inval.", "s.l.", and "s.str." and their orthographic variants.}

\item{$Infraspecific.rank}{Original infraspecific rank abbreviation as provided
in input, including "subsp.", "var.", "f.", and their orthographic variants.}

\item{$Infraspecific}{Original infraspecific epithet as provided in
input. If infra = FALSE, this is not shown.}

\item{$Authority}{Original author of taxon name as provided in input.}

\item{$ID}{The Plant List record ID of the matched taxon before resolving synonyms.}

\item{$Plant.Name.Index}{Logical. If TRUE the name is in TPL. If a taxon at
infraspecific level is not in TPL, Plant.Name.Index equals FALSE, except for
nominal infraspecies. Also compare \code{Higher.level}.}

\item{$TPL.version}{Version of TPL used.}

\item{$Taxonomic.status}{Taxonomic status of the matched taxon in TPL, either
'Accepted', 'Synonym', 'Unresolved', or 'Misapplied'.}

\item{$Family}{Family name, extracted from TPL for the valid form of the taxon.}

\item{$New.Genus}{Genus name, extracted from TPL for the valid form of the taxon.}

\item{$New.Hybrid.marker}{Hybrid marker, extracted from TPL for the valid form
of the taxon.}

\item{$New.Species}{Specific epithet, extracted from TPL for the valid form of the taxon.}

\item{$New.Infraspecific.rank}{Infraspecific rank abbreviation, extracted from TPL for
the valid form of the taxon, including "subsp.", "var." and "f.".}

\item{$New.Infraspecific}{Infraspecific epithet, extracted from TPL for the
valid form of the taxon.}

\item{$New.Authority}{Author of taxon name, extracted from TPL for the valid
form of the taxon.}

\item{$New.ID}{The Plant List record ID of the taxon, once synonyms have been
replaced by valid names. For accepted and unresolved names, this field will
be equivalent to \code{ID}.}

\item{$New.Taxonomic.status}{Taxonomic status of the resolved taxon in TPL,
once synonyms have been replaced by valid names. 'Accepted' or 'Unresolved'.}

\item{$Typo}{Logical. If TRUE there was a spelling error in the specific or
infraspecific epithet that has been corrected.}

\item{$WFormat}{Logical. If TRUE, fields in TPL had the wrong format for
information to be automatically extracted as they were not properly tabulated
or, alternatively, there was not a unique solution.}

\item{$Higher.level}{Logical. If TRUE, the input taxon is at infraspecific
level and does not occur in TPL, and the higher (species) level is provided in
the output instead. Also see \code{Plant.Name.Index}.}

\item{$Date}{Current date according to \code{\link[base]{Sys.Date}}.}
}
\references{
Cayuela, L., Granzow-de la Cerda, I., Albuquerque, F.S. and Golicher, J.D. 2012.
Taxonstand: An R package for species names standardization in vegetation databases.
\emph{Methods in Ecology and Evolution}, 3(6): 1078-1083.

Kalwij, J.M. 2012. Review of 'The Plant List, a working list of all plant species'.
\emph{Journal of Vegetation Science}, 23(5): 998-1002.
}
\author{
  Luis Cayuela & Anke Stein
}
\seealso{
\code{\link{TPLck}}.
}
\examples{
\dontrun{
data(bryophytes)

# Species names in full
r1 <- TPL(bryophytes$Full.name[1:20], corr = TRUE)
str(r1)

# A separate specification for genera, specific, and infraspecific
# epithets
r2 <- TPL(genus = bryophytes$Genus, species = bryophytes$Species,
infrasp = bryophytes$Intraspecific, corr = TRUE)
str(r2)

#----------------------------------
# An example using data from GBIF
#----------------------------------
# Download all records available in GBIF of all species within genus
# Liriodendron (GBIF table; note that a list of species can be also
# downloaded from GBIF for a defined geographical area)
require(dismo)
liriodendron <- gbif("Liriodendron", "*", geo = TRUE)

# Perform taxonomic standardization on plant names (TPL table)
sp.check <- TPL(unique(liriodendron$scientificName), infra = TRUE, 
corr = TRUE)
head(sp.check)

# Join GBIF table with TPL table
require(dplyr)
liriodendron.check <- liriodendron \%>\%
left_join(., sp.check, by = c("scientificName" = "Taxon"))

}
}
\keyword{ vegetation analysis }
