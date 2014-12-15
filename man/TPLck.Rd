\name{TPLck}
\alias{TPLck}
\title{
Connects to The Plant List (TPL) website and check for the validity of a plant name
}
\description{
Connects to TPL and validates the name of a single plant species name, replacing synonyms for accepted names and removing orthographical errors in plant names
}
\usage{
TPLck(sp, corr = TRUE, diffchar = 2, max.distance = 1, infra = TRUE, abbrev = TRUE, 
version = "1.1", encoding = "UTF-8")
}
\arguments{
  \item{sp}{
A character specifying the genus and specific epithet and, additionally, the infraspecific epithet.
}
  \item{corr}{
Logical. If 'TRUE' (default), then removal of spelling errors is performed (only) on the specific and infraspecific epithets prior to taxonomic standardization.
}
  \item{diffchar}{
A number indicating the maximum difference between number of characters in corrected and original species names. Not used if corr=FALSE.
}
  \item{max.distance}{
Maximum distance allowed for a match in \code{agrep} function when performing corrections of spelling errors in specific epithets. Not used if corr=FALSE.
}
  \item{infra}{
Logical. If 'TRUE' (default) then infraspecific epithets are used to validate the taxonomic status of species names in TPL.
}
  \item{abbrev}{
Logical. If 'TRUE' (default), abbreviations (aff., cf., subsp., var.) and their variants are removed prior to taxonomic standardization.
}
  \item{version}{
A character indicating whether to connect to the newest version of TPL (1.1) or to the older one (1.0). Defaults to "1.1".
}
  \item{encoding}{
Encoding to be assumed for input strings.  It is used to mark character strings as known to be in Latin-1 or UTF-8 (see 'Encoding'). It is not used to re-encode the input, but allows R to handle encoded strings in their native encoding (if one of those two). Defaults to "UTF-8" but it can be changed to "unknown" or "Latin-1".
}
}
\details{
The procedure searches for a species name in The Plant List (TPL) and provides its taxonomic status. If the status is either 'Accepted' or 'Unresolved' (i.e. names for which the contributing data sources did not contain sufficient evidence to decide whether they were 'Accepted' or 'Synonyms'), the function returns the species name unchanged. In those cases where the species name provided as input is recognised as a "Synonym", it is replaced with the current accepted name. Some data sets which contributed to The Plant List record not only how plant names should be used but also where in the published literature a given name may previously have been used inappropriately (to refer erroneously to another species). In those cases, the function returns the accepted name of the species to which this name has been previously and erroneously applied. Orthographic errors can be corrected only in specific epithets. By increasing arguments 'diffchar' and 'max.distance', larger differences can be detected in typos, but this also increases false positives (i.e. replacement of some names for others that do not really match), so some caution is recommended here. If 'infra=FALSE', then infraspecific epithets are neither considered for species name validation in TPL, nor returned in the output.
A new version of TPL was released in September 2013. Version 1.1 replaces Version 1.0, which still remains accessible. Version 1.1 includes new data sets, updated versions of the original data sets and improved algorithms to resolve logical conflicts between those data sets.
}

\value{
The function return an object of class \code{data.frame} with the following components:

\item{$Genus }{Original genus of species provided as input for taxonomic standardization.}
\item{$Species }{Original specific epithet of species provided as input for taxonomic standardization.}
\item{$Abbrev }{Standard annotation used in species epithet, including "cf.", "aff.", "s.l.", and "s.str." and their orthographic variants.}
\item{$Infraspecific }{Original intraspecific epithet of species provided as input for taxonomic standardization. If 'infra=FALSE', this is not shown. }
\item{$ID }{ The Plant List record 'id'.}
\item{$Plant.Name.Index }{Logical. If 'TRUE' the name is in TPL.}
\item{$Taxonomic.status}{Taxonomic status as in TPL, either 'Accepted', 'Synonym', 'Unresolved', or 'Misapplied'.}
\item{$TPL_version}{Version of TPL used.}
\item{$Family }{Family name, extracted from TPL for the valid form of the name.}
\item{$New.Genus }{Genus, extracted from TPL for the valid form of the name.}
\item{$New.Hybrid.marker }{Hybrid marker, extracted from TPL for the valid form of the name.}
\item{$New.Species }{Specific epithet, extracted from TPL for the valid form of the name.}
\item{$New.Infraspecific }{Infraspecific epithet, extracted from TPL for the valid form of the name.}
\item{$Authority }{A field designating the scientist(s) who first published the name, extracted from TPL for the valid form of the name.}
\item{$New.ID }{The Plant List record 'id' for the species name, once synonyms have been replaced by accepted names. In all other cases, this field will be equivalent to 'ID'.}
\item{$Typo }{Logical. If 'TRUE' there was a spelling error in the specific epithet that has been corrected.}
\item{$WFormat }{Logical. If 'TRUE', fields in TPL had the wrong format for information to be automatically extracted as they were not properly tabulated or, alternatively, there was not a unique solutions (see 'note').}
}
\references{
Cayuela, L., Granzow-de la Cerda, I., Albuquerque, F.S. and Golicher, J.D. 2012. Taxonstand: An R package for species names standardisation in vegetation databases. \emph{Methods in Ecology and Evolution}, 3(6): 1078-1083.

Kalwij, J.M. 2012. Review of 'The Plant List, a working list of all plant species'. \emph{Journal of Vegetation Science}, 23(5): 998-1002.

}
\author{
Luis Cayuela & Jari Oksanen
}
\note{
Various limitations have been identified to date in the implementation of this function. First, homonyms (i.e. a name for a taxon that is identical in spelling to another such name, but belongs to a different taxon) cannot be identified. Second, if the input infraspecific epithet does not match any of the infraspecific epithets provided in TPL, then the first accepted name with no infraspecific epithet is selected as the best match. If all names provided by TPL have infraspecific epithets and are all synonyms (e.g. Pottia starkeana), then there is no best match and the output will match the original name.
}
\seealso{
See also \code{\link{TPL}}.
}
\examples{
  \dontrun{
# An accepted name
sp1 <- TPLck("Amblystegium serpens juratzkanum")
sp1
# An unresolved name
sp2 <- TPLck("Bryum capillare cenomanicum")
sp2
# A synonym
sp3 <- TPLck("Pottia acutidentata")
sp3
# A misapplied name
sp4 <- TPLck("Colutea istria")
sp4
# A spelling error in specific epithet
sp5 <- TPLck("Pohlia longicola", corr=TRUE)
sp5
# A spelling error that is not corrected ('max.distance' defaults to 1)
sp6 <- TPLck("Microbryum curvicollum", corr=TRUE)
sp6
# If increasing 'max.distance', the spelling error is accounted for
sp7 <- TPLck("Microbryum curvicollum", corr=TRUE, max.distance=3)
sp7
}
}
\keyword{ vegetation analysis }
