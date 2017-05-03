\name{TPLck}
\alias{TPLck}
\title{
Standardize a plant name according to The Plant List.
}
\description{
Connects to The Plant List (TPL) website and validates the name of a single
plant taxon, replacing synonyms by accepted names and removing orthographical
errors in plant names.
}
\usage{
TPLck(sp, infra = TRUE, corr = TRUE, diffchar = 2, 
max.distance = 1, version = "1.1", encoding = "UTF-8", 
author = TRUE, drop.lower.level = FALSE)
}
\arguments{
\item{sp}{A character vector specifying the input taxon, i.e. genus and
specific epithet and, potentially, author name and infraspecific abbreviation
and epithet.}

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
}
\details{
The function searches for a taxon name on The Plant List (TPL) website and provides
its taxonomic status (\url{http://www.theplantlist.org}). If the status is either 
'Accepted' or 'Unresolved' (i.e. names for which the contributing 
data sources did not contain sufficient evidence to decide whether
they were 'Accepted' or 'Synonyms'), the function returns the taxon name
unchanged. In cases where the input taxon is recognised as a 'Synonym', the
according valid taxon is provided in the output, i.e. the current accepted
name or, in some cases, an unresolved name. Some data sets which contributed to
TPL record not only how plant names should be used but also where in the
published literature a given name may previously have been used inappropriately
(to refer erroneously to another species). In those cases, the taxonomic status
is 'Misapplied', and the function returns the name of the taxon to which this
name has been previously and erroneously applied.

If the author of the taxon name is provided in the input, it will be used for
taxon matching and distinction between homonyms. The best author match is based
on \code{\link[utils]{adist}}, i.e. spelling variations are taken into account.
Both the full authority and the current author only (i.e. dropping potential
basionym author or author preceding the word 'ex') are compared
(see \strong{Examples}). Warnings are given in case of imperfect author match
(except if differences are in spaces, dots, or '&' vs. 'et' only). The function
distinguishes between author name and infraspecific epithet based on
uppercase/lowercase only, so this procedure does not work properly if an
infraspecific epithet is capitalized or if an author name is spelled with
lowercase (but several lowercase author components like "auct", ex", "f./fil."
etc. are accounted for).
Note that if an input taxon at infraspecific level does not exist in TPL and
the species level is given in the output (with Plant.Name.Index = FALSE and
Higher.level = TRUE), the function so far does not use the species-level author
(if given in the input) for picking the corresponding species.

If TPL includes multiple entries for a given input taxon, the function tries to
match the correct taxon based on the author, infraspecific rank and hybrid
marker, if provided. Otherwise, preference is given to an Accepted name over a
Synonym, Misapplied, or Unresolved name (in this order). If the input taxon is
at the infraspecific level and does not exist on the TPL website with the given
infraspecific rank abbreviation, the taxon is matched to a different
infraspecific rank, if possible (e.g. 'subsp.' versus 'var.'). In case of multiple
synonyms, taxon names with a higher confidence level in TPL (regarding the
status of name records; see website) are preferred. Otherwise, the first entry
on the website that is not an Illegitimate or Invalid name (if possible) is
selected, and a warning is given.

Orthographic errors can be corrected (only) in specific and infraspecific
epithets. By increasing arguments 'diffchar' and 'max.distance', larger
differences can be detected in typos, but this also increases false positives
(i.e. replacement of some names for others that do not really match), so some
caution is recommended here. Note that if the correction would result in
multiple equally likely names (e.g. \emph{Acacia macrocantha} could be
corrected to \emph{A. macracantha} or \emph{A. microcantha}), no correction
will take place, and Plant.Name.Index is set to FALSE.
If a specific epithet cannot be matched after orthographic correction based on
\code{\link[base]{agrep}}, common endings are replaced by the equivalent
masculine/feminine/neuter or genitive form  (-a/-um/-us, -is/-e, -ii/-i) and
the taxon is tried to match again.

If 'infra = FALSE', then infraspecific epithets are neither considered for
species name validation in TPL, nor returned in the output.

The latest version of TPL, version 1.1, was released in September 2013.
Version 1.1 replaces version 1.0, which still remains accessible. Version 1.1
includes new data sets, updated versions of the original data sets and improved
algorithms to resolve logical conflicts between those data sets.

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
Taxonstand: An R package for species names standardisation in vegetation databases.
\emph{Methods in Ecology and Evolution}, 3(6): 1078-1083.

Kalwij, J.M. 2012. Review of 'The Plant List, a working list of all plant species'.
\emph{Journal of Vegetation Science}, 23(5): 998-1002.
}
\author{
Luis Cayuela, Anke Stein & Jari Oksanen
}
\seealso{
\code{\link{TPL}}.
}
\examples{
\dontrun{
# An accepted name
(TPLck("Amblystegium serpens juratzkanum"))

# An unresolved name
(TPLck("Bryum capillare cenomanicum"))

# A synonym
(TPLck("Pottia acutidentata"))

# A misapplied name
(TPLck("Colutea istria"))

# A name that is not in TPL
(TPLck("Hypochoeris balbisii"))

# A spelling error in the specific epithet
(TPLck("Pohlia longicola", corr = TRUE))

# A spelling error that is not corrected ('max.distance' defaults to 1)
(TPLck("Microbryum curvicollum", corr = TRUE))

# If increasing 'max.distance', the spelling error is accounted for
(TPLck("Microbryum curvicollum", corr = TRUE, max.distance = 3))

# A spelling error where the ending is changed to the 
# neuter/feminine form
(TPLck("Symphytum officinalis"))
(TPLck("Schinus terebinthifolium"))

# A spelling error that is not corrected because two different results
# are possible (see Details)
(TPLck("Acacia macrocantha"))

# A taxon matched through author name
(TPLck("Gladiolus communis L. subsp. byzantinus (Mill.) A.P.Ham."))

# If only the current author is provided (without the author of the
# basionym), the function still matches the correct taxon (even though
# adist returns a better overall match for the author of the homonym,
# Abies alba Mill.)
(TPLck("Abies alba Michx."))

# A difference between TPL versions 1.0 and 1.1
(TPLck("Fallopia japonica", version = "1.0"))
(TPLck("Fallopia japonica", version = "1.1"))

# Avoid illegitimate names when choosing between multiple synonyms
(TPLck("Anthemis altissima"))

# A nominal subspecies not in TPL (Higher.level == TRUE; 
# Plant.Name.Index == TRUE)
(TPLck("Callitriche brutia Petagna subsp. brutia"))

# A variety not in TPL (Higher.level == TRUE; Plant.Name.Index == FALSE)
(TPLck("Asplenium ruta-muraria var. lanceolum"))

# A taxon matched through infraspecific rank abbreviation
(TPLck("Heliopsis helianthoides subsp. scabra"))

# Drop variety and keep subspecies
(TPLck("Vicia sativa L. subsp. nigra (L.) Ehrh. var. minor (Bertol.) Gaudin",
drop.lower.level = TRUE))

# Drop nominal subspecies and keep variety
(TPLck("Anagallis arvensis subsp. arvensis var. caerulea", 
drop.lower.level = TRUE))
}
}
\keyword{ vegetation analysis }
