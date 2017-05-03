\name{NEWS}
\title{Taxonstand News}
\encoding{UTF-8}

\section{Changes in version 2.0}{

  \subsection{BUG FIXES}{
    \itemize{

      \item If no connection to TPL website can be established (temporarily),
      \code{TPL} now tries multiple times until connection is succesfully
      established.

      \item Fixed erroneous assignment of target taxa for synonyms. Previously,
      if there were multiple names listed on the website, \code{TPL}
      would pick the first entry instead of the designated target taxon.

      \item Fixed erroneous assignment if the genus of an input taxon has
      only a single species in TPL (e.g. under former version 
      \emph{Hypochoeris balbisii} was incorrectly resolved as 
      \emph{Hypochoeris argentina}).

      \item Fixed standardization for genus-level hybrids (previously, input 
      hybrid marker was recognized as genus name).

      \item Fixed issues with special characters (e.g. "×" causing problems
      under Windows).

    } %itemize
  } % bug fixes

  \subsection{NEW FEATURES}{
    \itemize{

      \item Use author name for matching taxon if provided in the input.

      \item Use infraspecific rank abbreviation and hybrid marker for 
      matching taxon if there are multiple otherwise identical synonyms.

      \item Avoid illegitimate and invalid names when picking one of
      multiple synonyms.

      \item Added some cleaning steps (i.e. replace ligatures and trema, 
      e.g. replace 'ﬃ' by 'ffi', and 'ë' by 'e'; add space before opening 
      parenthesis, if missing).

      \item Extended correction of spelling mistakes to account for potential
      different endings of specific epithet (e.g. -um vs. -a).

      \item Added option to drop one of the lower levels if multiple infraspecific
      levels are provided in the input.

      \item Return Plant.Name.Index == FALSE if taxon at infraspecific level 
      does not exist in \code{TPL} (except for nominal infraspecies).

      \item New field 'Taxon': complete input name.

      \item New field 'Hybrid.marker': hybrid marker of input name.

      \item New field 'Infraspecific.rank': rank abbreviation of input name; 
      field 'Abbrev' now only includes other abbreviations like "cf", "s.str."

      \item Previous field 'Authority' renamed to 'New.Authority'; new field
      'Authority': author as provided in input name.

      \item New field 'New.Infraspecific.rank': rank abbreviation of 
      valid/standardized name.

      \item New field 'New.Taxonomic.status': taxonomic status of 
      valid/standardized name (in most cases 'Accepted', but also 'Unresolved').

      \item New field 'Higher.level': indicates if infraspecific level does not
      exist in \code{TPL} and higher (species) level was returned.

      \item New field 'Date': current date.

      \item Show progress bar indicating the proportion of taxon names processed
      so far.

      \item Avoid Invalid names in addition to Illegitimate names when 
      choosing between multiple synonyms.

      \item Correct function for cases where only infraspecies occur 
      in TPL, e.g. \emph{Pilosella bauhinii}.

      \item Correct/remove some warnings.
, 
      \item Add "des" and "v." to authorstrings.

      \item Add/Change warning in case of multiple possible corrections
      of typo (e.g. Aloe x schoenlandi).

    } %itemize
  } % features

} % 2.0
