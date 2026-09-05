\name{alabaster.vcf News}
\title{alabaster.vcf News}
\encoding{UTF-8}

\section{Version 1.14.0}{\itemize{
\item Officially deprecate \code{stageObject()} and \code{loadObject()} in favor of \code{saveObject()} and \code{readObject()}, respectively.
The former have been soft-deprecated and scrubbed from documentation for quite some time, so it is unlikely that anyone is using them.
We also deprecate the separate load/stage methods for \code{VCFHeader} objects.
}}
