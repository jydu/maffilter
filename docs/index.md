MafFilter is a program dedicated to the analysis of genome alignments. It parses and manipulates [MAF files](https://genome.ucsc.edu/FAQ/FAQformat.html#format5) as well as more simple fasta files. Despite various filtering options and format conversion tools, MafFilter can compute a wide range of statistics (phylogenetic trees, nucleotide diversity, inferrence of selection, etc.). Current version is 1.2.1.


## What can MafFilter do?

!MafFilter applies a series of "filters" to a MAF file, in order to clean it, extract data and computer statistics while keeping track of the associated meta-data such as genome coordinates and quality scores.

* It can process the alignment to remove low-quality / ambiguous / masked regions.
* It can export data into a single or multiple alignment file in format such as Fasta or Clustal.
* It can read annotation data in GFF or GTF format, and extract the corresponding alignment.
* It can perform sliding windows calculations.
* It can reconstruct phylogeny/genealogy along the genome alignment.
* It can compute population genetics statistics, such as site frequency spectrum, number of fixed/polymorphic sites, etc.

## How can I get it?

MafFilter is built using the [Bio++ libraries](http://biopp.univ-montp2.fr), as well as the boost iostream library for handling of compressed files. [Debian](https://packages.debian.org/search?keywords=maffilter&searchon=names&suite=stable&section=all) and [RPM](https://download.opensuse.org/repositories/home:/jdutheil:/Bio++2.3.0/) packages are available.

For compiling the programs yourself, from the downloaded sources or from the git repository, please follow the instructions from the [Bio++ website](http://biopp.univ-montp2.fr/wiki/index.php/Installation).

## How do I use it? 

Several example data sets are distributed along with the source code of the package. A reference manual is also available [here](http://biopp.univ-montp2.fr/manual/html/maffilter/), or can be downloaded as [PDf](http://biopp.univ-montp2.fr/manual/pdf/maffilter/). Questions can be asked on the dedicated forum: [here](https://groups.google.com/forum/?hl=en#!forum/maffilter).
