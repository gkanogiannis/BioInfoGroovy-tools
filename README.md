[![DOI](https://zenodo.org/badge/428339664.svg)](https://zenodo.org/badge/latestdoi/428339664)

Bio-informatics tools written in the groovy programming/scripting language.

ExternalSort.groovy : Sorts reads of a fastq or fasta file according to their ID.

SynchronizePairReads.groovy : Synchronizes paired reads.

*In order to pass options to the JVM before executing a script, 
set them through the JAVA_OPTS environment variable.
For example, to increase the available RAM for the JVM to 4 GBytes, 
use 'export JAVA_OPTS="$JAVA_OPTS -Xms4G"' before running the groovy script.
