# laid-by-columns - 3rd prototype
OpenMPI/C implementation of LAID with dataset stored in hdf5 file

# General description
In this mode we don't write the disjoint matrix (DM) to file.
Everytime we need a line or column from the DM it's generated from the dataset in memory.

Each node root will open the original dataset file and read the contents to memory. This allows us to save memory in
each node, without sacrificing much performance, because we're not having to send data across nodes.

The node root(s) sort the dataset in memory, remove duplicates and add jnsqs bits if necessary.
The dataset is then shared with all the processes in the same node.

The dataset attributes (columns) are split among all processes, and each process calculates the totals for their attributes.
The best total from each process is shared using MPI and the best overall attribute is selected.

Each node root builds and shares the data column associated with the best attribute, that is used by all processes to
update the covered lines list and calculate the updated totals. The fact that all processes use the same covering data
column, even if they don't have the best attribute in its list, guarantees that the solution found applies to the
entire original dataset. But can be differente from the 1-process solution or even change between runs, depending on
the way the MPI_Allreduce selects the best attribute when several attributes have identical totals.

# TLDR:

* ALL:
  * Reads dataset attributes from hdf5 file
* NODE ROOT
  * Read dataset
  * Sort dataset
  * Remove duplicates
  * Add jnsqs
*ALL
  * Apply set covering algorithm

* ROOT
  * Show solution

## Set covering algorithm
ALL:
 - Calculate attributes totals

LOOP:
While there are lines to blacklist in the disjoint matrix:
ALL:
 - Calculate best local attribute
 - MPI_Allreduce to select best global attribute
 
 ROOT:
 - Mark best attribute as selected
 
 ALL:
 - Generate column for the best atribute
 - Update atributes totals based on the received best attribute column line and current covered line bit array
 - Update covered line array

ALL:
 - goto LOOP

SHOW_SOLUTION:
ROOT:
 - Display list of selected attributes
 