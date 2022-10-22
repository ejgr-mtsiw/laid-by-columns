# laid-c-hdf5-mpi - 4th prototype
OpenMPI/C implementation of LAID with dataset stored in hdf5 file

In this mode we don't write the disjoint matrix (DM).
Everytime we need a line or column from the DM it's generated from the dataset in memory.

Each node root will open the original dataset file and read the contents to memory. This allows us to save memory in
each node, without sacrificing much performance, because we're not having to send data across nodes.

The node root(s) sort the dataset in memory, remove duplicates and add jnsqs bits if necessary.
The dataset is then sorted by class, and shared with all the processes in the same node.

The dataset attributes are split among all processes, and each process calculates the totals for their attributes.
The best total from each attribute is shared using MPI and the best overall attribute is selected.

Each node root builds and shares the data column associated with the best attribute, that is used by all processes to
update the covered lines list and calculate the updated totals. The fact that all processes use the same covering data
column, even if they don't have the best attribute in its list, guarantees that the solution found applies to the
entire original dataset. But can be differente from the 1-process solution or even change between runs, depending on
the way the MPI_Allreduce selects attributes with identical totals.


TLDR:
Each node root
 - Reads dataset attributes from hdf5 file
 - Read dataset
 - Sort dataset
 - Remove duplicates
 - Add jnsqs
 - Sorts dataset by class

All processes
 - Divide dataset attributes between available processes
 - Apply set covering algorithm

Global root
 - Show solution

## Set covering algorithm
ALL:
 - Zero line covered array
 
ROOT:
 - Zero attribute covered array

ALL:
 - Calculate attributes totals

LOOP:
ALL:
 - Calculate best local attribute
 - MPI_Allreduce to select best global attribute
 
 ROOT:
 - Mark best attribute as selected in attribute covered array
 
 ALL:
 - If there are no more lines to blacklist:
  - goto SHOW_SOLUTION

NODE ROOTS:
 - Generate cover line for the best atribute
 - MPI_Bcast cover line to everyone else on the same node

PROCESSES WITHOUT THE BEST ATTRIBUTE:
 - Wait for cover line

ALL:
 - Update atributes totals based on the received cover line and current covered line array
 - Update covered line array

ALL:
 - goto LOOP

SHOW_SOLUTION:
ROOT:
 - Display list of selected attributes
 