# laid-c-hdf5-mpi
OpenMPI/C implementation of LAID with dataset stored in hdf5 file

In this mode we don't write the disjoint matrix (DM).
Everytime we need a line or column from the DM it's generated from the dataset in memory.

Each node root will open the original dataset file and read the contents to memory. This allows us to save memory in
each node, without sacrificing much performance, because we're not having to send data across nodes.

The node root(s) sort the dataset in memory, remove duplicates and  adds jnsqs bits if necessary.

Then they build a list of the steps needed to generate the disjoint matrix.
This list of steps allows us to generate any line or column of the disjoint matrix.

TLDR:
Each node root
 - Reads dataset attributes from hdf5 file
 - Read dataset
 - Sort dataset
 - Remove duplicates
 - Add jnsqs
 - Builds steps for matrix generation

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
 - Calculate best attribute
 - MPI_Allreduce to select best global attribute
 - If there are no more lines to blacklist (best attribute == -1):
  - goto SHOW_SOLUTION

ROOT:
 - Mark best attribute as selected in attribute covered array

NODE ROOTS:
 - Generate cover line for the best atribute
 - MPI_Bcast to everyone else on the same node

PROCESSES WITHOUT THE BEST ATTRIBUTE:
 - Wait for cover line

ALL:
 - Update atributes totals based on the received cover line
 - Update line covered array

ALL:
 - goto LOOP

SHOW_SOLUTION:
ROOT:
 - Display list of selected attributes
 
 