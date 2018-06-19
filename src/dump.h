/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    This file is from LAMMPS
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

#ifndef LMP_DUMP_H
#define LMP_DUMP_H

#include <mpi.h>
#include <stdio.h>
#include "pointers.h"
#include "sort_buffer.h"

namespace LAMMPS_NS {

class Dump : protected Pointers {

 friend class Info;

 public:
  char *id;                  // user-defined name of Dump
  char *style;               // style of Dump
  int igroup,groupbit;       // group that Dump is performed on

  int first_flag;            // 0 if no initial dump, 1 if yes initial dump
  int clearstep;             // 1 if dump invokes computes, 0 if not

  int comm_forward;          // size of forward communication (0 if none)
  int comm_reverse;          // size of reverse communication (0 if none)

  // static variable across all Dump objects

  static Dump *dumpptr;         // holds a ptr to Dump currently being used

  Dump(class LAMMPS *, int, char **);
  virtual ~Dump();
  void init();
  virtual void write();

  virtual int pack_comm(int, int *, double *, int, int *) {return 0;}
  virtual void unpack_comm(int, int, double *) {}
  virtual int pack_reverse_comm(int, int, double *) {return 0;}
  virtual void unpack_reverse_comm(int, int *, double *) {}

  void modify_params(int, char **);
  virtual bigint memory_usage();

 protected:
  int me,nprocs;             // proc info

  char *filename;            // user-specified file
  int compressed;            // 1 if dump file is written compressed, 0 no
  int binary;                // 1 if dump file is written binary, 0 no
  int multifile;             // 0 = one big file, 1 = one file per timestep

  int multiproc;             // 0 = proc 0 writes for all,
                             // else # of procs writing files
  int nclusterprocs;         // # of procs in my cluster that write to one file
  int filewriter;            // 1 if this proc writes a file, else 0
  int fileproc;              // ID of proc in my cluster who writes to file
  char *multiname;           // dump filename with % converted to cluster ID
  MPI_Comm clustercomm;      // MPI communicator within my cluster of procs

  int header_flag;           // 0 = item, 2 = xyz
  int flush_flag;            // 0 if no flush, 1 if flush every dump
  int append_flag;           // 1 if open file in append mode, 0 if not
  int buffer_allow;          // 1 if style allows for buffer_flag, 0 if not
  int buffer_flag;           // 1 if buffer output as one big string, 0 if not
  int padflag;               // timestep padding in filename
  int singlefile_opened;     // 1 = one big file, already opened, else 0

  char boundstr[9];          // encoding of boundary flags
  char *format_default;      // default format string
  char *format_user;         // format string set by user
  char *format;              // format string for the file write
  FILE *fp;                  // file to write dump to
  int size_one;              // # of quantities for one atom
  int nme;                   // # of atoms in this dump from me
  int nsme;                  // # of chars in string output from me

  double boxxlo,boxxhi;      // local copies of domain values
  double boxylo,boxyhi;      // lo/hi are bounding box for triclinic
  double boxzlo,boxzhi;
  double boxxy,boxxz,boxyz;

  bigint ntotal;             // total # of per-atom lines in snapshot

  int maxbuf;                // size of buf
  double *buf;               // memory for atom quantities

  int maxsbuf;               // size of sbuf
  char *sbuf;                // memory for atom quantities in string format

  SortBuffer *sortBuffer;   // class used for sorting buffers

  virtual void init_style() = 0;
  virtual void openfile();
  virtual int modify_param(int, char **) {return 0;}
  virtual void write_header(bigint) = 0;
  virtual int count();
  virtual void pack(int *) = 0;
  virtual int convert_string(int, double *) {return 0;}
  virtual void write_data(int, double *) = 0;
};

}

#endif

/* ERROR/WARNING messages:

E: Cannot dump sort when multiple procs write the dump file

UNDOCUMENTED

E: Cannot dump sort on atom IDs with no atom IDs defined

Self-explanatory.

E: Dump sort column is invalid

Self-explanatory.

E: Too many atoms to dump sort

Cannot sort when running with more than 2^31 atoms.

E: Too much per-proc info for dump

Number of local atoms times number of columns must fit in a 32-bit
integer for dump.

E: Cannot open gzipped file

LAMMPS was compiled without support for reading and writing gzipped
files through a pipeline to the gzip program with -DLAMMPS_GZIP.

E: Cannot open dump file

The specified file cannot be opened.  Check that the path and name are
correct. If the file is a compressed file, also check that the gzip
executable can be found and run.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot use dump_modify fileper without % in dump file name

UNDOCUMENTED

E: Cannot use dump_modify nfile without % in dump file name

UNDOCUMENTED

*/
