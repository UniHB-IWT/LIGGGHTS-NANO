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
    This file is from LAMMPS, but has been modified. Copyright for
    modification:

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz

    Copyright of original file:
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS

CommandStyle(read_data,ReadData)

#else

#ifndef LMP_READ_DATA_H
#define LMP_READ_DATA_H

#include <stdio.h>
#include "pointers.h"

namespace LAMMPS_NS {

class ReadData : protected Pointers {
 public:
  ReadData(class LAMMPS *);
  ~ReadData();
  void command(int, char **);

 private:
  int me;
  char *line,*keyword,*buffer;
  FILE *fp;
  int narg,maxarg,compressed;
  char **arg;

  int nfix;           // # of extra fixes that process/store info in data file
  int *fix_index;
  char **fix_header;
  char **fix_section;

  bigint nellipsoids;
  class AtomVecEllipsoid *avec_ellipsoid;
  bigint nlines;
  class AtomVecLine *avec_line;

  int add_to_existing;  
  bigint natoms_add;    

  void open(char *);
  void scan(int &, int &, int &, int &);
  int reallocate(int **, int, int);
  void header(int, int add = 1); 
  void parse_keyword(int, int);
  void skip_lines(int);
  void parse_coeffs(char *, const char *, int);

  void atoms();
  void velocities();
  void bonus(bigint, class AtomVec *, const char *);

  void bonds();
  void angles();
  void dihedrals();
  void impropers();

  void mass();
  void paircoeffs();
  void pairIJcoeffs();
  void bondcoeffs();
  void anglecoeffs(int);
  void dihedralcoeffs(int);
  void impropercoeffs(int);

  void fix(int, char *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot read_data after simulation box is defined

The read_data command cannot be used after a read_data,
read_restart, or create_box command.

E: Cannot run 2d simulation with nonperiodic Z dimension

Use the boundary command to make the z dimension periodic in order to
run a 2d simulation.

E: Fix ID for read_data does not exist

Self-explanatory.

E: Must read Atoms before Velocities

The Atoms section of a data file must come before a Velocities
section.

E: Invalid data file section: Ellipsoids

Atom style does not allow ellipsoids.

E: Must read Atoms before Ellipsoids

The Atoms section of a data file must come before a Ellipsoids
section.

E: Invalid data file section: Lines

Atom style does not allow lines.

E: Must read Atoms before Lines

The Atoms section of a data file must come before a Lines section.

E: Invalid data file section: Triangles

Atom style does not allow triangles.

E: Must read Atoms before Triangles

The Atoms section of a data file must come before a Triangles section.

E: Invalid data file section: Bodies

Atom style does not allow bodies.

E: Must read Atoms before Bodies

The Atoms section of a data file must come before a Bodies section.

E: Invalid data file section: Bonds

Atom style does not allow bonds.

E: Must read Atoms before Bonds

The Atoms section of a data file must come before a Bonds section.

E: Invalid data file section: Angles

Atom style does not allow angles.

E: Must read Atoms before Angles

The Atoms section of a data file must come before an Angles section.

E: Invalid data file section: Dihedrals

Atom style does not allow dihedrals.

E: Must read Atoms before Dihedrals

The Atoms section of a data file must come before a Dihedrals section.

E: Invalid data file section: Impropers

Atom style does not allow impropers.

E: Must read Atoms before Impropers

The Atoms section of a data file must come before an Impropers
section.

E: Must define pair_style before Pair Coeffs

Must use a pair_style command before reading a data file that defines
Pair Coeffs.

E: Must define pair_style before PairIJ Coeffs

UNDOCUMENTED

E: Invalid data file section: Bond Coeffs

Atom style does not allow bonds.

E: Must define bond_style before Bond Coeffs

Must use a bond_style command before reading a data file that
defines Bond Coeffs.

E: Invalid data file section: Angle Coeffs

Atom style does not allow angles.

E: Must define angle_style before Angle Coeffs

Must use an angle_style command before reading a data file that
defines Angle Coeffs.

E: Invalid data file section: Dihedral Coeffs

Atom style does not allow dihedrals.

E: Must define dihedral_style before Dihedral Coeffs

Must use a dihedral_style command before reading a data file that
defines Dihedral Coeffs.

E: Invalid data file section: Improper Coeffs

Atom style does not allow impropers.

E: Must define improper_style before Improper Coeffs

Must use an improper_style command before reading a data file that
defines Improper Coeffs.

E: Invalid data file section: BondBond Coeffs

Atom style does not allow angles.

E: Must define angle_style before BondBond Coeffs

Must use an angle_style command before reading a data file that
defines Angle Coeffs.

E: Invalid data file section: BondAngle Coeffs

Atom style does not allow angles.

E: Must define angle_style before BondAngle Coeffs

Must use an angle_style command before reading a data file that
defines Angle Coeffs.

E: Invalid data file section: MiddleBondTorsion Coeffs

Atom style does not allow dihedrals.

E: Must define dihedral_style before MiddleBondTorsion Coeffs

Must use a dihedral_style command before reading a data file that
defines MiddleBondTorsion Coeffs.

E: Invalid data file section: EndBondTorsion Coeffs

Atom style does not allow dihedrals.

E: Must define dihedral_style before EndBondTorsion Coeffs

Must use a dihedral_style command before reading a data file that
defines EndBondTorsion Coeffs.

E: Invalid data file section: AngleTorsion Coeffs

Atom style does not allow dihedrals.

E: Must define dihedral_style before AngleTorsion Coeffs

Must use a dihedral_style command before reading a data file that
defines AngleTorsion Coeffs.

E: Invalid data file section: AngleAngleTorsion Coeffs

Atom style does not allow dihedrals.

E: Must define dihedral_style before AngleAngleTorsion Coeffs

Must use a dihedral_style command before reading a data file that
defines AngleAngleTorsion Coeffs.

E: Invalid data file section: BondBond13 Coeffs

Atom style does not allow dihedrals.

E: Must define dihedral_style before BondBond13 Coeffs

Must use a dihedral_style command before reading a data file that
defines BondBond13 Coeffs.

E: Invalid data file section: AngleAngle Coeffs

Atom style does not allow impropers.

E: Must define improper_style before AngleAngle Coeffs

Must use an improper_style command before reading a data file that
defines AngleAngle Coeffs.

E: Unknown identifier in data file: %s

A section of the data file cannot be read by LAMMPS.

E: No atoms in data file

The header of the data file indicated that atoms would be included,
but they were not present.

E: Unexpected end of data file

LAMMPS hit the end of the data file while attempting to read a
section.  Something is wrong with the format of the data file.

E: No ellipsoids allowed with this atom style

Self-explanatory.  Check data file.

E: No lines allowed with this atom style

Self-explanatory.  Check data file.

E: No triangles allowed with this atom style

Self-explanatory.  Check data file.

E: No bodies allowed with this atom style

Self-explanatory.  Check data file.

E: System in data file is too big

See the setting for bigint in the src/lmptype.h file.

E: No bonds allowed with this atom style

Self-explanatory.  Check data file.

E: No angles allowed with this atom style

Self-explanatory.  Check data file.

E: No dihedrals allowed with this atom style

Self-explanatory.  Check data file.

E: No impropers allowed with this atom style

Self-explanatory.  Check data file.

E: Bonds defined but no bond types

The data file header lists bonds but no bond types.

E: Angles defined but no angle types

The data file header lists angles but no angle types.

E: Dihedrals defined but no dihedral types

The data file header lists dihedrals but no dihedral types.

E: Impropers defined but no improper types

The data file header lists improper but no improper types.

E: Did not assign all atoms correctly

Atoms read in from a data file were not assigned correctly to
processors.  This is likely due to some atom coordinates being
outside a non-periodic simulation box.

E: Invalid atom ID in Atoms section of data file

Atom IDs must be positive integers.

E: Too many lines in one body in data file - boost MAXBODY

MAXBODY is a setting at the top of the src/read_data.cpp file.
Set it larger and re-compile the code.

E: Bonds assigned incorrectly

Bonds read in from the data file were not assigned correctly to atoms.
This means there is something invalid about the topology definitions.

E: Angles assigned incorrectly

Angles read in from the data file were not assigned correctly to
atoms.  This means there is something invalid about the topology
definitions.

E: Dihedrals assigned incorrectly

Dihedrals read in from the data file were not assigned correctly to
atoms.  This means there is something invalid about the topology
definitions.

E: Impropers assigned incorrectly

Impropers read in from the data file were not assigned correctly to
atoms.  This means there is something invalid about the topology
definitions.

E: Molecular data file has too many atoms

These kids of data files are currently limited to a number
of atoms that fits in a 32-bit integer.

E: Needed topology not in data file

The header of the data file indicated that bonds or angles or
dihedrals or impropers would be included, but they were not present.

E: Needed bonus data not in data file

Some atom styles require bonus data.  See the read_data doc page for
details.

E: Cannot open gzipped file

LAMMPS is attempting to open a gzipped version of the specified file
but was unsuccessful.  Check that the path and name are correct.

E: Cannot open file %s

The specified file cannot be opened.  Check that the path and name are
correct.

*/
