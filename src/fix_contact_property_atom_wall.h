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

    Arno Mayrhofer (CFDEMresearch GmbH, Linz)
    Christoph Kloss (DCS Computing GmbH, Linz)

    Copyright 2016-     CFDEMresearch GmbH, Linz
    Copyright 2014-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(contactproperty/atom/wall,FixContactPropertyAtomWall)

#else

#ifndef LMP_FIX_CONTACT_PROPERTY_ATOM_WALL_H
#define LMP_FIX_CONTACT_PROPERTY_ATOM_WALL_H

#include "fix_contact_property_atom.h"
#include "fix_mesh_surface.h"

namespace LAMMPS_NS {

class FixContactPropertyAtomWall : public FixContactPropertyAtom {
    friend class Neighbor;
    friend class PairGran;

public:
    FixContactPropertyAtomWall(class LAMMPS *, int, char **);
    ~FixContactPropertyAtomWall();

    void clear();
    bool haveContact(const int iP, const int idTri, double *&history) const;

    FixMeshSurface *getMesh() const;

private:

    FixMeshSurface *fix_mesh_surface_;
    class FixPropertyAtom *fix_nneighs_;
    class PrimitiveWall *primitive_wall_;
};

}

#endif
#endif
