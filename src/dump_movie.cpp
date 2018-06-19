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

/* ----------------------------------------------------------------------
   Contributing author:  Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dump_movie.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DumpMovie::DumpMovie(LAMMPS *lmp, int narg, char **arg) :
  DumpImage(lmp, narg, arg)
{
  if (multiproc || compressed || multifile)
    error->all(FLERR,"Invalid dump movie filename");

  filetype = PPM;
  bitrate = 2000;
  framerate = 24;
  fp = NULL;
}

/* ---------------------------------------------------------------------- */

void DumpMovie::openfile()
{
  char moviecmd[1024];

  if ((comm->me == 0) && (fp == NULL)) {

#ifdef LAMMPS_FFMPEG
    sprintf(moviecmd,"ffmpeg -v error -y -r %.2f -f image2pipe -c:v ppm -i - "
            "-r 24.0 -b:v %dk %s ", framerate, bitrate, filename);
#else
    error->one(FLERR,"Cannot generate movie file");
#endif

#if defined(_WIN32)
    fp = _popen(moviecmd,"wb");
#else
    fp = popen(moviecmd,"w");
#endif

    if (fp == NULL) {
      char str[512];
      sprintf(str,"Failed to open FFmpeg pipeline to file %s",filename);
      error->one(FLERR,str);
    }
  }
}
/* ---------------------------------------------------------------------- */

void DumpMovie::init_style()
{
  // initialize image style circumventing multifile check

  multifile = 1;
  DumpImage::init_style();
  multifile = 0;
}

/* ---------------------------------------------------------------------- */

int DumpMovie::modify_param(int narg, char **arg)
{
  int n = DumpImage::modify_param(narg,arg);
  if (n) return n;

  if (strcmp(arg[0],"bitrate") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    bitrate = force->inumeric(FLERR,arg[1]);
    if (bitrate <= 0.0) error->all(FLERR,"Illegal dump_modify command");
    return 2;
  }

  if (strcmp(arg[0],"framerate") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    framerate = force->numeric(FLERR,arg[1]);
    if ((framerate <= 0.1) || (framerate > 24.0))
      error->all(FLERR,"Illegal dump_modify framerate command");
    return 2;
  }

  return 0;
}

