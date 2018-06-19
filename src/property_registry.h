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
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Richard Berger (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef PROPERTY_REGISTRY_H_
#define PROPERTY_REGISTRY_H_

#include "pointers.h"
#include "lammps.h"
#include "fix_property_global.h"
#include <map>
#include <set>
#include <string>
#include "properties.h"
#include "error.h"
#include "modify.h"

using namespace LAMMPS_NS;

class PropertyRegistry;

template<typename T>
class Property {
public:
  Property() : data(0) {}
  virtual ~Property(){}

  void connect(T & variable) {
    listeners.insert(&variable);
    variable = data;
  }

  void disconnect(T & variable) {
    typename std::set<T*>::iterator it = listeners.find(&variable);
    listeners.erase(it);
    variable = NULL;
  }

  void updateAll() {
    for(typename std::set<T*>::iterator it = listeners.begin(); it != listeners.end(); ++it) {
      *(*it) = data;
    }
  }

  T data;
  std::set<T*> listeners;

  void print_value(FILE* out) {
    fprintf(out, "%g", double(data));
  }
};

typedef Property<double> ScalarProperty;

class VectorProperty : public Property<double*>
{
public:
  int cols;

  VectorProperty(const int N) :
    cols(N)
  {
    data = new double[N];
    for(int col = 0; col < N; col++) {
      data[col] = 0.0;
    }
  }

  virtual ~VectorProperty(){
    delete [] data;
  }

  void print_value(FILE* out) {
    fprintf(out, "[");
    for(int col = 1; col < cols; col++) {
      fprintf(out, "%g", data[col]);
      if((col+1) < cols) fprintf(out, " ");
    }
    fprintf(out, "]");
  }
};

class MatrixProperty : public Property<double**>
{
public:
  int rows;
  int cols;

  MatrixProperty(const int N, const int M) :
    rows(N),
    cols(M)
  {
    double * array = new double[N*M];
    data = new double*[N];

    for(int row = 0; row < N; row++) {
      data[row] = &array[row*M];

      for(int col = 0; col < M; col++) {
        data[row][col] = 0.0;
      }
    }
  }

  virtual ~MatrixProperty() {
    delete [] data[0];
    delete [] data;
  }

  void print_value(FILE* out) {
    if (!out)  // out == NULL
      return;

    fprintf(out, "[");
    for(int row = 1; row < rows; row++) {
      for(int col = 1; col < cols; col++) {
        fprintf(out, "%g", data[row][col]);
        if((col+1) < cols) fprintf(out, " ");
      }
      if((row+1) < rows) fprintf(out, "; ");
    }
    fprintf(out, "]");
  }
};

// -------------------------------------------------------------------

typedef ScalarProperty* (*ScalarPropertyCreator)(PropertyRegistry & registry, const char * caller, bool sanity_checks);
typedef VectorProperty* (*VectorPropertyCreator)(PropertyRegistry & registry, const char * caller, bool sanity_checks);
typedef MatrixProperty* (*MatrixPropertyCreator)(PropertyRegistry & registry, const char * caller, bool sanity_checks);

class PropertyRegistry : protected Pointers {
  Properties properties;

public:
  PropertyRegistry(LAMMPS* lmp);
  ~PropertyRegistry();
  int max_type();
  double min_radius();
  double max_radius();
  LAMMPS * getLAMMPS();

  FixPropertyGlobal* getGlobalProperty(const char *varname, const char *style, const char *svmstyle, int len1, int len2, const char *caller);

  ScalarProperty * getScalarProperty(std::string varname,const char *caller);
  VectorProperty * getVectorProperty(std::string varname,const char *caller);
  MatrixProperty * getMatrixProperty(std::string varname,const char *caller);

  void registerProperty(std::string varname, ScalarPropertyCreator creator, bool sanity_checks = true);
  void registerProperty(std::string varname, VectorPropertyCreator creator, bool sanity_checks = true);
  void registerProperty(std::string varname, MatrixPropertyCreator creator, bool sanity_checks = true);

  void connect(std::string varname, double ** & variable, const char *caller);
  void connect(std::string varname, double * & variable, const char *caller);
  void connect(std::string varname, double & variable, const char *caller);

  void init();

  void print_all(FILE* out);

private:
  std::map<std::string, ScalarPropertyCreator> scalar_creators;
  std::map<std::string, VectorPropertyCreator> vector_creators;
  std::map<std::string, MatrixPropertyCreator> matrix_creators;

  std::map<std::string, ScalarProperty*> scalars;
  std::map<std::string, VectorProperty*> vectors;
  std::map<std::string, MatrixProperty*> matrices;

  std::map<std::string, bool> use_sanity_checks;
};

#endif /* PROPERTY_REGISTRY_H_ */
