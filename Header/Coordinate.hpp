#ifndef COORDINATE_H
#define COORDINATE_H

#include "AxialData.hpp"

/* 생성자 */
Coordinate::Coordinate ():dim (0) {
  this->comp = NULL;
}

Coordinate::Coordinate (int dms):dim (dms) {
  this->comp = new double[dms];

  if (!this->comp) {
    printf ("no memory for ");
    printf ("%d", this->Dimension ());
    printf ("D coordinate...\n");
    exit (1);
  }

  for (size_t i = 0; i < this->Dimension (); i++) this->comp[i] = ZeroValue;
}

Coordinate::Coordinate (double x, double y):dim (D2) {
  this->comp = new double[D2];

  if (!this->comp) {
    printf ("no memory for ");
    printf ("%d", this->Dimension ());
    printf ("D coordinate...\n");
    exit (1);
  }

  this->comp[0] = x;
  this->comp[1] = y;
}

Coordinate::Coordinate (double x, double y, double z):dim (D3) {
  this->comp = new double[D3];

  if (!this->comp) {
    printf ("no memory for ");
    printf ("%d", this->Dimension ());
    printf ("D coordinate...\n");
    exit (1);
  }

  this->comp[0] = x;
  this->comp[1] = y;
  this->comp[2] = z;
}

Coordinate::Coordinate (const Coordinate & src) {
  this->dim = src.Dimension ();

  this->comp = new double[this->Dimension ()];

  if (!this->comp) {
    printf ("no memory for ");
    printf ("%d", this->Dimension ());
    printf ("D coordinate...\n");
    exit (1);
  }

  for (size_t i = 0; i < this->Dimension (); i++) this->comp[i] = src[i];
}

Coordinate::~Coordinate () {
  if (this->comp) delete [] this->comp;
}

Coordinate & Coordinate::SetCoordinate (char coord, double value) {
  char errorMassage[1024];

  if (this->Dimension () == 2) {
    if (coord == 'x') {this->comp[0] = value; return *this;}
    if (coord == 'y') {this->comp[1] = value; return *this;}
    sprintf (errorMassage, "Coordinate::SetCoordinate error, dim = %d, coord = %c, Only x and y can be referenced.", this->Dimension (), coord);
    PrintError (errorMassage);
    exit (1);
  }

  if (this->Dimension () == 3) {
    if (coord == 'x') {this->comp[0] = value; return *this;}
    if (coord == 'y') {this->comp[1] = value; return *this;}
    if (coord == 'z') {this->comp[2] = value; return *this;}
    sprintf (errorMassage, "Coordinate::SetCoordinate error, dim = %d, coord = %c, Only x, y and z can be referenced.", this->Dimension (), coord);
    PrintError (errorMassage);
    exit (1);
  }

  sprintf (errorMassage, "Coordinate::SetCoordinate error, dim = %d, Please check dimension.", this->Dimension ());
  PrintError (errorMassage);
  exit (1);
}

Coordinate & Coordinate::SetCoordinate (double x, double y) {
  char errorMassage[1024];

  if (this->Dimension () != 2) {
    sprintf (errorMassage, "Coordinate::SetCoordinate error, dim = %d, Only x and y can be referenced.", this->Dimension ());
    PrintError (errorMassage);
    exit (1);
  }

  this->comp[0] = x;
  this->comp[1] = y;

  return *this;
}

Coordinate & Coordinate::SetCoordinate (double x, double y, double z) {
  char errorMassage[1024];

  if (this->Dimension () != 3) {
    sprintf (errorMassage, "Coordinate::SetCoordinate error, dim = %d, Only x, y and z can be referenced.", this->Dimension ());
    PrintError (errorMassage);
    exit (1);
  }

  this->comp[0] = x;
  this->comp[1] = y;
  this->comp[2] = z;
  return *this;
}

Coordinate & Coordinate::SetCoordinate (Coordinate & src) {
  this->dim = src.Dimension ();
  this->comp = new double[this->Dimension ()];

  if (!this->comp) {
    printf ("no memory for ");
    printf ("%d", this->Dimension ());
    printf ("D coordinate...\n");
    exit (1);
  }

  for (size_t i = 0; i < this->Dimension (); i++) this->comp[i] = src[i];

  return *this;
}

double & Coordinate::Value (char coord) const {
  char errorMassage[1024];

  if (this->Dimension () == 2) {
    if (coord == 'x') return this->comp[0];
    if (coord == 'y') return this->comp[1];
    sprintf (errorMassage, "Coordinate::Value error, dim = %d, coord = %c, Only x and y can be referenced.", this->Dimension (), coord);
    PrintError (errorMassage);
    exit (1);
  }

  if (this->Dimension () == 3) {
    if (coord == 'x') return this->comp[0];
    if (coord == 'y') return this->comp[1];
    if (coord == 'z') return this->comp[2];
    sprintf (errorMassage, "Coordinate::Value error, dim = %d, coord = %c, Only x, y and z can be referenced.", this->Dimension (), coord);
    PrintError (errorMassage);
    exit (1);
  }

  sprintf (errorMassage, "Coordinate::Value error, dim = %d, Please check dimension.", this->Dimension ());
  PrintError (errorMassage);
  exit (1);
}

double Coordinate::Distance (const Coordinate &src) {
  char errorMassage[1024];
  double value = ZeroValue;

  if (this->Dimension () != src.Dimension ()) {
    sprintf (errorMassage, "Coordinate::Distance error, this dimension = %d, but source dimension = %d", this->Dimension (), src.Dimension ());
    PrintError (errorMassage);
    exit (1);
  }

  for (size_t i = 0; i < this->Dimension (); i++) {
    value += (this->comp[i] - src[i]) * (this->comp[i] - src[i]);
  }
  value = sqrt (value);

  return value;
}


double & Coordinate::operator[] (const int idx) {
  char errorMassage[1024];

  if (idx >= 0 && idx < this->Dimension ()) return this->comp[idx];
  else {
    sprintf (errorMassage, "Coordiante operator[] error, dim = %d, Please check diminsion", this->Dimension ());
    PrintError (errorMassage);
    exit (1);
  }
}

double & Coordinate::operator[] (const int idx) const {
  char errorMassage[1024];

  if (idx >= 0 && idx < this->Dimension ()) return this->comp[idx];
  else {
    sprintf (errorMassage, "Coordiante operator[] error, dim = %d, Please check diminsion", this->Dimension ());
    PrintError (errorMassage);
    exit (1);
  }
}

Coordinate & Coordinate::operator= (const Coordinate & src) {

  if (this == &src) return *this;

  if (src.Dimension () != this->Dimension ()) {

    if (comp) {
      delete [] comp;
      comp = NULL;
    }

    this->dim = src.Dimension ();
    this->comp = new double[src.Dimension ()];

    if (!comp) {
      cout << "no memory for " << this->Dimension () << "D vector..." << endl;
      exit(1);
    }
  }

  for (size_t i = 0; i < src.Dimension (); i++) this->comp[i] = src[i];

  return *this;
}

double Coordinate::operator- (const Coordinate & src) {

  return this->Distance (src);
}

Coordinate & Coordinate::ShowContents (const char *front, const char *end) {

  printf ("%s ", front);

  printf ("[%dD]", this->Dimension ());

  printf ("<");
  for (size_t i = 0; i < this->Dimension (); i++) {
    printf ("%23.16e", this->comp[i]);
    if (i < this->Dimension () - 1) printf (", ");
  }
  printf (">");

  printf ("%s", end);

  return * this;
}

bool Coordinate::operator== (const Coordinate & src) {

  if (this->Dimension () != src.Dimension ()) return false;

  for (size_t i = 0; i < this->Dimension (); i++) if (!IsEqualDouble (this->comp[i], src[i])) return false;

  return true;
}

Coordinate & Coordinate::ShowContents (const char *str) {

  printf ("%s ", str);

  printf ("[%dD]", this->Dimension ());

  printf ("<");
  for (size_t i = 0; i < this->Dimension (); i++) {
    printf ("%23.16e", this->comp[i]);
    if (i < this->Dimension () - 1) printf (", ");
  }
  printf (">");
  printf ("\n");

  return * this;
}

#endif
