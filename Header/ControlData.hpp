#ifndef CONTROLDATA_H
#define CONTROLDATA_H

#include "util.hpp"

// control data 읽기
ControlData & ControlData::LoadCtrlData (string  AGL_output_file, string AGM_output_file, double initialTime, double terminalTime, int timeStep, double dt) {
  // output file 이름
  string output;

  // Axial line generator에서 만들어진 파일의 이름
  this->AxialFile  = AGL_output_file;
  // output file의 이름
  output     = AGM_output_file;
  // solution파일의 이름
  this->output_sol = output;
  // 미분파일의 이름
  this->output_del = output + "_del.dat";
  // 경계조건파일의 이름
  this->boundary = output + "_cond.dat";
  // x-축선파일의 이름
  this->xaxial = output + "_xaxial.dat";
  // y-축선파일의 이름
  this->yaxial = output + "_yaxial.dat";
  //
  this->initialTime = initialTime;
  //
  this->terminalTime = terminalTime;
  //
  this->timeStep = timeStep;
  //
  this->dt = dt;
  return *this;
}

// control data 출력
ControlData & ControlData::ShowCtrlData () {

  printf("======= Control data =======\n");
  printf("     AxialFile = %s\n", this->Axialfile ().c_str ());
  printf("    output_sol = %s\n", this->Output_Sol ().c_str ());
  printf("    output_del = %s\n", this->Output_Del ().c_str ());
  printf("   initialTime = %23.16e\n", this->InitialTime ());
  printf("  terminalTime = %23.16e\n", this->TerminalTime ());
  printf("            dt = %23.16e\n", this->Dt ());
  printf("      timeStep = %d\n", this->TimeStep ());
  printf("============================\n\n");

  return *this;
}

#endif
