#ifndef ISPD_PARSE_H
#define ISPD_PARSE_H

#include "parser_helper.h"
#include "fstream"
#include "iostream"
#include "map"
#include "vector"
#include "circuit_class.h"
#include "gate_class.h"

using namespace std;

void ispd_parse_verilog(string filename, circuit_class *ckt);
void ispd_parse_lib(string filename, map<string, double> &library );

void generate_placement( circuit_class *ckt, string benchmark, map<string, double> &library );
void write_nets_file(ofstream& fp, circuit_class *ckt);
void write_nodes_file(ofstream& fp, circuit_class *ckt, map<string, double> &library);
void write_wts_file(int ext, circuit_class *ckt , ofstream& fp);

#endif
