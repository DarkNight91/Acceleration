#ifndef PARSER_H
#define PARSER_H

#include "circuit_class.h"

int    parse_netlist(circuit_class* , const char*);
void parse_placement_file(const char* placement_file_name, circuit_class *ckt, float temp_max_min[]);


#endif
