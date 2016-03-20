//////////////////////////////////////////////////////////////////
//
//
//  Helper functions and classes to parse the ISPD 2013 contest
//  benchmark files.
//
//  This code is provided for description purposes only. The contest
//  organizers cannot guarantee that the provided code is free of
//  bugs or defects. !!!! USE THIS CODE AT YOUR OWN RISK !!!!!
//
//
//  The contestants are free to use these functions as-is or make
//  modifications. If the contestants choose to use the provided
//  code, they are responsible for making sure that it works as
//  expected.
//
//  The code provided here has no real or implied warranties.
//
//
////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <limits>
#include <cassert>
#include <cctype>

#include "Celllib.h"
#include "werc_102_vlsi.h"
#include "parser_helper.h"

bool is_special_char (char c) {

  static const char specialChars[] = {'(', ')', ',', ':', ';', '/', '#', '[', ']', '{', '}', '*', '\"', '\\'} ;
  
  for (int i=0; i < sizeof(specialChars); ++i) {
    if (c == specialChars[i])
      return true ;
  }

  return false ;
}

// Read the next line and return it as a list of tokens skipping white space and special characters
// The return value indicates success/failure.
bool read_line_as_tokens (istream& is, vector<string>& tokens,
			  bool includeSpecialChars) {

  tokens.clear() ;
  
  string line ;
  std::getline (is, line) ;

  while (is && tokens.empty()) {
    
    string token = "" ;

    for (int i=0; i < line.size(); ++i) {
      char currChar = line[i] ;
      bool isSpecialChar = is_special_char(currChar) ;

      if (std::isspace (currChar) || isSpecialChar) {

        if (!token.empty()) {
          // Add the current token to the list of tokens
          tokens.push_back(token) ;
          token.clear() ;
        }

	if (includeSpecialChars && isSpecialChar) {
	  tokens.push_back(string(1, currChar)) ;
	}

      } else {
        // Add the char to the current token
        token.push_back(currChar) ;
      }

    }

    if (!token.empty())
      tokens.push_back(token) ;
    

    if (tokens.empty())
      // Previous line read was empty. Read the next one.
      std::getline (is, line) ;    
  }

  //for (int i=0; i < tokens.size(); ++i)
  //  cout << tokens[i] << " " ;
  //cout << endl ;
  
  return !tokens.empty() ;
}



bool VerilogParser::read_module (string& moduleName) {

  vector<string> tokens ;
  bool valid = read_line_as_tokens (is, tokens) ;

  while (valid) {

    if (tokens.size() == 2 && tokens[0] == "module") {
      moduleName = tokens[1] ;

      break ;
    }

    valid = read_line_as_tokens (is, tokens) ;
  }

  // Read and skip the port names in the module definition
  // until we encounter the tokens {"Start", "PIs"}
  while (valid && !(tokens.size() == 2 && tokens[0] == "Start" && tokens[1] == "PIs")) {

    valid = read_line_as_tokens (is, tokens) ;
    assert (valid) ;
  }

  return valid ;
}

bool VerilogParser::read_primary_input (string& primaryInput) {

  primaryInput = "" ;
  
  vector<string> tokens ;
  bool valid = read_line_as_tokens (is, tokens) ;

  assert (valid) ;
  assert (tokens.size() == 2) ;

  if (valid && tokens[0] == "input") {
    primaryInput = tokens[1] ;

  } else {
    assert (tokens[0] == "Start" && tokens[1] == "POs") ;
    return false ;
  }

  return valid ; 
}



bool VerilogParser::read_primary_output (string& primaryOutput) {

  primaryOutput = "" ;
  
  vector<string> tokens ;
  bool valid = read_line_as_tokens (is, tokens) ;

  assert (valid) ;
  assert (tokens.size() == 2) ;

  if (valid && tokens[0] == "output") {
    primaryOutput = tokens[1] ;

  } else {
    assert (tokens[0] == "Start" && tokens[1] == "wires") ;
    return false ;
  }

  return valid ; 
}


bool VerilogParser::read_wire (string& wire) {

  wire = "" ;
  
  vector<string> tokens ;
  bool valid = read_line_as_tokens (is, tokens) ;

  assert (valid) ;
  assert (tokens.size() == 2) ;

  if (valid && tokens[0] == "wire") {
    wire = tokens[1] ;

  } else {
    assert (tokens[0] == "Start" && tokens[1] == "cells") ;
    return false ;
  }

  return valid ; 
}


bool VerilogParser::read_cell_inst (string& cellType, string& cellInstName,
                                    vector<std::pair<string, string> >& pinNetPairs) {

  cellType = "" ;
  cellInstName = "" ;
  pinNetPairs.clear() ;

  vector<string> tokens ;
  bool valid = read_line_as_tokens (is, tokens) ;

  assert (valid) ;

  if (tokens.size() == 1) {
    assert (tokens[0] == "endmodule") ;
    return false ;
  }

  assert (tokens.size() >= 4) ; // We should have cellType, instName, and at least one pin-net pair

  cellType = tokens[0] ;
  cellInstName = tokens[1] ;

  for (int i=2; i < tokens.size()-1; i += 2) {

    assert (tokens[i][0] == '.') ; // pin names start with '.'
    string pinName = tokens[i].substr(1) ; // skip the first character of tokens[i] 
    
    pinNetPairs.push_back(std::make_pair(pinName, tokens[i+1])) ;
  }

  return valid ;
}

// No need to parse the 3D LUTs, because they will be ignored
void LibParser::_skip_lut_3D () {

    std::vector<string> tokens ;

    bool valid = read_line_as_tokens (is, tokens) ;
    assert (valid) ;
    assert (tokens[0] == "index_1") ;
    assert (tokens.size() >= 2) ;
    int size1 = tokens.size() - 1 ;

    valid = read_line_as_tokens (is, tokens) ;
    assert (valid) ;
    assert (tokens[0] == "index_2") ;
    assert (tokens.size() >= 2) ;
    int size2 = tokens.size() - 1 ;  

    valid = read_line_as_tokens (is, tokens) ;
    assert (valid) ;
    assert (tokens[0] == "index_3") ;
    assert (tokens.size() >= 2) ;
    int size3 = tokens.size() - 1 ;  

    valid = read_line_as_tokens (is, tokens) ;
    assert (valid) ;
    assert (tokens.size() == 1 && tokens[0] == "values") ;

    for (int i=0; i < size1; ++i) {
        for (int j=0; j < size2; ++j) {

            valid = read_line_as_tokens (is, tokens) ;
            assert (valid) ;
            assert (tokens.size() == size3) ;
        }
    }

}

void LibParser::_begin_read_lut (LibParserLUT& lut) {

    std::vector<string> tokens ;
    bool valid = read_line_as_tokens (is, tokens) ;

    assert (valid) ;
    assert (tokens[0] == "index_1") ;
    assert (tokens.size() >= 2) ;

    int size1 = tokens.size()-1 ;
    lut.loadIndices.resize(size1) ;
    for (int i=0; i < tokens.size()-1; ++i) {

        lut.loadIndices[i] = std::atof(tokens[i+1].c_str()) ;
    }

    valid = read_line_as_tokens (is, tokens) ;

    assert (valid) ;
    assert (tokens[0] == "index_2") ;
    assert (tokens.size() >= 2) ;

    int size2 = tokens.size()-1 ;
    lut.transitionIndices.resize(size2) ;
    for (int i=0; i < tokens.size()-1; ++i) {

        lut.transitionIndices[i] = std::atof(tokens[i+1].c_str()) ;
    }

    valid = read_line_as_tokens (is, tokens) ;
    assert (valid) ;
    assert (tokens.size() == 1 && tokens[0] == "values") ;

    lut.tableVals.resize(size1) ;
    for (int i=0 ; i < lut.loadIndices.size(); ++i) {
        valid = read_line_as_tokens (is, tokens) ;
        assert (valid) ;
        assert (tokens.size() == lut.transitionIndices.size()) ;

        lut.tableVals[i].resize(size2) ;
        for (int j=0; j < lut.transitionIndices.size(); ++j) {
            lut.tableVals[i][j] = std::atof(tokens[j].c_str()) ;

        }
    }


}

void LibParser::_begin_read_timing_info (string toPin, LibParserTimingInfo& timing) {

    timing.toPin = toPin ;

    bool finishedReading = false ;

    std::vector<string> tokens ;  
    while (!finishedReading) {

        bool valid = read_line_as_tokens (is, tokens) ;
        assert (valid) ;
        assert (tokens.size() >= 1) ;

        if (tokens[0] == "cell_fall") {
            _begin_read_lut (timing.fallDelay) ;

        } else if (tokens[0] == "cell_rise") {
            _begin_read_lut (timing.riseDelay) ;

        } else if (tokens[0] == "fall_transition") {
            _begin_read_lut (timing.fallTransition) ;

        } else if (tokens[0] == "rise_transition") {
            _begin_read_lut (timing.riseTransition) ;

        } else if (tokens[0] == "fall_constraint") {

            _skip_lut_3D() ; // will ignore fall constraints      

        } else if (tokens[0] == "rise_constraint") {

            _skip_lut_3D() ; // will ignore rise constraints

        } else if (tokens[0] == "timing_sense") {
            timing.timingSense = tokens[1] ;

        } else if (tokens[0] == "related_pin") {

            assert (tokens.size() == 2) ;
            timing.fromPin = tokens[1] ;

        } else if (tokens[0] == "End") {

            assert (tokens.size() == 2) ;
            assert (tokens[1] == "timing") ;
            finishedReading = true ;

        } else if (tokens[0] == "timing_type") {
            // ignore data

        } else if (tokens[0] == "related_output_pin") {
            // ignore data

        } else {

            cout << "Error: Unknown keyword: " << tokens[0] << endl ;
            assert (false) ; // unknown keyword
        }

    }


}


void LibParser::_begin_read_pin_info (string pinName, LibParserCellInfo& cell, LibParserPinInfo& pin) {

    pin.name = pinName ;
    pin.isClock = false ;
    pin.maxCapacitance = std::numeric_limits<double>::max() ;

    bool finishedReading = false ;

    std::vector<string> tokens ;  
    while (!finishedReading) {

        bool valid = read_line_as_tokens (is, tokens) ;
        assert (valid) ;
        assert (tokens.size() >= 1) ;

        if (tokens[0] == "direction") {

            assert (tokens.size() == 2) ;
            if (tokens[1] == "input")
                pin.isInput = true ;
            else if (tokens[1] == "output")
                pin.isInput = false ;
            else
                assert (false) ; // undefined direction

        } else if (tokens[0] == "capacitance") {

            assert (tokens.size() == 2) ;
            pin.capacitance = std::atof(tokens[1].c_str()) ;

        } else if (tokens[0] == "max_capacitance") {

            assert (tokens.size() == 2) ;
            pin.maxCapacitance = std::atof(tokens[1].c_str()) ;


        } else if (tokens[0] == "timing") {

            cell.timingArcs.push_back(LibParserTimingInfo()) ; // add an empty TimingInfo object
            _begin_read_timing_info (pinName, cell.timingArcs.back()) ; // pass the empty object to the function to be filled

        } else if (tokens[0] == "clock") {

            pin.isClock = true ;

        } else if (tokens[0] == "End") {

            assert (tokens.size() == 2) ;
            assert (tokens[1] == "pin") ;
            finishedReading = true ;

        } else if (tokens[0] == "function") {

            // ignore data

        } else if (tokens[0] == "min_capacitance") {

            // ignore data

        } else if (tokens[0] == "nextstate_type") {

            // ignore data

        } else {
            cout << "Error: Unknown keyword: " << tokens[0] << endl ;      
            assert (false) ; // unknown keyword 

        }

    }


}

void LibParser::_begin_read_cell_info (string cellName, LibParserCellInfo& cell) {

    cell.name = cellName ;
    cell.isSequential = false ;
    cell.dontTouch = false ;

    bool finishedReading = false ;

    std::vector<string> tokens ;  
    while (!finishedReading) {

        bool valid = read_line_as_tokens (is, tokens) ;
        assert (valid) ;
        assert (tokens.size() >= 1) ;

        if (tokens[0] == "cell_leakage_power") {

            assert (tokens.size() == 2) ;
            cell.leakagePower = std::atof(tokens[1].c_str()) ;

        } else if (tokens[0] == "cell_footprint") {

            assert (tokens.size() == 2) ;
            cell.footprint = tokens[1] ;

        } else if (tokens[0] == "area") {

            assert (tokens.size() == 2) ;
            cell.area = std::atof(tokens[1].c_str()) ;

        } else if (tokens[0] == "clocked_on") {

            cell.isSequential = true ;

        } else if (tokens[0] == "dont_touch") {

            cell.dontTouch = true ;

        } else if (tokens[0] == "pin") {

            assert (tokens.size() == 2) ;

            cell.pins.push_back(LibParserPinInfo()) ; // add empty PinInfo object
            _begin_read_pin_info (tokens[1], cell, cell.pins.back()) ; // pass the new PinInfo object to be filled

        } else if (tokens[0] == "End") {

            assert (tokens.size() == 3) ;
            assert (tokens[1] == "cell") ;
            assert (tokens[2] == cellName) ;
            finishedReading = true ;

        } else if (tokens[0] == "cell_footprint") {

            // ignore data

        } else if (tokens[0] == "ff") {

            // ignore data

        } else if (tokens[0] == "next_state") {

            // ignore data

        } else if (tokens[0] == "dont_use") {

            // ignore data

        } else {

            cout << "Error: Unknown keyword: " << tokens[0] << endl ;
            assert (false) ; // unknown keyword
        }
    } 

}


// Read the default max_transition defined for the library.
// Return value indicates if the last read was successful or not.  
// This function must be called in the beginning before any read_cell_info function call.
bool LibParser::read_default_max_transition (double& maxTransition) {

    maxTransition = 0.0 ;
    vector<string> tokens ;

    bool valid = read_line_as_tokens (is, tokens) ;

    while (valid) {

        if (tokens.size() == 2 && tokens[0] == "default_max_transition") {
            maxTransition = std::atof(tokens[1].c_str()) ;
            return true ;
        }

        valid = read_line_as_tokens (is, tokens) ;
    }

    return false ;
}



// Read the next standard cell definition.
// Return value indicates if the last read was successful or not.  
bool LibParser::read_cell_info (LibParserCellInfo& cell) {

    vector<string> tokens ;
    bool valid = read_line_as_tokens (is, tokens) ;


    while (valid) {

        if (tokens.size() == 2 && tokens[0] == "cell") {
            _begin_read_cell_info (tokens[1], cell) ;

            return true ;
        }

        valid = read_line_as_tokens (is, tokens) ;
    }

    return false ;
}

ostream& operator<< (ostream& os, LibParserLUT& lut) {

    if (lut.loadIndices.empty() && lut.transitionIndices.empty() && lut.tableVals.empty())
        return os ;

    // We should have either all empty or none empty.
    assert (!lut.loadIndices.empty() && !lut.transitionIndices.empty() && !lut.tableVals.empty()) ;

    assert (lut.tableVals.size() == lut.loadIndices.size()) ;
    assert (lut.tableVals[0].size() == lut.transitionIndices.size()) ;

    cout << "\t" ;
    for (int i=0; i < lut.transitionIndices.size(); ++i) {
        cout << lut.transitionIndices[i] << "\t" ;
    }
    cout << endl ;


    for (int i=0; i < lut.loadIndices.size(); ++i) {
        cout << lut.loadIndices[i] << "\t" ;

        for (int j=0; j < lut.transitionIndices.size(); ++j)
            cout << lut.tableVals[i][j] << "\t" ;

        cout << endl ;

    }

    return os ;
}


ostream& operator<< (ostream& os, LibParserTimingInfo& timing) {

    cout << "Timing info from " << timing.fromPin << " to " << timing.toPin << ": " << endl ;
    cout << "Timing sense: " << timing.timingSense << endl ;

    cout << "Fall delay LUT: " << endl ;
    cout << timing.fallDelay ;

    cout << "Rise delay LUT: " << endl ;
    cout << timing.riseDelay ;

    cout << "Fall transition LUT: " << endl ;
    cout << timing.fallTransition ;

    cout << "Rise transition LUT: " << endl ;
    cout << timing.riseTransition ;

    return os ;
}


ostream& operator<< (ostream& os, LibParserPinInfo& pin) {

    cout << "Pin " << pin.name << ":" << endl ;
    cout << "capacitance: " << pin.capacitance << endl ;
    cout << "maxCapacitance: " << pin.maxCapacitance << endl ;
    cout << "isInput? " << (pin.isInput ? "true" : "false") << endl ;
    cout << "isClock? " << (pin.isClock ? "true" : "false") << endl ;
    cout << "End pin" << endl ;

    return os ;
}


ostream& operator<< (ostream& os, LibParserCellInfo& cell) {

    cout << "Library cell " << cell.name << ": " << endl ;

    cout << "Footprint: " << cell.footprint << endl ;
    cout << "Leakage power: " << cell.leakagePower << endl ;
    cout << "Area: " << cell.area << endl ;
    cout << "Sequential? " << (cell.isSequential ? "yes" : "no") << endl ;
    cout << "Dont-touch? " << (cell.dontTouch ? "yes" : "no") << endl ;      

    cout << "Cell has " << cell.pins.size() << " pins: " << endl ;
    for (int i=0; i < cell.pins.size(); ++i) {
        cout << cell.pins[i] << endl ;
    }

    cout << "Cell has " << cell.timingArcs.size() << " timing arcs: " << endl ;
    for (int i=0; i < cell.timingArcs.size(); ++i) {
        cout << cell.timingArcs[i] << endl ;
    }

    cout << "End of cell " << cell.name << endl << endl ;

    return os ;
}


// Example function that uses VerilogParser class to parse the given ISPD-13 verilog
// file. The extracted data is simply printed out in this example.
void test_verilog_parser (string filename) {

  VerilogParser vp (filename) ;

  string moduleName ;
  bool valid = vp.read_module (moduleName) ;
  assert (valid) ;
  
  cout << "Module " << moduleName << endl << endl ;

  do {
    string primaryInput ;
    valid = vp.read_primary_input (primaryInput) ;

    if (valid)
      cout << "Primary input: " << primaryInput << endl ;

  } while (valid) ;

  cout << endl ;

  do {
    string primaryOutput ;
    valid = vp.read_primary_output (primaryOutput) ;

    if (valid)
      cout << "Primary output: " << primaryOutput << endl ;

  } while (valid) ;

  cout << endl ;

  do {
    string net ;
    valid = vp.read_wire (net) ;

    if (valid)
      cout << "Net: " << net << endl ;

  } while (valid) ;


  cout << endl ;
  cout << "Cell insts: " << std::endl ;
  
  do {
    string cellType, cellInst ;
    vector<std::pair<string, string> > pinNetPairs ;
    
    valid = vp.read_cell_inst (cellType, cellInst, pinNetPairs) ;

    if (valid) {
      cout << cellType << " " << cellInst << " " ;
      for (int i=0; i < pinNetPairs.size(); ++i) {
        cout << "(" << pinNetPairs[i].first << " " << pinNetPairs[i].second << ") " ;
      }

      cout << endl ;
    }
    
  } while (valid) ;
  
  
}
// Example function that uses LibParser class to parse the given ISPD-13 lib
// file. The extracted data is simply printed out in this example.
void test_lib_parser (string filename) {

  LibParser lp (filename) ;

  double maxTransition = 0.0 ;
  bool valid = lp.read_default_max_transition(maxTransition) ;

  assert (valid) ;
  cout << "The default max transition defined is " << maxTransition << endl ;

  int readCnt = 0 ;
  do {
    LibParserCellInfo cell ;
    valid = lp.read_cell_info (cell) ;

    if (valid) {
      ++readCnt ;

      cout << cell << endl ;
    }
    
  } while (valid) ;

  cout << "Read " << readCnt << " number of library cells" << endl ;
}

VOID t_mmlparser::mp_parse_param_name(const UINT16D index)
{
    const CHAR *p = NULL;
    t_param tmp;//assert( mp != NULL ); 
    if (NULL == cur)
    {
        return;
    }

    if (p = strstr(cur + 1, "("))
    {
        //strncpy(mp->param[index].paramname, mp->cur+1, p - mp->cur -1 ); 
        //param[index].paramname.assign(cur+1, p - cur -1);
        tmp.paramname.assign(cur + 1, p - cur - 1);
        param.push_back(tmp);
    }
    else if (p = strstr(cur + 1, ","))
    {
        //strncpy(mp->param[index].paramname, mp->cur+1, p - mp->cur -1 ); 
        //param[index].paramname.assign(cur+1, p - cur -1 );
        tmp.paramname.assign(cur + 1, p - cur - 1);
        param.push_back(tmp);
        if (*(p + 1) == ' ')//delete the possible backward space  
        {
            p++;
        }
    }
    else if (p = strstr(cur + 1, ")"))
    {
        //strncpy(mp->param[index].paramname, mp->cur+1, p - mp->cur -1 ); 
        //param[index].paramname.assign(cur+1, p - cur -1);
        tmp.paramname.assign(cur + 1, p - cur - 1);
        param.push_back(tmp);
    }
    else if (p = strstr(cur + 1, "+"))
    {
        cur = p;
        p = strstr(cur + 1, ":");
        tmp.paramname.assign(cur + 1, p - cur - 1);
        //tmp.paramname.replace(" ", "");seems that atof will automatically replace the " "
        param.push_back(tmp);
        --p;//To move p one byte backward from ":", makes the next branch successful 
    }
    else if (p = strstr(cur + 1, ":"))
    {
        unsigned pos = p - mml.c_str();
        //get string from ":" to the end, wish no non-number character
        tmp.paramname = mml.substr(pos + 1);
        param.push_back(tmp);
    }
    else
    {
        H102_ASSERT(NULL == p);
        cur = p;
        return;
    }

    cur = p;
}


VOID t_mmlparser::mp_parse_param()
{
    UINT16D  index = 0;
    while (cur)
    {
        mp_parse_param_name(index);
        //mp_parse_param_value( mp, index ); 
        index++;
    }
}

VOID t_mmlparser::mp_parse_topo_name()
{
    const CHAR  *p = NULL;
    const CHAR  *n = NULL;

    /* Can only parse a specific line that before = or ( there is cmdname */
    if (n = p = strstr(mml.c_str(), "="))//n point to '='
    {
        if (*(p - 1) == ' ')//delete the possible forward space 
        {
            p--;
        }
        //strncpy(cmdname.data(), mml.c_str(), p - mml.c_str());
        cmdname.assign(mml.c_str(), p - mml.c_str());
        if (*(n + 1) == ' ')
        {
            p = n + 1;
        }//delete the possible backward space
        else
        {
            p = n;
        }
    }
    else if (p = strstr(mml.c_str(), "("))
    {
        //strncpy(cmdname.c_str(), mml.c_str(), p - mml.c_str());
        cmdname.assign(mml.c_str(), p - mml.c_str());
    }
    else
    {
        H102_ASSERT(NULL == p);
        cmdname = mml;
    }

    cur = p;
}

VOID t_mmlparser::mp_parse_topo(const CHAR* in_mml)
{
    mml = in_mml;
    mp_parse_topo_name();
    mp_parse_param();
}


VOID    begin_read_cell_param (istream &is, LIB_LUT& cell)
{
    //given a gate size, threshold voltage we can get (1)R (2)C (3) offset (4)power (5) area
#define MAX_CELL_CONTENT   5

    UINT8D   num = 0;
    UINT8D   index;
    bool    finishedReading = false ;
    LUT_PARAM   gate_param;
    cell.resize(VT_LEVEL_MAX);

    std::vector<string> tokens ;  
    while (!finishedReading)
    {
        bool valid = read_line_as_tokens (is, tokens);
        H102_ASSERT (valid) ;
        H102_ASSERT (tokens.size() >= 1) ;

        if (tokens[0] == "END")
        {
            H102_ASSERT(tokens.size() == 1) ;
            finishedReading = true ;
            H102_ASSERT(VT_LEVEL_MAX == num);
        }
        else
        {
            H102_ASSERT(MAX_CELL_CONTENT*PN_SIZE_MAX == tokens.size());
            for (index = 0; index < MAX_CELL_CONTENT*PN_SIZE_MAX;)
            {
                gate_param.cap = std::atof(tokens[index++].c_str());
                gate_param.resist = std::atof(tokens[index++].c_str());
                gate_param.offset = std::atof(tokens[index++].c_str());
                gate_param.leakpower = std::atof(tokens[index++].c_str());
                gate_param.area = std::atof(tokens[index++].c_str());
                cell[num].push_back(gate_param);
            }
            ++num;//totally there are three rows
        }
    }
    //cell_set[cellName] = cell;
}


#if 0
int main(int args, char** argv) {

  if (args != 3) {
    cout << "Usage: " << argv[0] << " <format> <file>" << endl ;
    exit(0) ;
  }

  if (argv[1] == string("verilog")) {
    cout << "Parsing verilog..." << endl ;
    test_verilog_parser(argv[2]) ;
    cout << "Finished parsing verilog." << endl ;    
  } 
  else if (argv[1] == string("sdc")) {
    cout << "Parsing sdc..." << endl ;  
    test_sdc_parser (argv[2]) ;
    cout << "Finished parsing sdc." << endl ;
  }
  else if (argv[1] == string("spef")) {
    cout << "Parsing spef..." << endl ;  
    test_spef_parser (argv[2]) ;
    cout << "Finished parsing spef." << endl ;
  }
  else if (argv[1] == string("timing")) {
    cout << "Parsing timing..." << endl ;  
    test_timing_parser (argv[2]) ;
    cout << "Finished parsing timing." << endl ;
  }
  else if (argv[1] == string("ceff")) {
    cout << "Parsing ceff..." << endl ;  
    test_ceff_parser (argv[2]) ;
    cout << "Finished parsing ceff." << endl ;
  }
  else if (argv[1] == string("lib")) {
    cout << "Parsing lib..." << endl ;
    test_lib_parser (argv[2]) ;
    cout << "Finished parsing lib." << endl ;
  }
  else {
    cout << "Unknown format " << argv[1] << endl ;

  }
  
  return 0 ;
}
#endif
