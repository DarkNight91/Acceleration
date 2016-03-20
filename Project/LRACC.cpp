#define _CRT_SECURE_NO_DEPRECATE
#include <iostream>
#include <string>
#include <cstdio>
#include "LRACC.h"
#include "Celllib.h"
#include "Circuit.h"
#include "parser_helper.h"
#include "werc_102_vlsi.h"
#include "ReadWriteFile.h"
#define MAX_GROUP 400
using namespace std;
using namespace boost;
int k_mux;
GMAP_INDEX gGate_Index;
//map the gate name in .bench to that used in .lib file
IJRR_NAME_MAP gGate_map[] =
{
    { "not", "in" },
    { "buf", "in" },
    { "nand", "na" },
    { "NAND", "na" },
    { "and", "na" },
    { "nor", "no" },
    { "or", "no" },
    { "xor", "oa" }
};


/************************************************************************/
/*         Initialize the private variables of LRACC struct             */
/************************************************************************/
bool    LRACC::init_lracc_database(CHAR* fpart_name)
{
    string filename(fpart_name);
    filename += ".curve";
    dest_file.open(filename.c_str());

    trade_off_r.initial(0.05);
    iter_lagrangian = 0;
    lm_a = 0.1;

    gGate_Index["in01"] = INV;
    gGate_Index["na02"] = NAND2;
    gGate_Index["na03"] = NAND3;
    gGate_Index["na04"] = NAND4;
    gGate_Index["na05"] = NAND5;
    gGate_Index["na06"] = NAND6;
    gGate_Index["na07"] = NAND7;
    gGate_Index["na08"] = NAND8;
    gGate_Index["na09"] = NAND9;
    gGate_Index["no02"] = NOR2;
    gGate_Index["no03"] = NOR3;
    gGate_Index["no04"] = NOR4;
    gGate_Index["no05"] = NOR5;
    gGate_Index["no06"] = NOR6;
    gGate_Index["no07"] = NOR7;
    gGate_Index["no08"] = NOR8;
    gGate_Index["no09"] = NOR9;
#ifdef IS_ISCAS85_TEST
    gGate_Index["oa22"] = XOR2;
#else
    gGate_Index["oa12"] = OA12;
    gGate_Index["oa22"] = OA22;
    gGate_Index["ao12"] = AO12;
    gGate_Index["ao22"] = AO22;
    gGate_Index["ms00"] = FFMS;
#endif

    return  TRUE;
}

/************************************************************************/
/*   To parse the net-list info, might later also parse the tile file   */
/************************************************************************/
#ifdef IS_ISCAS85_TEST// to parse the ISCAS85 netlist
bool    LRACC::load_graph_and_tiles(CHAR *filename)
{
    UINT16DDD  num = 0;
    UINT16DDD  pinput;
    string  name;
    string  line;
    GATE    *pGate;
    TPLGY_VDES  v;
    TPLGY_EDES  inter_edge;
    map<string, TPLGY_VDES> input_vertex;       //used for parsing the net-list file
    map<string, TPLGY_VDES> middle_vertex;      //used for parsing the net-list file
    map<string, TPLGY_VDES> output_vertex;      //used for parsing the net-list file

    ifstream myfile(filename);//string  filename("c17.bench");
    t_mmlparser *mp = new t_mmlparser; 
    map<string, TPLGY_VDES>::iterator source_iter;

    if (myfile.is_open())
    {
        while (getline(myfile,line))
        {
            mp->mp_parse_clear();//vector needs to be clear since I use vector::push_back
            mp->mp_parse_topo(line.c_str());
            if (0 == mp->param.size())
            {continue;}

            if (0 == mp->mp_get_cmd_name().compare("INPUT"))
            {
                pGate = new GATE;
                pGate->name = mp->mp_get_param_name(0);
                pGate->type = IN;
                BIT_SET(pGate->flag, INPUT_DRIVER); 
                v = add_vertex(circuit);
                circuit[v].info = pGate;
                input_vertex[pGate->name] = v;
            }
            else if (0 == mp->mp_get_cmd_name().compare("OUTPUT"))
            {
                pGate = new GATE;
                pGate->name = mp->mp_get_param_name(0);
                pGate->type = OUT;
                BIT_SET(pGate->flag, OUTPUT_LOAD);
                v = add_vertex(circuit);
                circuit[v].info = pGate;
                output_vertex[pGate->name] = v;
            }
            else
            {
                pGate = new GATE;
                pGate->name = mp->mp_get_cmd_name();//mp->mp_get_param_name(0);
                v = add_vertex(circuit);
                circuit[v].info = pGate;
                middle_vertex[mp->mp_get_cmd_name()] = v;

                pinput = 0;
                for (num = 1; num <= mp->param.size()-1; num++)
                {
                    name = mp->mp_get_param_name(num);
                    source_iter = input_vertex.find(name);
                    if (source_iter != input_vertex.end())
                    {
                        //first search the input vertex
                        inter_edge = add_edge(source_iter->second, v, circuit).first;
                        pinput++;
                        circuit[inter_edge].name = name;
                        BIT_SET(circuit[inter_edge].edge_flag, INPUT_EDGE);
                    }
                    else
                    {
                        source_iter = middle_vertex.find(name);
                        H102_ASSERT(source_iter != middle_vertex.end());
                        inter_edge = add_edge(source_iter->second, v, circuit).first;
                        circuit[inter_edge].name = name;
                    }
                }
                if (pinput == mp->param.size()-1)
                {
                    BIT_SET(pGate->flag, PRIMA_INGATE);
                }

                /************************************************************************/
                /*  To get the CELLIB index for a certain type gate                     */
                /************************************************************************/
                string  footprint;
                CHAR    inpin_num[MAX_GATE_NAME];
                for (num = 0; num < MAX_CP_TYPE; ++num)
                {
                    if (0 == mp->mp_get_param_name(0).compare(gGate_map[num].netlist_name))
                    {
                        footprint = gGate_map[num].stdlib_name;
                        if (0 == footprint.compare("oa"))
                        {
                            footprint.append("22");
                        }
                        else
                        {
                            sprintf(inpin_num, "%02u", mp->param.size()-1);
                            footprint.append(inpin_num);
                        }
                        break;
                    }
                }
                H102_ASSERT(MAX_CP_TYPE != num);
                map<string,COMPONENTTYPE>::iterator   gate_lut_iter = gGate_Index.find(footprint);
                H102_ASSERT(gate_lut_iter != gGate_Index.end());
                pGate->type = gate_lut_iter->second;
            }
        }
        myfile.close();

        /************************************************************************/
        /*             to connect the primary output edges                      */
        /************************************************************************/
        map<string,TPLGY_VDES>::iterator target_iter = output_vertex.begin();
        for (; target_iter != output_vertex.end(); ++target_iter)
        {
            source_iter = middle_vertex.find(target_iter->first);
            H102_ASSERT(source_iter != middle_vertex.end());
            inter_edge = add_edge(source_iter->second, target_iter->second, circuit).first;
            circuit[inter_edge].name = target_iter->first;
            BIT_SET(circuit[inter_edge].edge_flag, OUTPUT_EDGE);
        }
        cout << "Success! Topology file "<<filename<<" parsed!"<<endl;
    }
    else cout << "Failed to open file "<<filename<<endl;

    delete mp;
    return TRUE;
}
#else

#if 0
bool    LRACC::load_graph_and_tiles(string  &filename)
{
    UINT16DDD  num = 0;
    //UINT16DDD  pinput;
    GATE    *pGate;
    TPLGY_VDES  v;
    TPLGY_EDES  inter_edge;
    map<string, TPLGY_VDES> source_of_wire;       //used for parsing the net-list file
    map<string, TPLGY_VDES> newsource_wire;
    //map<string, TPLGY_VDES> nonsource_wire;     //used for parsing the net-list file
    vector<pair<string, TPLGY_VDES> >  nonsource_wire;
    map<string, TPLGY_VDES> target_of_wire;      //used for parsing the net-list file
    map<string, TPLGY_VDES>::iterator source_iter;

    GATEOPT         *pGate_opt;

    VerilogParser vp (filename) ;

    string moduleName ;
    bool valid = vp.read_module (moduleName) ;
    assert (valid) ;

    cout << "Module " << moduleName << endl << endl ;

    do {
        string primaryInput ;
        valid = vp.read_primary_input (primaryInput) ;

        if (valid) {
            cout << "Primary input: " << primaryInput << endl ;

        // skip the ispd_clk signal, don't create Vertex for it
        if (primaryInput == "ispd_clk")
            continue;

        pGate = new GATE;
        pGate->name = primaryInput;
        pGate->cellInst = "INPUT";
        pGate->type = IN;
        BIT_SET(pGate->flag, INPUT_DRIVER);
        v = add_vertex(circuit);
        circuit[v].info = pGate;
        source_of_wire[pGate->name] = v;

        /* config the solution for input
        pGate_opt = gateopt_cast(new GATEOPT);

        pGate_opt->set_solution(0);//for input driver the solution is the first
        pGate_opt->cap = 0;
        pGate_opt->resist.resize(ADAPTIVITY_1+1);//how many resistance are there???????????
        pGate_opt->resist[ADAPTIVITY_0] = 0;//set the input driver.push_back(
        pGate_opt->resist[ADAPTIVITY_1] = 0;
        pGate_opt->offset = 0;
        pGate_opt->area = 0;
        pGate_opt->off_power.resize(ADAPTIVITY_1+1);
        pGate_opt->off_power[ADAPTIVITY_0] = 0;//set the off_power of input driver
        pGate_opt->off_power[ADAPTIVITY_1] = 0;

        COMPONENT_SOLUTION  solution;
        solution.opt = pGate_opt;

        solution.adpt_num = ADAPTIVITY_0;
        circuit[v].candidates[ADAPTIVITY_0].push_back(solution);
        circuit[v].optimal_iter[ADAPTIVITY_0] = circuit[v].candidates[ADAPTIVITY_0].begin();
#if IS_SIMULTANEOUS
        solution.adpt_num = ADAPTIVITY_1;
        circuit[v].candidates[ADAPTIVITY_1].push_back(solution);
        circuit[v].optimal_iter[ADAPTIVITY_1] = circuit[v].candidates[ADAPTIVITY_1].begin();
#endif
        //circuit[*vg.first].foward_at = primary_iter->second[0];//set the primary input time*/
        }
    } while (valid) ;

    cout << endl ;

    do {
        string primaryOutput ;
        valid = vp.read_primary_output (primaryOutput) ;

        if (valid)
        {    cout << "Primary output: " << primaryOutput << endl ;

        pGate = new GATE;
        pGate->name = primaryOutput;
        pGate->cellInst = "OUTPUT";
        pGate->type = OUT;
        BIT_SET(pGate->flag, OUTPUT_LOAD);
        v = add_vertex(circuit);
        circuit[v].info = pGate;
        target_of_wire[pGate->name] = v;

        /* config the solution for output
        pGate_opt = gateopt_cast(new GATEOPT);

        pGate_opt->set_solution(0);//for output load, the solution is the first
        pGate_opt->cap = 3;//According to Rohit, he set the output load to 3
        pGate_opt->resist.resize(ADAPTIVITY_1+1);
        pGate_opt->resist[ADAPTIVITY_0] = 0;//set the input driver.push_back(
        pGate_opt->resist[ADAPTIVITY_1] = 0;
        pGate_opt->offset = 0;
        pGate_opt->area = 0;

        COMPONENT_SOLUTION  solution;
        solution.opt = pGate_opt;
        solution.propagates.c = pGate_opt->cap;
        //solution.propagates.raq = primary_iter->second[0];//set the primary output time

        solution.adpt_num = ADAPTIVITY_0;
        circuit[v].candidates[ADAPTIVITY_0].push_back(solution);
        circuit[v].optimal_iter[ADAPTIVITY_0] = circuit[v].candidates[ADAPTIVITY_0].begin();
#if IS_SIMULTANEOUS
        solution.adpt_num = ADAPTIVITY_1;
        circuit[v].candidates[ADAPTIVITY_1].push_back(solution);
        circuit[v].optimal_iter[ADAPTIVITY_1] = circuit[v].candidates[ADAPTIVITY_1].begin();
#endif
        // circuit[*vg.first].optimal_iter = 0;circuit[*vg.first].candidates[ADAPTIVITY_0].begin();*/
        }
    } while (valid) ;

    cout << endl ;

    do {
        string net ;
        valid = vp.read_wire (net) ;

        //if (valid)
        //    cout << "Net: " << net << endl ;

    } while (valid) ;


    cout << endl ;
    cout << "Cell insts: " << std::endl ;

    do {
        string cellType, cellInst ;
        vector<std::pair<string, string> > pinNetPairs ;

        valid = vp.read_cell_inst (cellType, cellInst, pinNetPairs) ;

        if (valid) {
            cout << cellType << " " << cellInst << " " << pinNetPairs.back().second << " ";

            pGate = new GATE;
            //to make sure Rohit use the output name as the cell's name
            pGate->name = pinNetPairs.back().second;//or cellInst instance's name?
            pGate->cellInst = cellInst;
            v = add_vertex(circuit);
            circuit[v].info = pGate;

            H102_ASSERT("o" == pinNetPairs.back().first);//make sure this is the output
            newsource_wire[pinNetPairs.back().second] = v;//the output wirename is the key to find the source/target Vertex

            //pinput = 0;
            for (int i=0; i < pinNetPairs.size()-1; ++i)//don't consider the output pin anymore
            {
                cout << "(" << pinNetPairs[i].first << " " << pinNetPairs[i].second << ") " ;
                if ("ck" == pinNetPairs[i].first)
                {
                    H102_ASSERT("ispd_clk" == pinNetPairs[i].second);
                    continue;
                }

                source_iter = source_of_wire.find(pinNetPairs[i].second);
                if (source_iter != source_of_wire.end())
                {
                    inter_edge = add_edge(source_iter->second, v, circuit).first;
                    //pinput++;
                    circuit[inter_edge].name = pinNetPairs[i].second;
                    BIT_SET(circuit[inter_edge].edge_flag, INPUT_EDGE);
                } 
                else
                {
                    nonsource_wire.push_back(make_pair(pinNetPairs[i].second, v));//[pinNetPairs[i].second] = v;

                    /*source_iter = target_of_wire.find(pinNetPairs[i].second);
                    if (source_iter != target_of_wire.end())
                    {
                        inter_edge = add_edge(source_iter->second, v, circuit).first;
                        circuit[inter_edge].name = pinNetPairs[i].second;
                        //BIT_SET(circuit[inter_edge].edge_flag, OUTPUT_EDGE);
                    } 
                    else
                    {
                        nonsource_wire.push_back(make_pair(pinNetPairs[i].second, v));//[pinNetPairs[i].second] = v;
                    }*/
                }
            }
            //if (pinput == pinNetPairs.size()-1)
            //{
            //    BIT_SET(pGate->flag, PRIMA_INGATE);
            //}

            // set the type ID for the cell
            string  footprint;
            footprint.assign(cellType, 0, MAX_FOOTPRINT_LEN);//just get the first 4 char as the footprint

            map<string,COMPONENTTYPE>::iterator   gate_lut_iter = gGate_Index.find(footprint);
            H102_ASSERT(gate_lut_iter != gGate_Index.end());
            pGate->type = gate_lut_iter->second;

            cout << endl ;
        }
    } while (valid) ;

    // to connect the nonsource edge during the parsing process
    for (num = 0; num < nonsource_wire.size(); ++num)
    {
        source_iter = newsource_wire.find(nonsource_wire[num].first);
        H102_ASSERT(source_iter != newsource_wire.end());
        inter_edge = add_edge(source_iter->second, nonsource_wire[num].second, circuit).first;
        circuit[inter_edge].name = source_iter->first;
    }
    /*target_iter = nonsource_wire.begin();
    for (; target_iter != nonsource_wire.end(); ++target_iter)
    {
        source_iter = newsource_wire.find(target_iter->first);
        H102_ASSERT(source_iter != newsource_wire.end());
        inter_edge = add_edge(source_iter->second, target_iter->second, circuit).first;
        circuit[inter_edge].name = target_iter->first;
    }*/

    map<string,TPLGY_VDES>::iterator target_iter;
    // to connect the primary output edge
    target_iter = target_of_wire.begin();
    for (; target_iter != target_of_wire.end(); ++target_iter)
    {
        source_iter = newsource_wire.find(target_iter->first);
        H102_ASSERT(source_iter != newsource_wire.end());
        inter_edge = add_edge(source_iter->second, target_iter->second, circuit).first;
        circuit[inter_edge].name = target_iter->first;
        BIT_SET(circuit[inter_edge].edge_flag, OUTPUT_EDGE);
    }

    cout << endl ;
    cout << "Success! Topology file "<<filename<<" parsed!"<<endl;
}
#endif

bool    LRACC::load_graph_and_tiles(string  &filename)
{
    UINT16D  num = 0;
    //UINT16DDD  pinput;
    GATE    *pGate;
    TPLGY_VDES  v;
    TPLGY_EDES  inter_edge;
    map<string, TPLGY_VDES> source_of_wire;       //used for parsing the net-list file
    map<string, TPLGY_VDES> newsource_wire;
    //map<string, TPLGY_VDES> nonsource_wire;     //used for parsing the net-list file
    vector<pair<string, TPLGY_VDES> >  nonsource_wire;
    map<string, TPLGY_VDES> target_of_wire;      //used for parsing the net-list file
    map<string, TPLGY_VDES>::iterator source_iter;

    GATEOPT         *pGate_opt;
    COMPONENTTYPE   type;
    {
        VerilogParser vp (filename) ;

        string moduleName ;
        bool valid = vp.read_module (moduleName) ;
        assert (valid) ;

        cout << "Module " << moduleName << endl << endl ;

        do {
            string primaryInput ;
            valid = vp.read_primary_input (primaryInput) ;

            if (valid) 
            {   //cout << "Primary input: " << primaryInput << endl ;

                // skip the ispd_clk signal, don't create Vertex for it
                if (primaryInput == "ispd_clk")
                    continue;

                pGate = new GATE;
                pGate->name = primaryInput;
                pGate->cellInst = "INPUT";
                pGate->type = IN1;
                pGate->tile_idx = MAX_TUNING_FFF;
                BIT_SET(pGate->flag, INPUT_DRIVER);
                v = add_vertex(circuit);
                circuit[v].info = pGate;
                source_of_wire[pGate->name] = v;
            }
        } while (valid) ;

        cout << endl ;

        do {
            string primaryOutput ;
            valid = vp.read_primary_output (primaryOutput) ;

            if (valid)
            {   // cout << "Primary output: " << primaryOutput << endl ;

                pGate = new GATE;
                pGate->name = primaryOutput;
                pGate->cellInst = "OUTPUT";
                pGate->type = OUT1;
                pGate->tile_idx = MAX_TUNING_FFF;
                BIT_SET(pGate->flag, OUTPUT_LOAD);
                v = add_vertex(circuit);
                circuit[v].info = pGate;
                target_of_wire[pGate->name] = v;
            }
        } while (valid) ;

        cout << endl ;

        do {
            string net ;
            valid = vp.read_wire (net) ;

            //if (valid)
            //    cout << "Net: " << net << endl ;

        } while (valid) ;


        cout << endl ;
        cout << "Cell insts: " << std::endl ;

        do {
            string cellType, cellInst ;
            vector<std::pair<string, string> > pinNetPairs ;

            valid = vp.read_cell_inst (cellType, cellInst, pinNetPairs) ;

            if (valid) {
                //cout << cellType << " " << cellInst << " " << pinNetPairs.back().second << " ";

                // set the type ID for the cell
                string  footprint;
                footprint.assign(cellType, 0, MAX_FOOTPRINT_LEN);//just get the first 4 char as the footprint

                map<string,COMPONENTTYPE>::iterator   gate_lut_iter = gGate_Index.find(footprint);
                H102_ASSERT(gate_lut_iter != gGate_Index.end());
                type = gate_lut_iter->second;

                if (FFMS != type)
                {
                    pGate = new GATE;
                    //to make sure Rohit use the output name as the cell's name
                    pGate->name = pinNetPairs.back().second;//or cellInst instance's name?
                    pGate->cellInst = cellInst;
                    pGate->type = type;
                    v = add_vertex(circuit);
                    circuit[v].info = pGate;

                    H102_ASSERT("o" == pinNetPairs.back().first);//make sure this is the output
                    newsource_wire[pinNetPairs.back().second] = v;//the output wirename is the key to find the source/target Vertex
                }
                else
                {
                    //treat the input of FF as PO and output of FF as PI
                    pGate = new GATE;
                    H102_ASSERT("o" == pinNetPairs.back().first);//make sure this is the output o
                    pGate->name = pinNetPairs.back().second;
                    pGate->cellInst = cellInst+"_o";
                    pGate->type = type;
                    pGate->tile_idx = MAX_TUNING_FFF;
                    BIT_SET(pGate->flag, FF_VIRTUAL);//don't set INPUT_DRIVER; treat it as a normal gate
                    v = add_vertex(circuit);
                    circuit[v].info = pGate;
                    source_of_wire[pGate->name] = v;

                    pGate = new GATE;
                    H102_ASSERT("d" == pinNetPairs[1].first);//make sure this is the input d
                    pGate->name = pinNetPairs[1].second;
                    pGate->cellInst = cellInst+"_d";//"OUTPUT";
                    pGate->type = type;
                    pGate->tile_idx = MAX_TUNING_FFF;
                    BIT_SET(pGate->flag, OUTPUT_LOAD|FF_VIRTUAL);
                    v = add_vertex(circuit);
                    circuit[v].info = pGate;
                    target_of_wire[pGate->name] = v;
                }

                //cout << endl ;
            }
        } while (valid) ;
    }

    {
        VerilogParser vp (filename) ;

        string moduleName ;
        bool valid = vp.read_module (moduleName) ;
        assert (valid) ;

        cout << "Module " << moduleName << endl << endl ;

        do {
            string primaryInput ;
            valid = vp.read_primary_input (primaryInput) ;

            //if (valid)
            //    cout << "Primary input: " << primaryInput << endl ;
        } while (valid) ;

        cout << endl ;

        do {
            string primaryOutput ;
            valid = vp.read_primary_output (primaryOutput) ;

            //if (valid)
            //    cout << "Primary output: " << primaryOutput << endl ;
        } while (valid) ;

        cout << endl ;

        do {
            string net ;
            valid = vp.read_wire (net) ;

            //if (valid)
            //    cout << "Net: " << net << endl ;

        } while (valid) ;


        cout << endl ;
        cout << "Cell insts: " << std::endl ;

        do {
            string cellType, cellInst ;
            vector<std::pair<string, string> > pinNetPairs ;

            valid = vp.read_cell_inst (cellType, cellInst, pinNetPairs) ;

            if (valid) {
                //cout << cellType << " " << cellInst << " " << pinNetPairs.back().second << " ";
                string  footprint;
                footprint.assign(cellType, 0, MAX_FOOTPRINT_LEN);//just get the first 4 char as the footprint

                map<string,COMPONENTTYPE>::iterator   gate_lut_iter = gGate_Index.find(footprint);
                H102_ASSERT(gate_lut_iter != gGate_Index.end());
                type = gate_lut_iter->second;

                if (FFMS != type)
                {
                    TPLGY_VDES  tNode;
                    tNode = newsource_wire[pinNetPairs.back().second];//first find self
                    for (int i=0; i < pinNetPairs.size()-1; ++i)//don't consider the output pin anymore
                    {
                        //cout << "(" << pinNetPairs[i].first << " " << pinNetPairs[i].second << ") " ;
                        if ("ck" == pinNetPairs[i].first)
                        {
                            H102_ASSERT("ispd_clk" == pinNetPairs[i].second);
                            continue;
                        }

                        source_iter = source_of_wire.find(pinNetPairs[i].second);
                        if (source_iter != source_of_wire.end())
                        {
                            inter_edge = add_edge(source_iter->second, tNode, circuit).first;
                            //pinput++;
                            circuit[inter_edge].name = pinNetPairs[i].second;
                            BIT_SET(circuit[inter_edge].edge_flag, INPUT_EDGE);
                        }
                        else
                        {
                            source_iter = newsource_wire.find(pinNetPairs[i].second);
                            H102_ASSERT(source_iter != newsource_wire.end());
                            inter_edge = add_edge(source_iter->second, tNode, circuit).first;
                            circuit[inter_edge].name = pinNetPairs[i].second;
                        }
                    }
                } 
                //if (pinput == pinNetPairs.size()-1)
                //{
                //    BIT_SET(pGate->flag, PRIMA_INGATE);
                //}
                //cout << endl ;
            }
        } while (valid) ;

        map<string,TPLGY_VDES>::iterator target_iter;
        // to connect the primary output node
        target_iter = target_of_wire.begin();
        for (; target_iter != target_of_wire.end(); ++target_iter)
        {
            source_iter = newsource_wire.find(target_iter->first);
            if (source_iter != newsource_wire.end())
            {
                inter_edge = add_edge(source_iter->second, target_iter->second, circuit).first;
                circuit[inter_edge].name = target_iter->first;
                BIT_SET(circuit[inter_edge].edge_flag, OUTPUT_EDGE);
            }
            else
            {
                source_iter = source_of_wire.find(target_iter->first);
                H102_ASSERT(source_iter != source_of_wire.end());
                inter_edge = add_edge(source_iter->second, target_iter->second, circuit).first;
                circuit[inter_edge].name = target_iter->first;
                BIT_SET(circuit[inter_edge].edge_flag, OUTPUT_EDGE);
            }
        }

        /*/ to connect the new created input node?????
        target_iter = target_of_wire.begin();
        for (; target_iter != target_of_wire.end(); ++target_iter)
        {
            
            //it's possible there is no connection between the input;
            if (source_iter != source_of_wire.end())
            {
                H102_ASSERT(FFMS == circuit[source_iter->second].info->type);
                inter_edge = add_edge(source_iter->second, target_iter->second, circuit).first;
                circuit[inter_edge].name = target_iter->first;
                BIT_SET(circuit[inter_edge].edge_flag, OUTPUT_EDGE);
            }
        }*/
    }

    cout << endl ;
    cout << "Success! Topology file "<<filename<<" parsed!"<<endl;
	return 1;
}
#endif


/************************************************************************/
/*          To parse the cell information from contest.lib              */
/************************************************************************/
bool LRACC::load_celllib(CHAR *filename)
{
    FLT32           bb_tmp;
    FLT32           delay_eql;
    string          old_name;//old footprint
    UINT8D           cell_pin;
    UINT8D           apt_index;
    UINT8D           sol_Cnt = 0;
    INT8D            tune_shift;
    UINT8D           size1;
    LibParserLUT    *pTimingarc;

    LibParser lp (filename);
    double maxTransition = 0.0 ;
    bool valid = lp.read_default_max_transition(maxTransition) ;
    H102_ASSERT(valid) ;
    cout << "The default max transition defined is " << maxTransition << endl ;
    cell_lib.resize(MAX_CP_TYPE);

    int readCnt = 0 ;
    do
    {
        LibParserCellInfo cell ;
        valid = lp.read_cell_info (cell) ;
        if (valid)
        {
            ++readCnt ;
            //cout << cell << endl ;//ostream& operator<< (ostream& os, LibParserCellInfo& cell) {

            map<string,COMPONENTTYPE>::iterator   gate_lut_iter = gGate_Index.find(cell.footprint);
            if (gate_lut_iter == gGate_Index.end())
            {
                //cout<<"Unknown type: "<<cell.footprint << " skip" << endl;
                continue;//ignore unknown type in contest.lib, e.g. ms00, vss, vcc
            }

            //cell already in local data lib
            if (old_name.compare(cell.footprint))
            {old_name = cell.footprint;sol_Cnt = 0;}
            else if ((++sol_Cnt%10) >= 7)//(++sol_Cnt >= MAX_SOL_VALID_NUM)
            {continue;}
            
            // simply choose the first pin assuming that all pin have the same input cap, pins[0] is o
            cell_pin = (FFMS == gate_lut_iter->second)?2:1;
            {
                GATEOPT *pGate_opt = gateopt_cast(new GATEOPT);
                pGate_opt->area = cell.area;
                //config the Vt, gate size(width and length) to the option
                pGate_opt->set_solution(sol_Cnt);

                /************************************************************************/
                /* set input cap, output resist and delay with no adaptivity               */
                /************************************************************************/
                pGate_opt->cap = cell.pins[cell_pin].capacitance;
                H102_ASSERT("o" != cell.pins[cell_pin].name);
                // simply choose the first timing arc thinking that the timing arc are the same and 
                // simply choose the fallDelay LUT without other computations. e.g cell+transition
                pTimingarc = &(cell.timingArcs[0].fallDelay);
                // simply choose the third column's first entry as the offset of the arc delay
                size1 = pTimingarc->tableVals.size()-1;
                pGate_opt->offset = 0;//pTimingarc->tableVals[0][INPUT_SLEW_INDEX];//offset[0]
                delay_eql = pTimingarc->tableVals[size1][INPUT_SLEW_INDEX] - pGate_opt->offset;

                H102_ASSERT(0 == pTimingarc->loadIndices[0]);
                // simply treat the delay and slew LUT is linear so we get equivalent resistance 
                pGate_opt->resist.push_back(delay_eql/pTimingarc->loadIndices[size1]);//resist[0]
                pGate_opt->off_power.push_back(cell.leakagePower);//off_power[0]
                H102_ASSERT((FFMS == gate_lut_iter->second) || cell.leakagePower);

                /************************************************************************/
                /*  generate the cap, resist and offset based on existing configuration  */
                /************************************************************************/
                for (apt_index = FBB_GRADE_1; apt_index < MAX_BB_GRADE; ++apt_index)
                {
                    GET_SHIFT_UNIT(apt_index, tune_shift);
                    bb_tmp = pGate_opt->get_leakpower_bb(tune_shift);
                    pGate_opt->off_power.push_back(bb_tmp);
                    bb_tmp = pGate_opt->get_resist_bb(tune_shift);
                    pGate_opt->resist.push_back(bb_tmp);
                }

                cell_lib[gate_lut_iter->second].push_back(pGate_opt);
            }
        }
    } while (valid) ;

    cout << "Read " << readCnt << " number of library cells" << endl ;
    //store the order and don't need to topological search again
    topological_sort(circuit, back_inserter(order));
    return TRUE;
}


/************************************************************************/
/*   To parse supplementary info for cells which are not found in contest.lib     */
/************************************************************************/
bool LRACC::load_supplement(CHAR *filename)
{
    UINT8D   index;
    UINT8D   num;
    UINT8D   apt_index;
    ifstream is(filename);//"supplement2.lib");
    vector<string> tokens ;
    BOOLEAN valid = read_line_as_tokens(is, tokens) ;
    LIB_LUT gate_lut;

    while (valid)
    {
        if (tokens.size() == 7 && tokens[1] == "name")
        {
            //tokens[2]; is the cell.footprint
            begin_read_cell_param (is, gate_lut);

            map<string,COMPONENTTYPE>::iterator   gate_lut_iter = gGate_Index.find(tokens[2]);
            if (gate_lut_iter == gGate_Index.end())
            {
                cout<<"Unknown supplement type: "<<tokens[2] << endl;
                //ignore unknown type in contest.lib, e.g. 
            }
            else//cell already in local data lib
            {
                for (num = 0; num < PN_SIZE_MAX; ++num)
                {
                    for (index = 0; index < VT_LEVEL_MAX; ++index)
                    {
                        GATEOPT *pGate_opt = gateopt_cast(new GATEOPT);

                        pGate_opt->set_solution(index*PN_SIZE_MAX + num);
                        pGate_opt->area = gate_lut[index][num].area;
                        pGate_opt->cap = gate_lut[index][num].cap;
                        pGate_opt->resist.push_back(gate_lut[index][num].resist);
                        pGate_opt->offset = gate_lut[index][num].offset;
                        pGate_opt->off_power.push_back(gate_lut[index][num].leakpower);

                        for (apt_index = ADAPTIVITY_0;; ++apt_index)
                        {
                            //waiting for Rohit to make this better
                            //pGate_opt->cap.push_back(pGate_opt->cap[0]);//cap[index]
                            pGate_opt->resist.push_back(pGate_opt->resist[0]);
                            //pGate_opt->offset.push_back(pGate_opt->offset[0]);
                            pGate_opt->off_power.push_back(pGate_opt->off_power[0]);
                            if (ADAPTIVITY_1 == apt_index)
                                break;
                        }

                        cell_lib[gate_lut_iter->second].push_back(pGate_opt);
                    }
                }
            }
        }
        valid = read_line_as_tokens (is, tokens) ;
    }

    return TRUE;
}


/************************************************************************/
/*  To parse the primary input timing/driver and output timing/load info */
/************************************************************************/
#ifdef IS_ISCAS85_TEST
bool LRACC::load_at_raq(CHAR *filename)
{
    //UINT8DD           num = 0;
    string          value;
    string          line;

    typedef map<string, vector<FLT32> > STR_VFLT32;
    STR_VFLT32      primary_at;//used for parsing the primary info, the vector<FLT32>
    STR_VFLT32      primary_rqt;// the first FLT32 is at or rat, the second FLT32 is R or C
    STR_VFLT32      *primary;

    ifstream myfile(filename);
    t_mmlparser     *mp = new t_mmlparser;
    vector<FLT32>   tmp;
    GATEOPT         *pGate_opt;
    //EDGE_SOL        edge_sol;

    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            mp->mp_parse_clear();//vector needs to be clear since I use vector::push_back
            mp->mp_parse_topo(line.c_str());
            tmp.clear();

            if (0 == mp->mp_get_cmd_name().compare("INPUT"))
            {
                primary = &primary_at;
            }
            else
            {
                H102_ASSERT(0 == mp->mp_get_cmd_name().compare("OUTPUT"));
                primary = &primary_rqt;
            }
            value = mp->mp_get_param_name(1);//get the at||rat
            tmp.push_back((FLT32)atof(value.c_str()));
            value = mp->mp_get_param_name(2);//get the R||C
            tmp.push_back((FLT32)atof(value.c_str()));
            (*primary)[mp->mp_get_param_name(0)] = tmp;
        }
        myfile.close();

        /* to get the R of input driver and the primary input time of the edge */
        pair<TPLGY_VITER, TPLGY_VITER> vg = vertices(circuit);
        STR_VFLT32::iterator primary_iter;
        for (; vg.first != vg.second; ++vg.first)
        {
            if (BIT_TEST(circuit[*vg.first].info->flag, INPUT_DRIVER))
            {
                primary_iter = primary_at.find(circuit[*vg.first].info->name);
                H102_ASSERT(primary_iter != primary_at.end());

                pGate_opt = gateopt_cast(new GATEOPT);

                pGate_opt->set_solution(0);//for input driver the solution is the first
                pGate_opt->cap = 0;
                pGate_opt->resist.resize(ADAPTIVITY_1+1);//how many resistance are there???????????
                pGate_opt->resist[ADAPTIVITY_0] = primary_iter->second[1];//set the input driver.push_back(
                pGate_opt->resist[ADAPTIVITY_1] = primary_iter->second[1];
                pGate_opt->offset = 0;
                pGate_opt->area = 0;
                pGate_opt->off_power.resize(ADAPTIVITY_1+1);
                pGate_opt->off_power[ADAPTIVITY_0] = 0;//set the off_power of input driver
                pGate_opt->off_power[ADAPTIVITY_1] = 0;

                COMPONENT_SOLUTION  solution;
                solution.opt = pGate_opt;

                solution.adpt_num = ADAPTIVITY_0;
                circuit[*vg.first].candidates[ADAPTIVITY_0].push_back(solution);
                circuit[*vg.first].optimal_iter[ADAPTIVITY_0] = circuit[*vg.first].candidates[ADAPTIVITY_0].begin();
#if IS_SIMULTANEOUS
                solution.adpt_num = ADAPTIVITY_1;
                circuit[*vg.first].candidates[ADAPTIVITY_1].push_back(solution);
                circuit[*vg.first].optimal_iter[ADAPTIVITY_1] = circuit[*vg.first].candidates[ADAPTIVITY_1].begin();
#endif
                //circuit[*vg.first].foward_at = primary_iter->second[0];//set the primary input time
            }
            else if (BIT_TEST(circuit[*vg.first].info->flag, OUTPUT_LOAD))
            {
                primary_iter = primary_rqt.find(circuit[*vg.first].info->name);
                H102_ASSERT(primary_iter != primary_rqt.end());

                pGate_opt = gateopt_cast(new GATEOPT);

                pGate_opt->set_solution(0);//for output load, the solution is the first
                pGate_opt->cap = primary_iter->second[1];//set the output load
                pGate_opt->resist.push_back(0);
                pGate_opt->offset = 0;
                pGate_opt->area = 0;

                COMPONENT_SOLUTION  solution;
                solution.opt = pGate_opt;
                solution.propagates.c = pGate_opt->cap;
                //solution.propagates.raq = primary_iter->second[0];//set the primary output time

                solution.adpt_num = ADAPTIVITY_0;
                circuit[*vg.first].candidates[ADAPTIVITY_0].push_back(solution);
                circuit[*vg.first].optimal_iter[ADAPTIVITY_0] = circuit[*vg.first].candidates[ADAPTIVITY_0].begin();
#if IS_SIMULTANEOUS
                solution.adpt_num = ADAPTIVITY_1;
                circuit[*vg.first].candidates[ADAPTIVITY_1].push_back(solution);
                circuit[*vg.first].optimal_iter[ADAPTIVITY_1] = circuit[*vg.first].candidates[ADAPTIVITY_1].begin();
#endif
                // circuit[*vg.first].optimal_iter = 0;circuit[*vg.first].candidates[ADAPTIVITY_0].begin();
                //cout<<"cap = "<<cap_resist.first.mean<<" resist = "<<cap_resist.second.mean<<endl;
            }
            else
            {
                break;//Based on the boost graph lib, the vertex is sequentially added with add_vertex
            }
        }
        cout << "Success! Primary input/output file " << filename <<" parsed!"<< endl;
    }
    else  cout << "Failed to open file " << filename << endl;

    delete mp;
    return TRUE;
}

#else
bool    LRACC::load_at_raq()
{
    GATEOPT         *pGate_opt;

    /* to get the R of input driver and the primary input time of the edge */
    pair<TPLGY_VITER, TPLGY_VITER> vg = vertices(circuit);
    for (; vg.first != vg.second; ++vg.first)
    {
        if (FFMS == circuit[*vg.first].info->type)
        {
            pGate_opt = gateopt_cast(cell_lib[FFMS][0]);
            COMPONENT_SOLUTION  solution;
            solution.opt = pGate_opt;
            solution.propagates.c = pGate_opt->cap;

            solution.adpt_num = ADAPTIVITY_0;
            circuit[*vg.first].candidates[ADAPTIVITY_0].push_back(solution);
            circuit[*vg.first].optimal_iter[ADAPTIVITY_0] = circuit[*vg.first].candidates[ADAPTIVITY_0].begin();
#if IS_SIMULTANEOUS
            solution.adpt_num = ADAPTIVITY_1;
            circuit[*vg.first].candidates[ADAPTIVITY_1].push_back(solution);
            circuit[*vg.first].optimal_iter[ADAPTIVITY_1] = circuit[*vg.first].candidates[ADAPTIVITY_1].begin();
#endif
        } 
        else if (BIT_TEST(circuit[*vg.first].info->flag, INPUT_DRIVER))
        {
            pGate_opt = gateopt_cast(new GATEOPT);

            pGate_opt->set_solution(0);//for input driver the solution is the first
            pGate_opt->cap = 0;
            pGate_opt->resist.resize(ADAPTIVITY_1+1);//how many resistance are there???????????
            pGate_opt->resist[ADAPTIVITY_0] = 0;//set the input driver.push_back(
            pGate_opt->resist[ADAPTIVITY_1] = 0;
            pGate_opt->offset = 0;
            pGate_opt->area = 0;
            pGate_opt->off_power.resize(ADAPTIVITY_1+1);
            pGate_opt->off_power[ADAPTIVITY_0] = 0;//set the off_power of input driver
            pGate_opt->off_power[ADAPTIVITY_1] = 0;

            COMPONENT_SOLUTION  solution;
            solution.opt = pGate_opt;

            solution.adpt_num = ADAPTIVITY_0;
            circuit[*vg.first].candidates[ADAPTIVITY_0].push_back(solution);
            circuit[*vg.first].optimal_iter[ADAPTIVITY_0] = circuit[*vg.first].candidates[ADAPTIVITY_0].begin();
#if IS_SIMULTANEOUS
            solution.adpt_num = ADAPTIVITY_1;
            circuit[*vg.first].candidates[ADAPTIVITY_1].push_back(solution);
            circuit[*vg.first].optimal_iter[ADAPTIVITY_1] = circuit[*vg.first].candidates[ADAPTIVITY_1].begin();
#endif
            //circuit[*vg.first].foward_at = primary_iter->second[0];//set the primary input time
        }
        else if (BIT_TEST(circuit[*vg.first].info->flag, OUTPUT_LOAD))
        {
            pGate_opt = gateopt_cast(new GATEOPT);

            pGate_opt->set_solution(0);//for output load, the solution is the first
            pGate_opt->cap = 3;//set the output load
            pGate_opt->resist.push_back(0);
            pGate_opt->offset = 0;
            pGate_opt->area = 0;

            COMPONENT_SOLUTION  solution;
            solution.opt = pGate_opt;
            solution.propagates.c = pGate_opt->cap;
            //solution.propagates.raq = primary_iter->second[0];//set the primary output time

            solution.adpt_num = ADAPTIVITY_0;
            circuit[*vg.first].candidates[ADAPTIVITY_0].push_back(solution);
            circuit[*vg.first].optimal_iter[ADAPTIVITY_0] = circuit[*vg.first].candidates[ADAPTIVITY_0].begin();
#if IS_SIMULTANEOUS
            solution.adpt_num = ADAPTIVITY_1;
            circuit[*vg.first].candidates[ADAPTIVITY_1].push_back(solution);
            circuit[*vg.first].optimal_iter[ADAPTIVITY_1] = circuit[*vg.first].candidates[ADAPTIVITY_1].begin();
#endif
            // circuit[*vg.first].optimal_iter = 0;circuit[*vg.first].candidates[ADAPTIVITY_0].begin();
            //cout<<"cap = "<<cap_resist.first.mean<<" resist = "<<cap_resist.second.mean<<endl;
        }
        else
        {
            //break;//Based on the boost graph lib, the vertex is sequentially added with add_vertex
        }
    }
    cout << "Success! Primary input/output parsed!"<< endl;
	return 1;
}
#endif

/************************************************************************/
/*      Load the tunning unit information                               */
/************************************************************************/
bool    LRACC::load_tuners()
{
#define     GRADE_LOW       (0)
#define     GRADE_HIGH      (1)
    GRADE   grad_tmp;
    TILE    *pTile;
    UINT8D   num;

    gGrid_num = no_of_adaptive_blocks;//MAX_TUNING_UNIT;//later this number might change
    H102_ASSERT(gGrid_num);
    cout<<"Total number of grid = "<<gGrid_num<<endl;
    pTile = new TILE;
    for (num = 0; num < gGrid_num; ++num)
    {
        pTile->is_adaptable = NO;
        grad_tmp.vbb = GRADE_LOW;
        pTile->tuner.allstates.push_back(grad_tmp);
        grad_tmp.vbb = GRADE_HIGH;
        pTile->tuner.allstates.push_back(grad_tmp);

        pTile->tuner.state = GRADE_LOW;
        tile.push_back(*pTile);
    }

    return  TRUE;
}


/************************************************************************/
/*      To update the gate size info for Rohit's ssta                   */
/************************************************************************/
FLT32   LRACC::update_gates_info()
{
    FLT32    area_sum = 0;
    vector<float>   param;// The map key is "Gate/Net name" and value is the "Size, R, C, Vth"
    VDES_CONTAINER::iterator ii = order.begin();
    TILE        *pTile;
    GATEOPT     *pGate_opt;
    UINT8D       adpt_num;
    GATEOPT     tmp_opt;
    tmp_opt.resist.resize(2);
    tmp_opt.off_power.resize(2);

    param.resize(SSTA_PARAM_MAX);
    //cout << "A reverse topological ordering: " << endl;
    for (; ii != order.end(); ++ii)
    {
        if (BIT_TEST(circuit[*ii].info->flag,OUTPUT_LOAD)
            ||BIT_TEST(circuit[*ii].info->flag,INPUT_DRIVER))//IN == circuit[*ii].info->type)
        {
            continue;//skip the INPUT and OUTPUT gate
        }
        else
        {
            //initialize the gate sizes for SSTA
            if (0 == iter_lagrangian)
            {
                pGate_opt = gateopt_cast(cell_lib[circuit[*ii].info->type][0]);
            } 
            //update the optimized options for SSTA
            else
            {
#if IS_SIMULTANEOUS
                pTile = &(tile[circuit[*ii].info->tile_idx]);
                adpt_num = pTile->is_adaptable;//whether there is a adaptivity circuit
                pGate_opt = gateopt_cast(circuit[*ii].optimal_iter[adpt_num]->opt);                
#else
                pGate_opt = gateopt_cast(circuit[*ii].optimal_iter[ADAPTIVITY_0]->opt);
#endif
            }

            area_sum += pGate_opt->area;
            param[SSTA_PARAM_IO1] = pGate_opt->w;
            //always give the resistance of the Low Grade
            param[SSTA_PARAM_IO2] = pGate_opt->resist[ZERO_GRADE];
            param[SSTA_PARAM_IO3] = pGate_opt->cap;
            param[SSTA_PARAM_IO4] = pGate_opt->vth;
            param[SSTA_PARAM_IO5] = pGate_opt->off_power[ZERO_GRADE];
            //H102_ASSERT(param[SSTA_PARAM_IO5]);
            gate_sizes[circuit[*ii].info->name] = param;
        }
    }


    //update the adaptivity block info
    UINT8D   num;
    for (num = 0; num < gGrid_num; ++num)
    {
        adpt_block_exist[num] = tile[num].is_adaptable;
    }
    return area_sum;
}

/************************************************************************/
/*       To run the DFS algorithm and mark the redundant edge           */
/************************************************************************/
void    LRACC::load_plca_mark_edge()
{
    //store the order and don't need to topological search again
    //topological_sort(circuit, back_inserter(order));

    //set the gate position and block info based on the PLCA I/O 
    VDES_CONTAINER::reverse_iterator ii;
    GATE    *pGate;
    for (ii = order.rbegin(); ii != order.rend(); ++ii)
    {
        pGate = gate_cast(circuit[*ii].info);
        if ((BIT_TEST(pGate->flag, INPUT_DRIVER))
			|| (BIT_TEST(pGate->flag, OUTPUT_LOAD)) || (BIT_TEST(pGate->flag, FF_VIRTUAL)))// IN == pGate->type
        {
            pGate->tile_idx = MAX_TUNING_FFF;
            continue;
        }
        else
        {
            map<string, vector<float> >::iterator   iter_plac = plac_info.find(pGate->name);
            pGate->positionX = iter_plac->second[PLAC_PARAM_IO1];
            pGate->positionY = iter_plac->second[PLAC_PARAM_IO2];
            pGate->tile_idx = (int)(iter_plac->second[PLAC_PARAM_IO3]);
        }
    }

    //set the redundant flag for certain edge
    UINT32D  num;
    vector<TPLGY_EDES> has_cycle;//store the EDGE descriptor that need to be marked
    redundant_detector vis(&has_cycle);
    dfs_print_visitor prt;
    stack_gates_in_tile finish(&tile);

    //boost::depth_first_search(circuit, visitor(prt));
    cout<<endl;
    boost::depth_first_search(circuit, visitor(vis));
    cout<<endl << "The DFS mark finished " << "Has redundant edge? " << has_cycle.size() << endl;
    boost::depth_first_search(circuit, visitor(finish));
    cout<<endl << "finished vertex print over " << endl;

    for (num = 0; num < has_cycle.size(); ++num)
    {
        BIT_SET(circuit[has_cycle[num]].edge_flag, REDUNDANT_EDGE);
    }

    //final_result.tile_adpt.resize(gGrid_num);
    //update the overhead of the adaptivity block
    for (num = 0; num < gGrid_num; ++num)
    {
        tile[num].set_tune_overhead(IS_NO_ADAPTIVITY);
        cout<<"Tile["<<num<<"] area_overhead = "<<tile[num].tuner.area_overhead;
        cout<<"\t power_overhead = "<<tile[num].tuner.power_overhead;
        cout<<"\t gate#= "<<tile[num].components.size()<<endl;
    }
}

/************************************************************************/
/*         Compute the variation of L and W at very beginning           */
/************************************************************************/
#ifdef  IS_OPTIMAL_CODE
void    LRACC::opt_gate_solution()
{
    COMPONENTTYPE   type;
    GATEOPT *pGate_opt;
    UINT16DDD  sol_num;
    FLT32   tmp;
    UINT8DD   adpt_num;
    FLT32   gate_r[2];
    UINT16DDD  num;
    UINT16DDD  middle = 0;
    GATE    *pGate;
    FLT32   delta_L = 0;
    FLT32   delta_W = 0;
    pair<FLT32,FLT32>   overhead;

    VDES_CONTAINER::iterator ii = order.begin();
    for (; ii != order.end(); ++ii)
    {
        pGate = gate_cast(circuit[*ii].info);
        if (!BIT_TEST(circuit[*ii].info->flag,INPUT_DRIVER) && !(BIT_TEST(circuit[*ii].info->flag,OUTPUT_LOAD)))
        {
            type = circuit[*ii].info->type;
            for (sol_num = 0; sol_num < cell_lib[type].size(); ++sol_num)
            {
                pGate_opt = gateopt_cast(cell_lib[type][sol_num]);
                gate_r[ADAPTIVITY_0] = pGate_opt->resist[ZERO_GRADE];
                gate_r[ADAPTIVITY_0+1] = pGate_opt->resist[FBB_GRADE_1];

            /*
            r = (ADAPTIVITY_0 == adpt_num)?
                o.resist[ZERO_GRADE]:o.resist[FBB_GRADE_1];
        
            //get the dR/dL()
            dRdL = r*0.15/3;// RESIST_PARAM_P * pow(w,-1);
            delta_L = this->info->deltaL * pow(dRdL, 2);

            //get the dR/dW()
            dRdW = r*0.08/3;//-RESIST_PARAM_P * pow(w,-2);
            delta_W = this->info->deltaW * pow(dRdW, 2);

            gate_delay += 3*sqrt(delta_L + delta_W)*c;*/

            map<string, vector<double> >::iterator  iter_pca = pca_params.find(pGate->name);
            H102_ASSERT(iter_pca != pca_params.end());

            delta_L = 0;
            delta_W = 0;
            middle = iter_pca->second.size()/2;
            for (num = 0; num < middle; ++num)
            {
                delta_L += pow(iter_pca->second[num],2);
            }

            for (num = middle; num < 2*middle; ++num)
            {
                delta_W += pow(iter_pca->second[num],2);
            }

            //pGate->deltaL = delta_L;
            //pGate->deltaW = delta_W;

            for (adpt_num = ADAPTIVITY_0; ; ++adpt_num)
            {
                tmp = delta_L * pow(gate_r[adpt_num]*0.15/3, 2) + delta_W * pow(gate_r[adpt_num]*0.08/3, 2);
                    pGate_opt->three_sqrt_delta_LW[adpt_num] = 3*sqrt(tmp);

                if ((ADAPTIVITY_0 + 1) == adpt_num)
                    break;
            }
        }
        }

        overhead = circuit[*ii].get_gate_overhead(*ii, circuit, tile);
        pGate->overhead_area = overhead.first;
        pGate->overhead_power = overhead.second;
    }
}

#endif


/************************************************************************/
/*          Compute the objective value of the Sub problem              */
/************************************************************************/
vector<float> LRACC::_get_inner_lagrangian_obj()
{
    FLT32   leak_sum = 0;
    FLT32   power_sum = 0;
    FLT32   area_sum = 0;
    FLT32   area_lagran = 0;
    FLT32   delay_lagran = 0;
    FLT32   delay_sum = 0;

    TPLGY_VDES  fanoutv;
    UINT8D   is_adapt;
    UINT8D   adpt_num;//whether or not there is adaptivity block
    UINT16D  tile_idx;
    pair<FLT32, FLT32>  overhead;
    FLT32   fanout_c;
    FLT32   fdelay;
    pair<FLT32, FLT32>  fpower;
    FLT32   farea;
    TILE    *pTile;
    GATE    *pGate;

    GATEOPT *pGate_opt;
    GATEOPT *pGate_fout;//the optimal solution for fanout gate
    vector<float>   obj_item;

    VDES_CONTAINER::iterator ii = order.begin();
    for (; ii != order.end(); ++ii)
    {
        if (BIT_TEST(circuit[*ii].info->flag,OUTPUT_LOAD))  continue;//skip the OUTPUT gate

        //to compute the fanout capacitance which is needed in delay and power
        fanout_c = 0;
        pair<TPLGY_AVITER, TPLGY_AVITER> adv_iter = adjacent_vertices(*ii, circuit);
        for (; adv_iter.first != adv_iter.second; ++adv_iter.first)
        {
            fanoutv = *(adv_iter.first);
            if (BIT_TEST(circuit[fanoutv].info->flag, OUTPUT_LOAD))
            {
                adpt_num = ADAPTIVITY_0;//always choose the LOW level
                pGate_fout = gateopt_cast(circuit[fanoutv].optimal_iter[adpt_num]->opt);
            }
            else
            {
#if IS_SIMULTANEOUS
                tile_idx = circuit[fanoutv].info->tile_idx;
                pTile = &(tile[tile_idx]);
                list<COMPONENT_SOLUTION>::iterator  iter_sol = circuit[fanoutv].optimal_iter[pTile->is_adaptable];
                pGate_fout = gateopt_cast(iter_sol->opt);
#else
                pGate_fout = gateopt_cast(circuit[fanoutv].optimal_iter[ADAPTIVITY_0]->opt);
#endif
            }
            fanout_c += pGate_fout->cap;
        }

        //to get the adaptivity block information
        if (BIT_TEST(circuit[*ii].info->flag, INPUT_DRIVER))
        {
            adpt_num = ADAPTIVITY_0;
            is_adapt = 0;
        } 
        else
        {
			if (!BIT_TEST(circuit[*ii].info->flag, FF_VIRTUAL)){
            tile_idx = circuit[*ii].info->tile_idx;
            pTile = &(tile[tile_idx]);

#if IS_SIMULTANEOUS
            adpt_num = pTile->is_adaptable;
#else
            adpt_num = ADAPTIVITY_0;
#endif
            is_adapt = pTile->is_adaptable;
			}
        }

        //computer the area, power and delay for each gate;
        pGate_opt = gateopt_cast(circuit[*ii].optimal_iter[adpt_num]->opt);
        farea = pGate_opt->area;
        fpower = circuit[*ii].calc_gate_power(*pGate_opt, is_adapt, fanout_c, block_prob);
		map<string, vector<double> >::iterator  iter_pca = pca_params.find(circuit[*ii].info->name);
		fdelay = circuit[*ii].calc_gate_delay(*pGate_opt, is_adapt, fanout_c, pca_params, iter_pca) + pGate_opt->offset;
#ifdef  IS_OPTIMAL_CODE
        overhead.first = circuit[*ii].info->overhead_area;
        overhead.second = circuit[*ii].info->overhead_power;
#else
        overhead = circuit[*ii].get_gate_overhead(*ii, circuit, tile);
#endif
        //cout<<"G(" << circuit[*ii].info->name << ") A("<< farea <<") P(" << fpower <<")"<< endl;
        /*if (!BIT_TEST(circuit[*ii].info->flag, INPUT_DRIVER))
        {
            cout<<"G(" << circuit[*ii].info->name<<")";
            cout<<"Size="<<pGate_opt->w<<"\t Vt="<<pGate_opt->vth<<"\t lm="<<circuit[*ii].lm<<endl;
        }*/

        //then sum them to the sum value
        area_sum += farea + is_adapt*overhead.first;
        power_sum += (fpower.first + fpower.second + is_adapt*overhead.second);
        leak_sum += fpower.second;
        delay_sum += fdelay;
        delay_lagran += circuit[*ii].lm * fdelay;
    }
    area_lagran = lm_a*(area_sum-Area_Limit);//MAX_AREA_LIMIT

    //save date to the history struct
    obj_item.resize(SUBPP_MAX_DATA_NUM);
    obj_item[SUBPP_AREA_SUM] = area_sum;
    obj_item[SUBPP_POWER_SUM] = power_sum;
    obj_item[SUBPP_LEAKP_SUM] = leak_sum;
    obj_item[SUBPP_DELAY_SUM] = delay_sum;
    obj_item[SUBPP_AREA_LM] = area_lagran;
    obj_item[SUBPP_DELAY_LM] = delay_lagran;

    return obj_item;
}

/************************************************************************/
/*            Compute the improvement of the Subproblem                 */
/************************************************************************/
bool LRACC::_is_no_improvement(ITER_JUDGE &trade_off_r, vector<FLT32> &new_object)
{
    UINT8D   judge_rtn;
    FLT32   oldst_value = 0;
    FLT32   newst_value = 0;
    vector<FLT32>   old_obj;

    old_obj = trade_off_r.object.back();
    oldst_value += old_obj[SUBPP_AREA_LM];//area
    oldst_value += old_obj[SUBPP_POWER_SUM];//power
    oldst_value += old_obj[SUBPP_DELAY_LM];//delay
    
    newst_value += new_object[SUBPP_AREA_LM];//area
    newst_value += new_object[SUBPP_POWER_SUM];//power
    newst_value += new_object[SUBPP_DELAY_LM];//delay

    if (0 == oldst_value)
    {
        judge_rtn = NO;
    }
    else if (abs(oldst_value-newst_value)/oldst_value > trade_off_r.threshold)
    {
        judge_rtn = NO;
    }
    else
    {
        judge_rtn = YES;
    }

    trade_off_r.object.push_front(new_object);//insert new object
    trade_off_r.object.pop_back();//remove the oldest object

    return  judge_rtn;
}



/* used to relax the consistence */
void    LRACC::_pp_consistence_relax(BOOLEAN is_refine, int _index, int level)
{
	UINT8D       is_adapt;
	UINT8D       adpt_num;
	UINT16D      fan_num;
	pair<FLT32, FLT32>    overhead;
	FLT32       fresist;
	pair<FLT32, FLT32>    fpower;
	FLT32       fdelay;
	FLT32       farea;
	FLT32       fobj_opt;
	FLT32       fobj_pru;
	FLT32       fobj_self;
	FLT32       min_obj;
	FLT32       min_pru;
	list<FANOUT_SUM*>::iterator iter_h;

	UINT8D       lib_index;
	GRADE       adpt_grade;
	GATE        *pGate;
	GATEOPT     *pGate_opt;
	list<FANOUT_SUM*>    fanout_set;//store the combination of the fanout gates
	list<COMPONENT_SOLUTION>    *pfanin_can;//point to the candidates of the fanin NODECOMPONENT
	COMPONENT_SOLUTION  sol_tmp;

	COMPONENT_SOLUTION  multi_opt[2];

	//topological_sort(circuit, back_inserter(c));
	VDES_CONTAINER::iterator ii, iend;

	if (!is_refine)
	{
		 //cout << "CRelax: reverse topological ordering: "<<endl;
		 //cout << level << endl;
	}
	else
	{
		//cout << "CRefine: ====" <<endl;//relaxation or refinement
	}
	
	if (rNodeSetV[level].size() == 0)
	return;
	int si1 = rNodeSetV[level].size() / 4;

	if (_index==0){
	ii = rNodeSetV[level].begin();
	iend = rNodeSetV[level].begin() + si1;
	}
	else if(_index==1){
	ii = rNodeSetV[level].begin() + si1;
	iend = rNodeSetV[level].begin() + 2 * si1;
	}

	else if (_index == 2){
	ii = rNodeSetV[level].begin() + 2 * si1;
	iend = rNodeSetV[level].begin() + 3 * si1;
	}

	else if (_index == 3){
	ii = rNodeSetV[level].begin() + 3 * si1;
	iend = rNodeSetV[level].end();
	}
	/*
	else if (_index == 4){
	ii = rNodeSetV[level].begin() + 4 * si1;
	iend = rNodeSetV[level].begin() + 5 * si1;
	}
	else if (_index == 5){
	ii = rNodeSetV[level].begin() + 5 * si1;
	iend = rNodeSetV[level].begin() + 6 * si1;
	}
	else if (_index == 6){
	ii = rNodeSetV[level].begin() + 6 * si1;
	iend = rNodeSetV[level].begin() + 7 * si1;
	}
	else if (_index == 7){
	ii = rNodeSetV[level].begin() + 7 * si1;
	iend = rNodeSetV[level].end();
	}
	*/
	/*
	//while (k_mux < rNodeSetV[level].size()){
	//boost::mutex i_mux;
	i_mux.lock();
	ii = rNodeSetV[level].begin() + k_mux;
	
	if (k_mux >= rNodeSetV[level].size())
	{	
		ii = rNodeSetV[level].end()-1;
	}
	k_mux++;
	/*
	if (k_mux + 50 < rNodeSetV[level].size()){
		iend = rNodeSetV[level].begin() + k_mux + 50;
		k_mux = k_mux + 50;
	}
	else{
		iend = rNodeSetV[level].end();
		k_mux = rNodeSetV[level].size();
	}
	*
	i_mux.unlock();
	*/
	for (; ii != iend; ++ii)
	//if (ii != rNodeSetV[level].end())
	{
		//cout<< "Test #" << ++tmp << endl;
		if (BIT_TEST(circuit[*ii].info->flag, OUTPUT_LOAD))
			//||BIT_TEST(circuit[*ii].info->flag,INPUT_DRIVER))
		{
			circuit[*ii].touched = true;
			continue;//skip the INPUT and OUTPUT gate
		}
		else
		{
			circuit[*ii].touched = true;
			//During refinement, the optimal solution should be kept
			if (is_refine && 1 < in_degree(*ii, circuit))
			{
				multi_opt[ADAPTIVITY_0] = *circuit[*ii].optimal_iter[ADAPTIVITY_0];
#if IS_SIMULTANEOUS
				multi_opt[ADAPTIVITY_1] = *circuit[*ii].optimal_iter[ADAPTIVITY_1];
#endif
			}
			// Candidates got after each RELAX are invalid and should be cleaned;
			// Besides, since REFINE and RELAX share this function, the optimal solution needed 
			// in REFINE should be kept before cleaning as above and restore at the end
			if (!BIT_TEST(circuit[*ii].info->flag, INPUT_DRIVER))
			{
				circuit[*ii].candidates[ADAPTIVITY_0].clear();
#if IS_SIMULTANEOUS
				circuit[*ii].candidates[ADAPTIVITY_1].clear();
#endif
			}

			//Fanout_merging(*ii, is_refine);
			circuit[*ii].merge_and_prune(*ii, circuit, fanout_set, is_refine);
			//H102_ASSERT(fanout_set[0].size()!=0);
			pGate = gate_cast(circuit[*ii].info);
			//enumerate all the solutions of the fanin gate, not finish
#ifdef  IS_OPTIMAL_CODE
			overhead.first = pGate->overhead_area;
			overhead.second = pGate->overhead_power;
#else
			overhead = circuit[*ii].get_gate_overhead(*ii, circuit, tile);
#endif
			lib_index = circuit[*ii].info->type;
			for (adpt_num = ADAPTIVITY_0;; ++adpt_num)
			{
				map<string, vector<double> >::iterator  iter_pca = pca_params.find(circuit[*ii].info->name);
				if (BIT_TEST(circuit[*ii].info->flag, INPUT_DRIVER))
				{
#if IS_SIMULTANEOUS
					is_adapt = adpt_num;
#else
					is_adapt = 0;
#endif
					min_obj = MAX_OBJECTIVE;
					pfanin_can = &(circuit[*ii].candidates[adpt_num]);
					pGate_opt = gateopt_cast(pfanin_can->begin()->opt);
					//fresist = pGate->r_all(*pGate_opt, adpt_num); //pGate_opt->resist[0];
					iter_h = fanout_set.end();
					list<FANOUT_SUM*>::iterator iter_out = fanout_set.begin();
					for (; iter_out != fanout_set.end(); ++iter_out)
					{
						fdelay = circuit[*ii].calc_gate_delay(*pGate_opt, is_adapt, (*iter_out)->c, pca_params, iter_pca) + pGate_opt->offset;//fresist * iter_out->c
						//pGate_opt->off_power[0] + 0.5*pGate_opt->cap[0];
						fpower = circuit[*ii].calc_gate_power(*pGate_opt, is_adapt, (*iter_out)->c, block_prob);
						farea = pGate_opt->area;
						fobj_self = fpower.first + fpower.second + lm_a*farea + circuit[*ii].lm * fdelay;
						fobj_self += is_adapt*(lm_a*overhead.first + overhead.second);
						fobj_opt = fobj_self + (*iter_out)->obj_opt;
						fobj_pru = fobj_self + (*iter_out)->obj_pru;
						//fobj_pru += adpt_num*(lm_a*overhead.first + overhead.second);

						if (fobj_pru < min_obj)
						{
							min_obj = fobj_opt;
							iter_h = iter_out;
							min_pru = fobj_pru;
						}
					}

					//   H102_ASSERT(iter_h != fanout_set[adpt_num].end());
					pfanin_can->begin()->fanouts = *(*iter_h);
					pfanin_can->begin()->propagates.obj_opt = min_obj;
					pfanin_can->begin()->propagates.obj_pru = min_pru;
				}
				else
				{
#if IS_SIMULTANEOUS
					is_adapt = adpt_num;
#else
					if (!BIT_TEST(circuit[*ii].info->flag, FF_VIRTUAL))
						is_adapt = tile[circuit[*ii].info->tile_idx].is_adaptable;
#endif
					for (fan_num = 0; fan_num < cell_lib[lib_index].size(); ++fan_num)
					{
						min_obj = MAX_OBJECTIVE;
						pGate_opt = gateopt_cast(cell_lib[lib_index][fan_num]);
						//fresist = pGate->r_all(*pGate_opt, adpt_num);
						iter_h = fanout_set.end();
						list<FANOUT_SUM*>::iterator iter_out = fanout_set.begin();

						for (; iter_out != fanout_set.end(); ++iter_out)
						{
							fdelay = circuit[*ii].calc_gate_delay(*pGate_opt, is_adapt, (*iter_out)->c, pca_params, iter_pca) + pGate_opt->offset;//fresist * iter_out->c
							//pGate_opt->off_power[0] + 0.5*pGate_opt->cap[0];
							fpower = circuit[*ii].calc_gate_power(*pGate_opt, is_adapt, (*iter_out)->c, block_prob);
							farea = pGate_opt->area;
							fobj_self = fpower.first + fpower.second + lm_a*farea + circuit[*ii].lm * fdelay;
							fobj_self += is_adapt*(lm_a*overhead.first + overhead.second);
							fobj_opt = fobj_self + (*iter_out)->obj_opt;
							fobj_pru = fobj_self + (*iter_out)->obj_pru;

							if (fobj_pru < min_obj)
							{
								min_obj = fobj_opt;
								iter_h = iter_out;
								min_pru = fobj_pru;
							}
						}

						//      H102_ASSERT(iter_h != fanout_set[adpt_num].end());
#if IS_SIMULTANEOUS
						sol_tmp.adpt_num = adpt_num;
#else
						sol_tmp.adpt_num = is_adapt;
#endif
						sol_tmp.fanouts = *(*iter_h);
						sol_tmp.opt = pGate_opt;
						sol_tmp.propagates.obj_opt = min_obj;
						sol_tmp.propagates.obj_pru = min_pru;
						sol_tmp.propagates.c = pGate_opt->cap;//??????????
						circuit[*ii].candidates[adpt_num].push_back(sol_tmp);
					}
				}

				if (ADAPTIVITY_1 == adpt_num)
					break;
			}

			// Merge the candidate solutions of the current gate level by level
			circuit[*ii].solution_pruning(*ii, circuit);

			// During REFINE, the optimal solution got from last RESTORE should be kept
			if (is_refine && 1 < in_degree(*ii, circuit))
			{
				circuit[*ii].candidates[ADAPTIVITY_0].push_front(multi_opt[ADAPTIVITY_0]);
				circuit[*ii].optimal_iter[ADAPTIVITY_0] = circuit[*ii].candidates[ADAPTIVITY_0].begin();
#if IS_SIMULTANEOUS
				circuit[*ii].candidates[ADAPTIVITY_1].push_front(multi_opt[ADAPTIVITY_1]);
				circuit[*ii].optimal_iter[ADAPTIVITY_1] = circuit[*ii].candidates[ADAPTIVITY_1].begin();
#endif
			}
		}
	}

//}
}


/* used to restore the consistence */
void    LRACC::_pp_consistence_restore(int _index, int level)
{
    UINT16D  tile_num;
    UINT8D   is_adapt;
    UINT8D   adpt_num;
    UINT8D   num;
    list<COMPONENT_SOLUTION>::iterator  iter_sol;
    list<COMPONENT_SOLUTION>::iterator  iter_h;
    FANOUT_SUM  *fanout_sum;
    TILE        *pTile;
    pair<FLT32,FLT32>       overhead;

    FLT32   min_obj;
    FLT32   fanin_pad;//fanin power + area + delay
    FLT32   fobj_pru;
    FLT32   fdelay;
    pair<FLT32,FLT32>   fpower;
    FLT32   farea;
    GATE    *pGate;
    GATEOPT *pGate_opt;//the optimal solution for fanin gate
    GATEOPT *pGate_fout;//the optimal solution for fanout gate

    list<FANIN_SUM>     fanin_sum[2];
    TPLGY_VDES  faninv;

    //topological_sort(circuit, back_inserter(c));
   // cout << "CRestore: topological ordering: " << endl;//Consistency_restoration

#if IS_SIMULTANEOUS
    for (tile_num = 0; tile_num < gGrid_num; ++tile_num)
    {
        pTile = &(tile[tile_num]);
        pTile->obj_value[ADAPTIVITY_0] = 0;
        pTile->obj_value[ADAPTIVITY_1] = 0;
    }
#endif
	VDES_CONTAINER::iterator ii, iend;

	if (NodeSetV[level].size() == 0)
		return;
	int si1 = NodeSetV[level].size() / 2;

	if (_index == 0){
		ii = NodeSetV[level].begin();
		iend = NodeSetV[level].begin() + si1;
	}
	else if (_index == 1){
		ii = NodeSetV[level].begin() + si1;
		iend = NodeSetV[level].end();
	}
	/*
	else if (_index == 2){
	ii = rNodeSetV[level].begin() + 2 * si1;
	iend = rNodeSetV[level].end();
	}

	else if (_index == 3){
	ii = rNodeSetV[level].begin() + 3 * si1;
	iend = rNodeSetV[level].end();
	}

	else if (_index == 4){
	ii = rNodeSetV[level].begin() + 4 * si1;
	iend = rNodeSetV[level].begin() + 5 * si1;
	}
	else if (_index == 5){
	ii = rNodeSetV[level].begin() + 5 * si1;
	iend = rNodeSetV[level].begin() + 6 * si1;
	}
	else if (_index == 6){
	ii = rNodeSetV[level].begin() + 6 * si1;
	iend = rNodeSetV[level].begin() + 7 * si1;
	}
	else if (_index == 7){
	ii = rNodeSetV[level].begin() + 7 * si1;
	iend = rNodeSetV[level].end();
	}
	*/
	//ii = order.rbegin();
	//iend = order.rend();

    //VDES_CONTAINER::reverse_iterator ii = order.rbegin();
    for (; ii != iend; ++ii)
    {
        if (BIT_TEST(circuit[*ii].info->flag, OUTPUT_LOAD)
            || BIT_TEST(circuit[*ii].info->flag, INPUT_DRIVER))
        {
            continue;//skip the INPUT and OUTPUT gate
        }

		if (!BIT_TEST(circuit[*ii].info->flag, FF_VIRTUAL))
		 pTile = &(tile[circuit[*ii].info->tile_idx]);
        pGate = gate_cast(circuit[*ii].info);
        //list and merge the fanin optimum solutioin set
        circuit[*ii].list_and_merge(*ii, circuit, fanin_sum);

        /*/for single fanin gate, inherit its solution from the fanin optimal solution
        if (1 == in_degree(*ii, circuit))
        {
            pair<TPLGY_IEITER, TPLGY_IEITER> vi;
            vi = in_edges(*ii, circuit);
            faninv = source(*vi.first, circuit);

            //to find out the current vertex's index in the fanout sum
            pair<TPLGY_AVITER, TPLGY_AVITER> adv_iter = adjacent_vertices(faninv, circuit);
            for (num = 0; adv_iter.first != adv_iter.second; ++adv_iter.first, ++num)
            {
                if (*ii == *(adv_iter.first))
                    break;
            }
            H102_ASSERT(adv_iter.first != adv_iter.second);

            overhead = circuit[faninv].get_gate_overhead(faninv, circuit, tile);
            for (adpt_num = ADAPTIVITY_0; ;++adpt_num)
            {
                pGate_opt = gateopt_cast(fanin_sum[adpt_num].begin()->fanin_opt->opt);
                fanout_sum = &(fanin_sum[adpt_num].begin()->fanin_opt->fanouts);
                pGate_fout = gateopt_cast(fanout_sum->fanout[num].first);

                fobj_pru = 0;
                fdelay = circuit[faninv].calc_gate_delay(*pGate_opt,
                    adpt_num, pGate_fout->cap, pca_params)+ pGate_opt->offset;
                //(pGate_opt->resist[0] * pGate_fout->cap )
                fpower = circuit[faninv].calc_gate_power(*pGate_opt,
                    adpt_num, pGate_fout->cap, probability);//pGate_opt->off_power[0] + 0.5*pGate_opt->cap;//0;
                farea = pGate_opt->area;
                fobj_pru += (fpower + circuit[faninv].lm*fdelay + lm_a*farea);
                fobj_pru += adpt_num*(lm_a*overhead.first+ overhead.second);

                //to find out the corresponding solution index in the fanout vertex as optimal_index
                iter_sol = circuit[*ii].candidates[adpt_num].begin();
                for (; iter_sol != circuit[*ii].candidates[adpt_num].end(); ++iter_sol)
                {
                    if (pGate_fout == iter_sol->opt)
                        break;
                }
                H102_ASSERT(iter_sol != circuit[*ii].candidates[adpt_num].end());

                fanin_pad = fobj_pru;
                circuit[*ii].optimal_iter[adpt_num] = iter_sol;
                pTile->obj_value[adpt_num] += fanin_pad;

                if (ADAPTIVITY_1 == adpt_num)
                    break;
            }
        }
        else*/
        {
            for (adpt_num = ADAPTIVITY_0;; ++adpt_num)
            {
#if IS_SIMULTANEOUS
                is_adapt = adpt_num;
#else
				if (!BIT_TEST(circuit[*ii].info->flag, FF_VIRTUAL))
					is_adapt = pTile->is_adaptable;
				else
					is_adapt = 0;
#endif
                iter_sol = circuit[*ii].candidates[adpt_num].begin();
                fanin_pad = 0;
                min_obj = MAX_OBJECTIVE;
                iter_h = circuit[*ii].candidates[adpt_num].end();
                for (; iter_sol != circuit[*ii].candidates[adpt_num].end(); ++iter_sol)
                {
                    fobj_pru = 0;
                    pGate_fout = gateopt_cast(iter_sol->opt);
                    list<FANIN_SUM>::iterator iter_in = fanin_sum[adpt_num].begin();
                    for (; iter_in != fanin_sum[adpt_num].end(); ++iter_in)
                    {
                        faninv = iter_in->faninv;
#ifdef  IS_OPTIMAL_CODE
                        overhead.first = pGate->overhead_area;
                        overhead.second = pGate->overhead_power;
#else
                        overhead = circuit[faninv].get_gate_overhead(faninv, circuit, tile);
#endif
                        pGate_opt = gateopt_cast(iter_in->fanin_opt->opt);
						map<string, vector<double> >::iterator  iter_pca = pca_params.find(circuit[*ii].info->name);
                        fdelay = circuit[faninv].calc_gate_delay(*pGate_opt,
							is_adapt, pGate_fout->cap, pca_params, iter_pca) + pGate_opt->offset;
                        //pGate_opt->resist[0] * pGate_fout_opt->cap
                        fpower = circuit[faninv].calc_gate_power(*pGate_opt,
                            is_adapt, pGate_fout->cap, block_prob);//pGate_opt->off_power[0] + 0.5*pGate_opt->cap[0];
                        farea = pGate_opt->area;
                        fobj_pru += (fpower.first + fpower.second + circuit[faninv].lm*fdelay + lm_a*farea);
                        fobj_pru += is_adapt*(lm_a*overhead.first + overhead.second);
                    }

                    // to choose the optimize solution, order by the obj + sum of fanin(power + area + delay)
                    if (fobj_pru + iter_sol->propagates.obj_pru < min_obj)
                    {
                        min_obj = fobj_pru + iter_sol->propagates.obj_pru;
                        fanin_pad = fobj_pru;
                        iter_h = iter_sol;
                    }
                }

               // H102_ASSERT(iter_h != circuit[*ii].candidates[adpt_num].end());
                circuit[*ii].optimal_iter[adpt_num] = iter_h;
#if IS_SIMULTANEOUS
                pTile->obj_value[adpt_num] += fanin_pad;
#endif
                if (ADAPTIVITY_1 == adpt_num)
                    break;
            }
        }

#if IS_SIMULTANEOUS
        //after finish the last gate in the adaptivity block, judge whether adaptivity unit is better 
        if (*ii == pTile->components[0])
        {
            if (pTile->obj_value[ADAPTIVITY_0] > pTile->obj_value[ADAPTIVITY_1])
            {
                pTile->is_adaptable = YES;
            }
            else
            {
                pTile->is_adaptable = NO;
            }
        }
#endif
    }
}

/******************************************************************************/
/* Initialize the Timing Analysis which parses the required files etc.        */
/******************************************************************************/
double output_req = 100 ;
void LRACC::init_timing_analysis(string benchmark, string placement, double output_required_time)
{
	
	project = new ocl_opt;
	//string benchfile = "E:/Project/ocl_opt/ocl_opt/benchmarks/fft.v";
	bool valid = FALSE;
	valid = project->init_lracc_database();
	valid = project->load_graph_and_tiles(benchmark);

	int no_of_gates;
	no_of_gates = igraph_vcount(&(project->circuitGraph));
	int no_of_edges;
	no_of_edges = igraph_vector_size(&project->edges);

	//project.load_celllib(const_cast<char *>("E:/Project/ocl_opt/ocl_opt/benchmarks/contest.lib"));
	//project.load_at_raq();
	//project.update_gates_info();
	
	//string place_bench = "E:/Project/ocl_opt/ocl_opt/benchmarks/fft.pl";

	map<string, pair<double, double> > mu_sigma_local;
	map<string, pair<double, double> > arrival_time_slack_local;
	project->init_timing(mu_sigma_local, arrival_time_slack_local, placement, output_required_time, gate_sizes, adpt_block_exist, no_of_adaptive_blocks, false); //adpt_blk_exist
	this->mu_sigma = mu_sigma_local;
	this->arrival_time_slack = arrival_time_slack_local;



	map<string, vector<double> > pca_param_local;
	map<string, vector<float> >  plac_info_local;
	project->pass_pca_plac_info(pca_param_local, plac_info_local);
	this->pca_params = pca_param_local;
	this->plac_info = plac_info_local;

	map<int, vector<float> >   block_prob_local;
	project->calculate_prob(block_prob_local);
	this->block_prob = block_prob_local;
		//int actual_blocks;
	//	int x,y;

//		x = 4; y = 4;//x = 6; y =6;

   // Timing_Analysis::Init_Timing_Analysis(benchmark, placement, gate_sizes, output_required_time,x,y,&actual_blocks, true /*true = only FBB*/);


		//no_of_adaptive_blocks = x*y;
		output_req = output_required_time;
    //LRACC_ssta();
    //Timing_Analysis::get_pca_params_plac_info(pca_param_local, plac_info_local);
    //get_probability();
}

void LRACC::dsta()
{
		LRACC_ssta();
    critical_path_delay = Timing_Analysis::DSTA();
		cout << "THE CRITICAL PATH DELAY = " << critical_path_delay << endl;\
		//int wait;
		//cin >> wait;
}

void LRACC::get_power(){
	//power_mu_sigma = Timing_Analysis::get_power(gate_sizes, adpt_block_exist, for_get_power);
	//cout << "THE power distribution = " << power_mu_sigma.first << "   " << power_mu_sigma.second << endl;
	//int wait;
	//cin >> wait;
}

void LRACC::get_timing_yield(){
		double timing_yield;
		timing_yield = Timing_Analysis::get_timing_yield_temp(adpt_block_exist, output_req, gate_sizes);
		cout << "TIMING YIELD  = " << timing_yield;
		//int wait;
		//cin >> wait;
}

void LRACC::LRACC_ssta()
{
    //cout<< "LRACC Primal finished"<< endl;
	map<string, pair<double, double> > mu_sigma_local;
	map<string, pair<double, double> > arrival_time_slack_local;
	project->SSTA(mu_sigma_local,  arrival_time_slack_local, gate_sizes, adpt_block_exist);

//    Timing_Analysis::SSTA(mu_sigma_local, arrival_time_slack_local,gate_sizes, output_req, adpt_block_exist);

    this->mu_sigma = mu_sigma_local;
    this->arrival_time_slack = arrival_time_slack_local;


	get_probability();
	/*	
		cout << "Printint hte values "  << endl;
     for(std::map<string,pair<double,double> >::iterator it = arrival_time_slack_local.begin();
                                                      it != arrival_time_slack_local.end();  ++it){
      cout << it->first  << "  "<< (it->second).first << "   " << (it->second).first << endl;
     }
 
     int wait;
     cin >> wait;
*/

}

void LRACC::get_probability(){
    //map<int, double> probability_local;
		//map<int, string> prob_type_local;
    map<int, vector<float> >   block_prob_local;
	project->calculate_prob(block_prob_local);
    //LRACC_ssta();
    //Timing_Analysis::get_probability(block_prob_local);//probability_local, prob_type_local);
    //this->probability = probability_local;
		//this->prob_type = prob_type_local;
    this->block_prob = block_prob_local;
}

void LRACC::calc_prob()
{
		map<int,bool> adpt_tmp = adpt_block_exist;
    for(TPLGY_VITER viter=boost::vertices(circuit).first;viter!=boost::vertices(circuit).second;viter++)
    {
        if(BIT_TEST(circuit[*viter].info->flag,INPUT_DRIVER)||BIT_TEST(circuit[*viter].info->flag,OUTPUT_LOAD)) continue;
        std::string nm=circuit[*viter].info->name;
        
        GATEOPT   *pGate = NULL;
        if (1 != iter_lagrangian)
        {
            pGate = gateopt_cast(circuit[*viter].optimal_iter[0]->opt);
        }
        else
        {
            H102_ASSERT(1 == iter_lagrangian);
            pGate = gateopt_cast(cell_lib[circuit[*viter].info->type][0]);
        }

        std::vector<float> tmp;
        tmp.push_back(pGate->w);
        tmp.push_back(pGate->resist[0]);
        tmp.push_back(pGate->cap);
        tmp.push_back(pGate->vth);
        gate_sizes[nm] = tmp;
    }
    for (UINT8D num = 0; num < MAX_TUNING_UNIT; ++num)
    {
        adpt_block_exist[num] = 0;//all use 0 adaptivity
    }
    LRACC_ssta();	
    //*
    block_prob.clear();
    for(int i=0;i<tile.size();i++)
    {
        float min_slack = 1e100;
        for(int j=0;j<tile[i].components.size();j++)
        {
            string nm = circuit[tile[i].components[j]].info->name;
            if(arrival_time_slack[nm].second<min_slack) min_slack = arrival_time_slack[nm].second;
        }
        std::vector<float> tmp(3,0);
        if(min_slack>=40)
        {
            tmp[ZERO_GRADE] = 1;
            tmp[FBB_GRADE_1] = 0;
            tmp[RBB_GRADE_1] = 0;
        }
        else if(min_slack>0&&min_slack<40)
        {
            tmp[ZERO_GRADE] = (40-min_slack)/40;
            tmp[RBB_GRADE_1] = 1-tmp[ZERO_GRADE];
            tmp[FBB_GRADE_1] = 0;
        }
        else if(min_slack==0)
        {
            tmp[ZERO_GRADE] = 1;
            tmp[RBB_GRADE_1] = 0;
            tmp[FBB_GRADE_1] = 0;
        }
        else if(min_slack>-40&&min_slack<0)
        {
            tmp[ZERO_GRADE] = (40+min_slack)/40;
            tmp[RBB_GRADE_1] = 0;
            tmp[FBB_GRADE_1] = 1-tmp[ZERO_GRADE];	
        }
        else
        {
            tmp[ZERO_GRADE] = 0;
            tmp[RBB_GRADE_1] = 0;
            tmp[FBB_GRADE_1] = 1;
        }
        block_prob[i] = tmp;
    }//*/
		adpt_block_exist = adpt_tmp;
}

int LRACC::incr_lag_num()
{
    return  ++iter_lagrangian;
}

VOID    LRACC::LRACC_primal()
{
    BOOLEAN refine_iter = FALSE;
    vector<FLT32>   object;
    UINT16D    num;
    FLT32     Min_obj;
    trade_off_r.initial(0.05);//initialize the record history
	clock_t t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;

	//gate level schedule...
	rscheduling();
	scheduling();
	boost::thread_group a;
    while (TRUE) 
    {
		cout << "Begin relax!" << endl;
		t1 = clock();
		
		for (int i = 0; i < rNodeSetV.size(); ++i){
			for (int j = 0; j < rNodeSetV[i].size(); ++j){
				circuit[rNodeSetV[i][j]].touched = false;
			}
			k_mux = 0;
			//a.create_thread(boost::bind(&LRACC::_pp_consistence_relax, this, refine_iter, 0, i));
			//a.create_thread(boost::bind(&LRACC::_pp_consistence_relax, this, refine_iter, 1, i));
			//a.create_thread(boost::bind(&LRACC::_pp_consistence_relax, this, refine_iter, 2, i));
			//a.create_thread(boost::bind(&LRACC::_pp_consistence_relax, this, refine_iter, 3, i));
		//	a.create_thread(boost::bind(&LRACC::_pp_consistence_relax, this, refine_iter, 4, i));
		//	a.create_thread(boost::bind(&LRACC::_pp_consistence_relax, this, refine_iter, 5, i));
		//	a.create_thread(boost::bind(&LRACC::_pp_consistence_relax, this, refine_iter, 6, i));
		//	a.create_thread(boost::bind(&LRACC::_pp_consistence_relax, this, refine_iter, 7, i));

		boost::thread a1 (&LRACC::_pp_consistence_relax, this, refine_iter, 0, i);
		
		boost::thread a2(&LRACC::_pp_consistence_relax, this, refine_iter, 1, i);
		//a1.join();
		//a2.join();
		boost::thread a3(&LRACC::_pp_consistence_relax, this, refine_iter, 2, i);
		//a3.join();
		boost::thread a4(&LRACC::_pp_consistence_relax, this, refine_iter, 3, i);
		//a4.join();
		//boost::thread a5(&LRACC::_pp_consistence_relax, this, refine_iter, 4, i);
		//a5.join();
		//boost::thread a6(&LRACC::_pp_consistence_relax, this, refine_iter, 5, i);
		//a6.join();
		//boost::thread a7(&LRACC::_pp_consistence_relax, this, refine_iter, 6, i);
		//a7.join();
		//boost::thread a8(&LRACC::_pp_consistence_relax, this, refine_iter, 7, i);
		
		a1.join();
		a2.join();
		a3.join();
		a4.join();
		//a5.join();
		//a6.join();
		//a7.join();
		//a8.join();
		//boost::thread a2(&LRACC::_pp_consistence_relax, this, refine_iter, 1);
//		a2.join();
			//a.join_all();
			cout << "Level over, k is: " << k_mux << endl;
			for (int j = 0; j < rNodeSetV[i].size(); ++j){
				if (!circuit[rNodeSetV[i][j]].touched)
				{
					cout << j << endl;
					getchar();
				}
			}
			getchar();
		}
		

		//_pp_consistence_relax(refine_iter, 0, 0);
		t2 = clock();
	//	std::cout << "Joining threads....." << endl;
		float d = (float)t2 - (float)t1;
		cout << "Test for one iter!! = " << d / CLOCKS_PER_SEC << endl;
		
		getchar();
		//boost::thread b1 (&LRACC::_pp_consistence_restore, this, 0);
		//b1.join();
		//boost::thread b2(&LRACC::_pp_consistence_restore, this, 1);
		//b2.join();

		for (int i = 0; i < NodeSetV.size(); ++i){
			boost::thread a1(&LRACC::_pp_consistence_restore, this, 0, i);
			a1.join();
			boost::thread a2(&LRACC::_pp_consistence_restore, this, 1, i);
			
			a2.join();
		}

        object = _get_inner_lagrangian_obj();
        /* the sequence of the vector<FLT32> is area, power and delay */
        cout<< "Area(" <<object[SUBPP_AREA_SUM] <<") \t";
        cout<< "Power(" <<object[SUBPP_POWER_SUM]<<") \t";
        cout<< "Leak(" <<object[SUBPP_LEAKP_SUM]<<") \t";
        cout<< "Delay(" <<object[SUBPP_DELAY_SUM]<<") \t";
        cout<< "Area_lm(" <<object[SUBPP_AREA_LM]<<") \t";
        cout<< "Delay_lm(" <<object[SUBPP_DELAY_LM]<<") \t";

        Min_obj = object[SUBPP_POWER_SUM]+object[SUBPP_DELAY_LM]+object[SUBPP_AREA_LM];
        cout<< "Min=" << Min_obj <<endl;

        for (num = 0; num < tile.size(); ++num)
        {
            //cout<<"Tile("<<num<<"): adapt= "<<tile[num].is_adaptable <<endl;
        }

        if (_is_no_improvement(trade_off_r, object))
        {
            break;//use threshold to break
        }
        else
        {
            refine_iter = TRUE;
        }
    }
    //later, update solution set here

#if 0
    while (TRUE)
    {
        _pp_consistence_relax(NO);
        _pp_consistence_restore();

        while (TRUE)
        {
            _pp_consistence_relax(YES);
            _pp_consistence_restore();

            if (_is_no_improvement(&trade_off_r, object))
                break;//use threshold to break later
        }
        break;//update later
    }
#endif

    cout<< "LRACC Primal finished; lag= \n"<<iter_lagrangian<< endl;
}

void LRACC::check_gate_solution()
{
    FLT32   a = 0;
    UINT8D   adpt_num = 0;
    UINT16D  tile_num = 0;
    UINT32D  gate_num = 0;
    GATEOPT *pGate_opt;
    TILE    *pTile;
    TPLGY_VDES  fgate;

    for (; tile_num<gGrid_num; ++tile_num)
    {
        pTile = &(tile[tile_num]);
        adpt_num = pTile->is_adaptable;
        for (gate_num = 0; gate_num < pTile->components.size(); ++gate_num)
        {
            fgate = pTile->components[gate_num];
            pGate_opt = gateopt_cast(circuit[fgate].optimal_iter[adpt_num]->opt);
            a += pGate_opt->area;

            if (arrival_time_slack[circuit[fgate].info->name].second<0)
            {cout<<"G("<<circuit[fgate].info->name <<") timing infeasible!"<<endl;}
        }

        a+= adpt_num*pTile->tuner.area();
    }

    if (a>Area_Limit)
    {cout<<"Solutions Area infeasible!"<<endl;}
}

/*
//this function is to record the best gate solution till now
bool    LRACC::record_gate_solution()
{
    UINT8D   is_adapt;
    UINT8D   adpt_num;//whether or not there is adaptivity block
    FLT32   power_sum = 0;
    FLT32   area_sum = 0;
    FLT32   leak_sum = 0;
    UINT32D  tmp_flag = 0;
    FLT32   adpt_area = 0;
    FLT32   adpt_power = 0;

    TPLGY_VDES  fanoutv;
    FLT32   fanout_c;
    FLT32   fslack;
    pair<FLT32,FLT32>   fpower;
    FLT32   farea;

    UINT16D  tile_num;
    TILE    *pTile;

    GATEOPT *pGate_opt;
    GATEOPT *pGate_fout;//the optimal solution for fanout gate
    FLT32   min_slack = 1e6;

    map<int, bool>        tile_tmp;
    GATE_FINAL          gate_tmp;
    list<GATE_FINAL>    fgate_tmp;
    //tile_tmp.resize(gGrid_num);

    VDES_CONTAINER::iterator ii = order.begin();
    for (; ii != order.end(); ++ii)
    {
        if (BIT_TEST(circuit[*ii].info->flag,OUTPUT_LOAD))  continue;//skip the OUTPUT gate

        //to compute the fanout capacitance which is needed in delay and power
        fanout_c = 0;
        pair<TPLGY_AVITER, TPLGY_AVITER> adv_iter = adjacent_vertices(*ii, circuit);
        for (; adv_iter.first != adv_iter.second; ++adv_iter.first)
        {
            fanoutv = *(adv_iter.first);
            if (BIT_TEST(circuit[fanoutv].info->flag, OUTPUT_LOAD))
            {
                adpt_num = ADAPTIVITY_0;//always choose the LOW level
                pGate_fout = gateopt_cast(circuit[fanoutv].optimal_iter[adpt_num]->opt);
            }
            else if (1 != iter_lagrangian)
            {
#if IS_SIMULTANEOUS
                pTile = &(tile[circuit[fanoutv].info->tile_idx]);
                list<COMPONENT_SOLUTION>::iterator  iter_sol = circuit[fanoutv].optimal_iter[pTile->is_adaptable];
                pGate_fout = gateopt_cast(iter_sol->opt);
#else
                pGate_fout = gateopt_cast(circuit[fanoutv].optimal_iter[ADAPTIVITY_0]->opt);
#endif
            }
            else
            {
                pGate_fout = gateopt_cast(cell_lib[circuit[fanoutv].info->type][0]);
            }
            fanout_c += pGate_fout->cap;
        }

        //to get the solution of the current gate
        pGate_opt = NULL;
        if (BIT_TEST(circuit[*ii].info->flag, INPUT_DRIVER))
        {
            is_adapt = 0;
            adpt_num = ADAPTIVITY_0;
            pGate_opt = gateopt_cast(circuit[*ii].optimal_iter[adpt_num]->opt);
        } 
        else if (1 != iter_lagrangian)
        {
            pTile = &(tile[circuit[*ii].info->tile_idx]);
            is_adapt = pTile->is_adaptable;
#if IS_SIMULTANEOUS
            adpt_num = pTile->is_adaptable;
#else
            adpt_num = ADAPTIVITY_0;
#endif
            pGate_opt = gateopt_cast(circuit[*ii].optimal_iter[adpt_num]->opt);
        }
        else
        {
            //for the first iteration, no adaptivity blk
            is_adapt = 0;//pTile->is_adaptable;
            adpt_num = ADAPTIVITY_0;
            pGate_opt = gateopt_cast(cell_lib[circuit[*ii].info->type][0]);
        }

        //computer the area, power and delay for current gate;
        farea = pGate_opt->area;
        fpower = circuit[*ii].calc_gate_power(*pGate_opt, is_adapt, fanout_c, block_prob);
        fslack = arrival_time_slack[circuit[*ii].info->name].second;
        if (min_slack > fslack) min_slack = fslack;

        //cout<<"G(" << circuit[*ii].info->name << ") A("<< farea <<") P(" << fpower <<")"<< endl;
        //if (!BIT_TEST(circuit[*ii].info->flag, INPUT_DRIVER))
        {
            //cout<<"G(" << circuit[*ii].info->name<<")\t";
            //cout<<"Size="<<pGate_opt->w<<"\t Vt="<<pGate_opt->vth<<"\t Area="<<farea;//
            //cout<<"\t Power="<<fpower.first+fpower.second<<"\t Slack="<<fslack<<endl;
            gate_tmp.fgate = *ii;
            gate_tmp.pGate = pGate_opt;
            gate_tmp.is_adpt = is_adapt;
            fgate_tmp.push_back(gate_tmp);
            if (fslack<0)
            {
                BIT_SET(tmp_flag, TIME_INFEASIBLE);
            }
        }

        //then sum them to the sum value
        area_sum += farea;// + adpt_num*overhead.first;
        power_sum += (fpower.first + fpower.second);// + adpt_num*overhead.second;
        leak_sum += fpower.second;//sum the leakage power
    }

    // Update the area and power info caused by adding the adaptivity block
    for (tile_num = 0 ; tile_num < gGrid_num; ++tile_num)
    {
        pTile = &(tile[tile_num]);
        tile_tmp[tile_num] = pTile->is_adaptable;
        if (tile_tmp[tile_num])
        {
            adpt_area += pTile->tuner.area();
            area_sum += pTile->tuner.area();
            adpt_power += pTile->tuner.power();
            power_sum += pTile->tuner.power();
        }
    }
    if (area_sum > Area_Limit)
    {
        cout<<"Solutions Area infeasible!"<<endl;
        BIT_SET(tmp_flag, AREA_INFEASIBLE);
    }
    cout<<"Total Area = "<<area_sum<<"; Total Power = "<<power_sum<<"; Total leak = "<<leak_sum<<"; Min slack = "<<min_slack<<endl;
		total_a = area_sum;

////////////////////cutting plane//////////////////////////////////
		total_p = power_sum;
////////////////////////////////////////////////////////////////////

    // compare the data with the history record, update it if current solution is better
    BOOLEAN     need_update = NO;
    if (INIT_INFEASIBLE == final_result.solution_flag)
    {
        need_update = YES;
    }
    else
    {
        if (BIT_TEST(final_result.solution_flag, TIME_INFEASIBLE))
        {
            if (min_slack > final_result.min_slack)
            {
                need_update = YES;
            } 
            else
            {
                need_update = NO;
            }
        } 
        else
        {
            if (BIT_TEST(tmp_flag, TIME_INFEASIBLE))
            {
                need_update = NO;
            } 
            else
            {
                if (BIT_TEST(final_result.solution_flag, AREA_INFEASIBLE))
                {
                    if (area_sum < final_result.area_sum)
                    {
                        need_update = YES;
                    }
                    else
                    {
                        need_update = NO;
                    }
                } 
                else
                {
                    if (BIT_TEST(tmp_flag, AREA_INFEASIBLE))
                    {
                        need_update = NO;
                    } 
                    else
                    {
                        if (power_sum < final_result.power_sum)
                        {
                            need_update = YES;
                        } 
                        else
                        {
                            need_update = NO;
                        }
                    }
                }
            }
        }
    }

    dest_file<<iter_lagrangian<<"\t";
    dest_file<<power_sum<<"\t";
    dest_file<<leak_sum<<"\t";
    dest_file<<area_sum<<"\t";
    dest_file<<min_slack<<"\t";
    int dest_adpt_num = 0;
    for (tile_num = 0; tile_num < gGrid_num; ++tile_num)
        dest_adpt_num += tile_tmp[tile_num];
    dest_file<<dest_adpt_num<<"\t";
    dest_file<<adpt_power<<"\t";
    dest_file<<adpt_area<<endl;

    if (need_update)
    {
        cout<<"The #"<<iter_lagrangian<<endl;
        for (tile_num = 0; tile_num < gGrid_num; ++tile_num)
        {
            final_result.tile_prob[tile_num] = block_prob[tile_num];
            final_result.tile_adpt[tile_num] = tile_tmp[tile_num];
            cout<<"Record Adapt("<<tile_num<<") is_Adpt"<<tile_tmp[tile_num]<<";";
            for (UINT8D apt_index = ZERO_GRADE; apt_index < MAX_BB_GRADE; ++apt_index)
                cout<<"P["<<(UINT16D)apt_index<<"]="<<block_prob[tile_num][apt_index]<<endl;
        }
        cout<<endl;

        final_result.iter_lm = iter_lagrangian;
        final_result.solution_flag = tmp_flag;
        final_result.area_sum = area_sum;
        final_result.power_sum = power_sum;
        final_result.leak_sum = leak_sum;
        final_result.tile_adpt = tile_tmp;
        final_result.fopt_set = fgate_tmp;
        final_result.min_slack = min_slack;
        //final_result.max_slack = max_slack;
        final_result.adpt_area = adpt_area;
        final_result.adpt_power = adpt_power;
    }

    return need_update;
}
*

void    LRACC::print_lracc_result()
{
    vector<float>   param;
    param.resize(SSTA_PARAM_MAX);

    //if (!final_result.solution_flag)
    {
        cout<<endl;
        cout<<"Find solution at iteration#= "<< final_result.iter_lm<<endl;
        cout<<"solution: power="<< final_result.power_sum<<"\t leak="<<final_result.leak_sum;
        cout<<"\t area="<<final_result.area_sum<<"\t min_slack="<<final_result.min_slack;
        cout<<"\t adpt_area="<<final_result.adpt_area<<"\t adpt_power="<<final_result.adpt_power<<endl;
        for (UINT16D tile_num = 0; tile_num < gGrid_num; ++tile_num)
        {
            cout<<"Tile("<<tile_num<<") #gate="<<tile[tile_num].components.size()<<"\t adaptivity = "<<final_result.tile_adpt[tile_num]<<endl;

            //set the adpt_block and probability info for the final timing yield
            adpt_block_exist[tile_num] = final_result.tile_adpt[tile_num];

            //for (UINT8DD apt_index = FBB_GRADE_1; apt_index < MAX_BB_GRADE; ++apt_index)
            for_get_power[tile_num] = final_result.tile_prob[tile_num];
        }

        list<GATE_FINAL>::iterator final_iter = final_result.fopt_set.begin();
        for (;final_iter != final_result.fopt_set.end(); ++final_iter)
        {
            //cout<<"G("<<circuit[final_iter->fgate].info->name<< ")\t Size="<<final_iter->pGate->w<<"\t Vt="<<final_iter->pGate->vth<<endl;
        
            //update the gate size info for the final timing yield
            if (BIT_TEST(circuit[final_iter->fgate].info->flag, INPUT_DRIVER)
                ||BIT_TEST(circuit[final_iter->fgate].info->flag, OUTPUT_LOAD)) continue;

            param[SSTA_PARAM_IO1] = final_iter->pGate->w;
            //always give the resistance of the Low Grade
            param[SSTA_PARAM_IO2] = final_iter->pGate->resist[ZERO_GRADE];
            param[SSTA_PARAM_IO3] = final_iter->pGate->cap;
            param[SSTA_PARAM_IO4] = final_iter->pGate->vth;
            param[SSTA_PARAM_IO5] = final_iter->pGate->off_power[ZERO_GRADE];
            //H102_ASSERT(param[SSTA_PARAM_IO5]);
            gate_sizes[circuit[final_iter->fgate].info->name] = param;
        }

        //get_power();
        //get_timing_yield();
    }

    cout <<endl;
    if (0 == final_result.solution_flag)
    {cout<<"Solution feasible"<<endl;}
    else
    {
        if (BIT_TEST(final_result.solution_flag, AREA_INFEASIBLE))
				{cout<<"Area Infeasible"<<endl;}
        if (BIT_TEST(final_result.solution_flag, TIME_INFEASIBLE))
        {cout<<"Timing Infeasible"<<endl;}
    }

    write_solution_to_file();

    dest_file<<final_result.iter_lm<<"\t";
    dest_file<<final_result.iter_lm<<"\t";
    dest_file<<final_result.iter_lm<<"\t";
    dest_file<<final_result.iter_lm<<"\t";
    dest_file<<final_result.iter_lm<<"\t";
    dest_file<<final_result.iter_lm<<"\t";
    dest_file<<final_result.iter_lm<<"\t";
    dest_file<<final_result.iter_lm<<endl;
}

*
VOID    LRACC::write_solution_to_file()
{
    string  filename("solution.tmp");
    UINT32D  num = 0;
    UINT32D  max_num = 0;
    COMPONENTTYPE type;
    GATEOPT *pGate_opt = NULL;

    ofstream tmp;
    tmp.open(filename.c_str(),ios::trunc);
    tmp.close();

    list<GATE_FINAL>::iterator final_iter = final_result.fopt_set.begin();
    for (;final_iter != final_result.fopt_set.end(); ++final_iter)
    {
        if (BIT_TEST(circuit[final_iter->fgate].info->flag, OUTPUT_LOAD)
            ||BIT_TEST(circuit[final_iter->fgate].info->flag, INPUT_DRIVER)) continue;

        pGate_opt = final_iter->pGate;
        type = circuit[final_iter->fgate].info->type;

        Store(filename, final_iter->fgate);
        max_num = cell_lib[type].size();
        for (num = 0; num < max_num; ++num)
        {
            if (pGate_opt == cell_lib[type][num])
                break;
        }
        H102_ASSERT(max_num != num);

        Store(filename, num);
        //cout<<"name = "<<circuit[final_iter->fgate].info->name<<"\t opt = "<< num<<endl;
        /*Store(filename, pGate_opt->w);
        Store(filename, pGate_opt->vth);
        Store(filename, pGate_opt->cap);
        Store(filename, pGate_opt->area);
        Store(filename, pGate_opt->off_power);
        Store(filename, pGate_opt->resist);*
    }
}


FLT32   LRACC::read_solution_by_file()
{
    vector<FLT32>   param;
    FLT32   area_sum = 0;
    GATEOPT Gate_tmp;
    COMPONENT_SOLUTION  sol_tmp;
    GATEOPT *pGate_opt = &Gate_tmp;
    UINT8D   adpt_num;
    UINT64D  fopt_num;
    UINT32D  num = 0;
    UINT32D  max_num = 0;
    string filename("solution.tmp");
    TPLGY_VDES  fves;

    param.resize(SSTA_PARAM_MAX);
    VDES_CONTAINER::iterator ii = order.begin();
    for (; ii != order.end(); ++ii)
    {
        if (BIT_TEST(circuit[*ii].info->flag,OUTPUT_LOAD)
            ||BIT_TEST(circuit[*ii].info->flag,INPUT_DRIVER))  continue;//skip the IN/OUTPUT gate
        
        Read(filename, ++num, fves);
        H102_ASSERT(fves == *ii);// used benchmark file not match with input file
        Read(filename, ++num, fopt_num);
        max_num = cell_lib[circuit[*ii].info->type].size();
        H102_ASSERT(max_num > fopt_num);
        pGate_opt = gateopt_cast(cell_lib[circuit[*ii].info->type][fopt_num]);
        area_sum += pGate_opt->area;

        for (adpt_num = ADAPTIVITY_0; ; ++adpt_num)
        {
            sol_tmp.adpt_num = adpt_num;
            //sol_tmp.fanouts = *iter_h;
            sol_tmp.opt = pGate_opt;
            //sol_tmp.propagates.obj_opt = min_obj;
            //sol_tmp.propagates.obj_pru = min_pru;
            sol_tmp.propagates.c = pGate_opt->cap;
            circuit[*ii].candidates[adpt_num].push_back(sol_tmp);

            if (ADAPTIVITY_1 == adpt_num)
                break;
        }
        cout<<"name = "<<circuit[*ii].info->name<<"\t opt = "<<fopt_num<<endl;
        param[SSTA_PARAM_IO1] = pGate_opt->w;
        //always give the resistance of the Low Grade
        param[SSTA_PARAM_IO2] = pGate_opt->resist[ZERO_GRADE];
        param[SSTA_PARAM_IO3] = pGate_opt->cap;
        param[SSTA_PARAM_IO4] = pGate_opt->vth;
        param[SSTA_PARAM_IO5] = pGate_opt->off_power[ZERO_GRADE];
        //H102_ASSERT(param[SSTA_PARAM_IO5]);
        gate_sizes[circuit[*ii].info->name] = param;
    }

    return  area_sum;
}

*/

VOID LRACC::rscheduling(){

	VDES_CONTAINER::iterator nodeiter; //reversed traversal
	pair<TPLGY_IEITER, TPLGY_IEITER> edgeiter;
	VDES_CONTAINER tmpV;
	TPLGY_VDES  faninv;
	//---------find the primary input-----------
	for (nodeiter = order.begin(); nodeiter != order.end(); ++nodeiter){
		circuit[*nodeiter].visit_times = 0;
		if (BIT_TEST(circuit[*nodeiter].info->flag, OUTPUT_LOAD)){
			circuit[*nodeiter].ready = true;
			tmpV.push_back(*nodeiter);
		}
		else{
			circuit[*nodeiter].ready = false;
		}
	}
	
	
	//----schedule-------
	while (!tmpV.empty()){

		if (tmpV.size() <= MAX_GROUP){
			
			for (nodeiter = tmpV.begin(); nodeiter != tmpV.end(); ++nodeiter)
			{
				for (edgeiter = in_edges(*nodeiter, circuit); edgeiter.first != edgeiter.second; ++edgeiter.first){
					
					faninv = source(*edgeiter.first, circuit);
					circuit[faninv].visit_times++;
				}
			}
			rNodeSetV.push_back(tmpV);
			tmpV.clear();
		}
		else{
			
			/*
			int tmp_max = -1;
			int max_idx = 0;
			TPLGY_VDES tmpN;
			for (int i = 0; i < tmpV.size(); ++i){
				for (int j = i; j < tmpV.size(); ++j){
					//cout << circuit[tmpV[j]].candidates[0].size() << endl;
					if (circuit[tmpV[j]].candidates[0].size() > tmp_max){
						max_idx = j;
						tmp_max = circuit[tmpV[j]].candidates[0].size();
					}
				}
				tmpN = tmpV[max_idx];
				tmpV[max_idx] = tmpV[i];
				tmpV[i] = tmpN;
				tmp_max = -1;
			}
			*/
			for (nodeiter = tmpV.begin(); nodeiter < tmpV.begin() + MAX_GROUP; ++nodeiter)
			{
				for (edgeiter = in_edges(*nodeiter, circuit); edgeiter.first != edgeiter.second; ++edgeiter.first){
					faninv = source(*edgeiter.first, circuit);
					circuit[faninv].visit_times++;
				}
			}

			VDES_CONTAINER TT(tmpV.begin(), tmpV.begin() + MAX_GROUP);
			rNodeSetV.push_back(TT);
			tmpV.erase(tmpV.begin(), tmpV.begin() + MAX_GROUP);
			TT.clear();
		}

		//int kk = 0;
		for (nodeiter = order.begin(); nodeiter != order.end(); ++nodeiter){
			int a = out_degree(*nodeiter, circuit);
			int b = circuit[*nodeiter].visit_times;
			bool c = circuit[*nodeiter].ready;
			//bool d = BIT_TEST(circuit[*nodeiter].info->flag, OUTPUT_LOAD);
			if ((a <= b) && !c){
				circuit[*nodeiter].ready = true;
				tmpV.push_back(*nodeiter);
				//kk++;
			}
		}

		//cout << "Finish one rbatch! Having " << rNodeSetV.back().size() << " Nodes! "<< kk <<" new nodes.."<<endl;
	}

	/* ======topological order test============
	int kk = 0;
	for (int i = 0; i < rNodeSetV.size(); ++i){

		for (int j = 0; j < rNodeSetV[i].size(); ++j){
			circuit[rNodeSetV[i][j]].batch_no = i;
			kk++;
		}

	}
	cout << kk << endl;

	for (int i = 0; i < rNodeSetV.size(); ++i){
		for (int j = 0; j < rNodeSetV[i].size(); ++j){
			for (edgeiter = in_edges(rNodeSetV[i][j], circuit); edgeiter.first != edgeiter.second; ++edgeiter.first){
				faninv = source(*edgeiter.first, circuit);
				int a = (circuit[faninv].batch_no), b = (circuit[rNodeSetV[i][j]].batch_no);
				assert( a > b);
			}
		}
	
	}
	*/
}

VOID LRACC::scheduling(){

	VDES_CONTAINER::iterator nodeiter; //reversed traversal
	pair<TPLGY_OEITER, TPLGY_OEITER> edgeiter;
	VDES_CONTAINER tmpV;
	TPLGY_VDES  fanoutv;
	//---------find the primary input-----------
	for (nodeiter = order.begin(); nodeiter != order.end(); ++nodeiter){
		circuit[*nodeiter].visit_times = 0;
		if (BIT_TEST(circuit[*nodeiter].info->flag, INPUT_DRIVER)){
			circuit[*nodeiter].ready = true;
			tmpV.push_back(*nodeiter);
		}
		else{
			circuit[*nodeiter].ready = false;
		}
	}


	//----schedule-------
	while (!tmpV.empty()){

		if (tmpV.size() <= MAX_GROUP){

			for (nodeiter = tmpV.begin(); nodeiter != tmpV.end(); ++nodeiter)
			{
				for (edgeiter = out_edges(*nodeiter, circuit); edgeiter.first != edgeiter.second; ++edgeiter.first){

					fanoutv = target(*edgeiter.first, circuit);
					circuit[fanoutv].visit_times++;
				}
			}
			NodeSetV.push_back(tmpV);
			tmpV.clear();
		}
		else{
			for (nodeiter = tmpV.begin(); nodeiter < tmpV.begin() + MAX_GROUP; ++nodeiter)
			{
				for (edgeiter = out_edges(*nodeiter, circuit); edgeiter.first != edgeiter.second; ++edgeiter.first){
					fanoutv = target(*edgeiter.first, circuit);
					circuit[fanoutv].visit_times++;
				}
			}

			VDES_CONTAINER TT(tmpV.begin(), tmpV.begin() + MAX_GROUP);
			NodeSetV.push_back(TT);
			tmpV.erase(tmpV.begin(), tmpV.begin() + MAX_GROUP);
			TT.clear();
		}
		//int kk = 0;
		for (nodeiter = order.begin(); nodeiter != order.end(); ++nodeiter){
			int a = in_degree(*nodeiter, circuit);
			int b = circuit[*nodeiter].visit_times;
			bool c = circuit[*nodeiter].ready;
		//	bool d = BIT_TEST(circuit[*nodeiter].info->flag, OUTPUT_LOAD);
			if ((a <= b) && !c){
				circuit[*nodeiter].ready = true;
				tmpV.push_back(*nodeiter);
				//kk++;
			}
		}

	//	cout << "Finish one batch! Having " << NodeSetV.back().size() << " Nodes! "<< kk <<" new nodes.."<<endl;
	}

	/*
	// ======topological order test============
	int kk = 0;
	for (int i = 0; i < NodeSetV.size(); ++i){

	for (int j = 0; j < NodeSetV[i].size(); ++j){
	circuit[NodeSetV[i][j]].batch_no = i;
	kk++;
	}

	}
	cout << kk << endl;

	for (int i = 0; i < NodeSetV.size(); ++i){
	for (int j = 0; j < NodeSetV[i].size(); ++j){
	for (edgeiter = out_edges(NodeSetV[i][j], circuit); edgeiter.first != edgeiter.second; ++edgeiter.first){
	fanoutv = target(*edgeiter.first, circuit);
	int a = (circuit[fanoutv].batch_no), b = (circuit[NodeSetV[i][j]].batch_no);
	assert( a > b);
	}
	}

	}
	*/
}