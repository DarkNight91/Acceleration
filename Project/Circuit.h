#ifndef CIRCUIT_H
#define CIRCUIT_H

#include <vector>
#include <string>
#include "boost/graph/topological_sort.hpp"
#include "boost/graph/depth_first_search.hpp"//dfs_visitor
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/graph_traits.hpp"
#include "boost/graph/visitors.hpp"
#include "Celllib.h"
#include <utility>

class COMPONENT;
class GATE;
class WIRE;

class FANOUT_SUM;
class FANIN_SUM;
class COMPONENT_SOLUTION;
class SOLUTIONS;

class PROPAGATES;
class NODEPROPERTY;
class EDGEPROPERTY;

class GRADE;
class STATE;
class TUNEUNIT;
class TILE;

//////////////////////////////////////////////////////////////////
//Description:
//	TPLGY is a graph that contains topology of the circuit,
//	type,placement,variation of each gate, tile each gate belongs to,
//	wire length, 
//	candidate solutions for each gate/wire, including whether to add 
//	adaptability,and the granularity to add adaptability
//Author: Hao He
//////////////////////////////////////////////////////////////////
typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::bidirectionalS,NODEPROPERTY,EDGEPROPERTY> TPLGY;
typedef boost::graph_traits<TPLGY>::vertex_descriptor TPLGY_VDES;
typedef boost::graph_traits<TPLGY>::vertex_iterator TPLGY_VITER;
typedef boost::graph_traits<TPLGY>::edge_descriptor TPLGY_EDES;
typedef boost::graph_traits<TPLGY>::edge_iterator TPLGY_EITER;
typedef boost::graph_traits<TPLGY>::out_edge_iterator TPLGY_OEITER;
typedef boost::graph_traits<TPLGY>::in_edge_iterator TPLGY_IEITER;
typedef boost::graph_traits<TPLGY>::adjacency_iterator TPLGY_AVITER;
typedef std::vector<TPLGY_VDES> VDES_CONTAINER;
typedef std::vector<TILE> TILES;


//////////////////////////////////////////////////////////////////
//Description:
//	parameters of tuning unit: tuning grades,
//	inverse of characteristic tuning function,
//	power consumption (mostly static), area
//Author: Hao He
//Note: tuning is performed once for a long period, 
//	characteristic tuning function is assumed to be staircase function
//	actual implementation may be more complex 
//////////////////////////////////////////////////////////////////
class GRADE
{
public:
	float vbb;
};

class TUNEUNIT
{
public:
	std::vector<GRADE> allstates;
  float   area_overhead;
  float   power_overhead;
	int     state;
    float area() {return area_overhead;};//area and power may differ w.r.t. configuration
    float power() {return power_overhead;};//4000000
  TUNEUNIT(): area_overhead(0), power_overhead(0){};
};

//////////////////////////////////////////////////////////////////
//Description:
//	tiles is a block within which variation is assumed the same
//	every tile has a potential tuner that can be dropped
//Author: Hao He
//////////////////////////////////////////////////////////////////
class TILE
{
public:
	float minX,minY,maxX,maxY;//range of positions
	std::vector<TPLGY_VDES> components;//idx to components this tile contains
	TUNEUNIT tuner;
	TPLGY_VDES tplgcl_first,tplgcl_last;//topological first and last gate in this tile
	bool is_adaptable;//temporarily store the solution while forward refinement
#if !IS_SIMULTANEOUS
	bool is_adaptable_dup;
	bool is_adaptable_mux;
#endif
	//calculate state probability by input range of different states and Gaussian distribution of path delay based on previous solution	
	//call calc_gate_delay/calc_wire_delay for each node and sum them up
	std::vector<float> calc_tuner_states_prob();
    FLT32   obj_value[2];//used to store the objective value of the whole adaptivity block
    TILE(): minX(0), minY(0), maxX(0), maxY(0) {obj_value[0] = 0; obj_value[1] = 0;}

    void  set_tune_overhead(UINT8D is_no_adaptivity)
    {
        {
            //set the area overhead according to the gate# in a certain block 20;
            this->tuner.area_overhead = this->components.size()*2+150;//0.13 + 87;
            this->tuner.power_overhead = this->components.size()*0.1+800;
        }
    }
};

//the following belongs to the private field of TPLGY,
//so making all their field public won't cause leakage
class COMPONENT
{
public:
#define INPUT_DRIVER    0x01        /* the input driver gate, store resistance */
#define OUTPUT_LOAD     0x02        /* the output load gate, store capacitance */
#define PRIMA_INGATE    0x04        /* primary input gates, special use */
#define FF_VIRTUAL      0x08        /* the virtual node created for Flip-Flop */

#ifdef IS_OPTIMAL_CODE
    //FLT32   deltaL;
    //FLT32   deltaW;
    FLT32   overhead_area;
    FLT32   overhead_power;
#endif

    std::string  cellInst;
    std::string  name;
    int     flag;
	int     tile_idx;
	COMPONENTTYPE type;//the index of the CELLIB
	float positionX;
	float positionY;
	virtual ~COMPONENT() = 0;//abstract class
    COMPONENT(): flag(0), tile_idx(0), positionX(0), positionY(0){};
};


// The DFS visitor is used to mark the edge which might cause the redundant objective computation
struct redundant_detector : public boost::dfs_visitor<>
{
    redundant_detector( std::vector<TPLGY_EDES> *has_cycle) 
        : _has_cycle(has_cycle) { }

    template <class TPLGY_EDES, class TPLGY>
    void forward_or_cross_edge(TPLGY_EDES e, TPLGY& g){
        //not allow to BIT_SET(g[e].edge_flag, REDUNDANT_EDGE) inside the visitor; back_edge
        _has_cycle->push_back(e);//
    }
protected:
    std::vector<TPLGY_EDES> *_has_cycle;
};

struct dfs_print_visitor : public boost::dfs_visitor<>
{
    template <class TPLGY_VDES, class TPLGY>
    void discover_vertex(TPLGY_VDES v, TPLGY& g) {
        std::cout << g[v].info->name<< " ";
    }
};

struct stack_gates_in_tile : public boost::dfs_visitor<>
{
    stack_gates_in_tile(TILES *tile_out)
        : _tile(tile_out){}

    template <class TPLGY_VDES, class TPLGY>
    void finish_vertex(TPLGY_VDES v, TPLGY& g) {
        //std::cout << g[v].info->name << " ";
        UINT16D  tile_idx;
        
        if (!BIT_TEST(g[v].info->flag, INPUT_DRIVER) && !BIT_TEST(g[v].info->flag, OUTPUT_LOAD) && !BIT_TEST(g[v].info->flag, FF_VIRTUAL))
        {
            tile_idx = g[v].info->tile_idx;
            (*_tile)[tile_idx].components.push_back(v);
        }
    }
protected:
    TILES *_tile;
};


#define gate_cast(x) (dynamic_cast<GATE *>(x))
class GATE: public COMPONENT
{
public:
    float r_all(GATEOPT& o, UINT8D adpt_num);
	float r_sd(GATEOPT& o,GRADE g);//standard deviation of r is the sum of dr2dl*l_variation+dr2dw*w_variation+...
    float dr2dl(GATEOPT& o,GRADE g);
    float dr2dw(GATEOPT& o,GRADE g);
	//derivative to process parameters
    float c_sd(GATEOPT o,GRADE g);
	float dc2dl(GATEOPT o,GRADE g);
	float dc2dw(GATEOPT o,GRADE g);
	//process variations in terms of pca coef
	std::vector<float> l_variation;
	std::vector<float> w_variation;
	//thermal and aging variation can be represented the same way
	//by derivative to certain parameters * pca coefs
};

#define wire_cast(x) (dynamic_cast<WIRE *>(x))
class WIRE: public COMPONENT
{
public:
	float wire_length;
};
 
class PROPAGATES
{
public:
	float c;//capacitance seen from fanins of this component	
	float obj_opt;//objective value for downstream circuits
    float obj_pru;//objective value for pruning downstream
    //float at;primary input
    //float raq;primary output
    PROPAGATES() : c(0), obj_opt(0), obj_pru(0){};//, at(0), raq(0)
};

typedef std::vector<int> DIGIT_TAG;
class NODEPROPERTY
{
public:
	//propagate solutions from fanouts, split and merge solutions when enter and leave a tile, operate within each group when revisit a tile
	//merge_and_prune calls merge_one, and it dynamically explore whether a node is a gate or a wire before calling merge_one
	//merge_and_prune updates both candidates of this node and SOLUTIONS
    void list_and_merge(TPLGY_VDES, TPLGY&, std::list<FANIN_SUM> *);
    std::pair<float, float> get_gate_overhead(TPLGY_VDES, TPLGY& , TILES& );
    void merge_and_prune(TPLGY_VDES, TPLGY&, std::list<FANOUT_SUM*>& fanout_set, BOOLEAN);
    void solution_pruning(TPLGY_VDES, TPLGY&);
	//propagate solution from fanins, either inherit a solution if it's single fanin node, or enumerate all possible solutions if it's multiple fanin node
	//once a solution is selected after backwards refinement, each tile gets a fixed adaptability and is stored in tile.is_adaptable
	//inherit_or_enumerate calls fanin_obj of all its fanins to calculate fanins' contribution to the cumulative objective, it dynamically dispatches
	void inherit_or_enumerate();

	//information of this component, need to be casted to GATE or WIRE according to info.type
	COMPONENT   *info;
	int visit_times;
	int batch_no;
	bool ready;
    float   lm;//Lagrangian multiplier for delay
	bool touched;
		
    std::list<COMPONENT_SOLUTION>::iterator optimal_iter[2];
    //FLT32   foward_at;

    std::list<COMPONENT_SOLUTION>   candidates[2];//for ADAPTIVITY_0 and ADAPTIVITY_1
	//record solution candidates(propagated from downstream in backward refinement) of multi-fanin node for forward refinement
	//in forward refinement of multi-fanin, the group whose digit_tag matches the settled tile adaptability config is kept and enumerated
	//then only a best solution is kept in candidates
	//std::list<std::pair<std::list<COMPONENT_SOLUTION>, DIGIT_TAG> > candidates;

    GRADE   tile_grade();
	void tile_split(int digit);//called at the beginning of merge_and_prune when entering a tile
	void tile_merge_and_prune(int digit);// called in the end of merge_and_prune when leaving a tile(visit the last node)
	//calculations for gate
    //delay depends on variation parameters, downstream capacitance and current option
	float calc_gate_delay(GATEOPT&, UINT8D, float, std::map<std::string, std::vector<double> > &, std::map<std::string, std::vector<double> >::iterator&);
	std::pair<float, float> calc_gate_power(GATEOPT&, UINT8D, float, std::map<int, std::vector<float> > &);

	float calc_gate_expected_power(GATEOPT o,TUNEUNIT t,PROPAGATES p);//depends on probability of tuner states and power of current gate in each state
	float calc_gate_timing_constraints(GATEOPT o,TUNEUNIT t,PROPAGATES p);
	float calc_gate_area_constraints(GATEOPT o,TUNEUNIT t,PROPAGATES p);
	//overloading for gate called by merge_and_prune/inherit_or_enumerate
	PROPAGATES merge_one(GATEOPT o,std::vector<PROPAGATES> fanout_solutions);//update fanin_c,obj based on fanout solutions and opts,is_adaptable
	float fanin_obj(GATEOPT o,PROPAGATES p);//calculate objective of an assigned node given one fanout, used for forward refinement, called by fanout nodes
	//overloading for wire called by merge_and_prune/inherit_or_enumerate
	PROPAGATES merge_one(WIREOPT o,std::vector<PROPAGATES> fanout_solutions);
	float fanin_obj(WIREOPT o,PROPAGATES p); 

    NODEPROPERTY() : lm(0) {}
};

class EDGEPROPERTY
{
public:
#define INPUT_EDGE      0x01        /* connect to the input driver gate */
#define OUTPUT_EDGE     0x02        /* connect to the output load gate */
#define REDUNDANT_EDGE  0x04        /* connect to the node which has been computed the objective */

    std::string  name;
    int edge_flag;
    float lm;
    EDGEPROPERTY() : edge_flag(0),lm(10) {}
};

//////////////////////////////////////////////////////////////////
//Description:
//	FANOUT_SUM is used to store the fanout gates' combinational
//	solution, and the corresponding C and objective value 
//Author: Jiafan Wang
//////////////////////////////////////////////////////////////////
class FANOUT_SUM
{
public:
    std::vector<std::pair<COMPONENTOPT *, bool> > fanout;
    float   c;
    float   obj_opt;
    float   obj_pru;
    VOID clear() {fanout.clear(); c = 0; obj_opt = 0; obj_pru = 0;};
    FANOUT_SUM(): c(0),obj_opt(0),obj_pru(0){};
};

class FANIN_SUM
{
public:
    TPLGY_VDES              faninv;
    COMPONENT_SOLUTION*     fanin_opt;
    INT16D       tile_idx;
    FLT32       fanin_pad;
    FANIN_SUM(): faninv(0), fanin_opt(NULL), tile_idx(0), fanin_pad(0){};
};

//////////////////////////////////////////////////////////////////
//Description:
//	SOLUTIONS is a data structure that stores the propagated,
//	conflicted, and grouped solutions 
//Author: Hao He
//////////////////////////////////////////////////////////////////
class COMPONENT_SOLUTION//solution for a single component
{
public:
	//gate solution is a tuple of gate size, threshold voltage and adaptability
	//wire solution is a wire width
    bool    adpt_num;
	COMPONENTOPT * opt;//pointer to celllib
	//results of this solution, including the 
	PROPAGATES  propagates;
    //the fanout gates' combinational solution
    FANOUT_SUM  fanouts;
    COMPONENT_SOLUTION():fanouts(){};
};

/* this structure is used to improve the performance */
class ITER_JUDGE
{
public:
    FLT32   threshold;          /* IJRR_TRADEOFF_N/IJRR_TRADEOFF_R */
    std::list<std::vector<FLT32> >   object;/* [MAX_LATEST_ITERATION];the first object is the newest */

    /* the sequence of the vector<FLT32> is area, power and delay */
    VOID    initial(FLT32 in_threshold)
    {
        UINT8D   num = 0;
        std::vector<FLT32>   area_p_d(8, 0);//SUBPP_MAX_DATA_NUM should < 8
        //area_p_d.push_back(0);
        //area_p_d.push_back(0);
        //area_p_d.push_back(0);
        //area_p_d.push_back(0);
        //area_p_d.push_back(0);

        object.clear();
        for (num = 0; num < MAX_LATEST_ITERATION; ++num)
        {
            object.push_back(area_p_d);
        }
        threshold = in_threshold;
    }
};

#define GET_WIDTH(_GATE_WIDTH, _GATE_SIZE)\
    do{\
    if(_GATE_SIZE < 4)\
    (_GATE_WIDTH) = (_GATE_SIZE)+1;\
    else if(_GATE_SIZE < 7)\
    (_GATE_WIDTH) = (_GATE_SIZE)*2-2;\
    else if (_GATE_SIZE < 10)\
    (_GATE_WIDTH) = (_GATE_SIZE-6)*20;\
    else\
    H102_ASSERT(FALSE);\
    }\
    while(0)

typedef struct  GATE_FINAL
{
    TPLGY_VDES  fgate;
    GATEOPT     *pGate;
    BOOLEAN     is_adpt;
};

class FINAL_SOL
{
public:
#define AREA_INFEASIBLE     0x01
#define TIME_INFEASIBLE     0x02
#define INIT_INFEASIBLE     0x04
    std::map<int, std::vector<float> >  tile_prob;

    UINT16D  iter_lm;
    int     solution_flag;
    std::list<GATE_FINAL>  fopt_set;
    std::map<int, bool>   tile_adpt;
    FLT32   leak_sum;
    FLT32   area_sum;
    FLT32   power_sum;
    FLT32   min_slack;
    FLT32   max_slack;
    FLT32   adpt_area;
    FLT32   adpt_power;
    FINAL_SOL(): iter_lm(0), solution_flag(INIT_INFEASIBLE), area_sum(0), power_sum(0), leak_sum(0), min_slack(0), max_slack(0), adpt_area(0), adpt_power(0) {};
};

#if 0
//conflicted solution of whole circuit
//because the solution may be conflicted for the first backward refinement, component_solution is stored as a vector
typedef std::map<TPLGY_VDES, std::vector<COMPONENT_SOLUTION> > CIRCUIT_SOLUTION;
typedef std::vector<int> DIGIT_TAG;
//vectorization of multi-dim array, each element correspond is a solution set of a given tile adaptability config and a tag encoding the tile adaptability config
//a config is encoded by a binary number whose no. of digits==no. of tiles, which digit a tile is associated is recorded in digit2tile 
//example: solutions[i].second = 010 means no adaptability for the 1st and 3rd digit tile
class SOLUTIONS
{
public:
	std::list< std::pair< std::list< CIRCUIT_SOLUTION >, DIGIT_TAG > > data;
	DIGIT_TAG digit2tile;
};
#endif



#endif
