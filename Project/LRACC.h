#ifndef LRACC_H
#define LRACC_H

#include <vector>
#include <string>
#include <fstream>
#include "Celllib.h"
#include "Circuit.h"
#include "werc_102_vlsi.h"
#include "timing_analysis.h"
#include "ocl_opt.h"
#include "TimingAnalysis.h"
#include <boost/thread.hpp>

class ocl_opt;
/*------------------------------------------------------------------------*/
/*                           Macro Definitions                            */
/*------------------------------------------------------------------------*/
#define     IS_SIMULTANEOUS     (NO)//YES = simultaneously Adaptivity +Varibility Only, NO = Sensitive method
//#define     IS_ISCAS85_TEST     (YES)//YES = ISCAS85; NO = ISPD Test
#define     IS_NO_ADAPTIVITY    (NO)//next time, to use environment variable
//#define     IS_OPTIMAL_CODE     (YES)//YES = optimal code, NO = original code

#define     SSTA_PARAM_IO1      (0)//Gate size
#define     SSTA_PARAM_IO2      (1)//Resistance
#define     SSTA_PARAM_IO3      (2)//Capacitance
#define     SSTA_PARAM_IO4      (3)//Vth voltage
#define     SSTA_PARAM_IO5      (4)//Leakage Power
#define     SSTA_PARAM_MAX      (5)

#define     PLAC_PARAM_IO1      (0)//X-position
#define     PLAC_PARAM_IO2      (1)//Y-position
#define     PLAC_PARAM_IO3      (2)//adaptivity_block

#define     SUBPP_AREA_SUM      (0)
#define     SUBPP_POWER_SUM     (1)
#define     SUBPP_LEAKP_SUM     (2)
#define     SUBPP_DELAY_SUM     (3)
#define     SUBPP_AREA_LM       (4)
#define     SUBPP_DELAY_LM      (5)
#define     SUBPP_MAX_DATA_NUM  (6)

#define     MAX_TUNING_UNIT     (4)//current tuning unit
#define     MAX_TUNING_FFF      (0xFFFF)//Invalid tuning unit

#define     ADAPTIVITY_0        (0)//no adaptivity circuit
#if IS_SIMULTANEOUS
#define     ADAPTIVITY_1        (1)//adaptivity circuit exist
#else
#define     ADAPTIVITY_1        (ADAPTIVITY_0)//adaptivity circuit doesn't exist
#endif
#define     MAX_AREA_LIMIT      (1000)//this might change for different circuit
#define     RESIST_PARAM_P      (0.2)//R=P*L/W

class CHANGES
{
public:
	CHANGES(TPLGY_VDES a, float b, int c, bool d = false) :gate(a), sensi(b), tileidx(c), block_change(d) {};
	TPLGY_VDES gate;
	bool block_change;
	int tileidx;
	float sensi;
};

typedef std::list<CHANGES> CHANGE_LIST;

class LRACC : public Timing_Analysis{
private:
	TPLGY circuit;
	//SOLUTIONS sltn;
	TILES tile;
    VDES_CONTAINER  order;//std::vector<TPLGY_VDES>
    //std::vector<float> lm;Lagrangian multiplier
    //UINT16DDD  iter_lagrangian;//Lagrangian loop number

	CELLLIB cell_lib;//only used as a cache before storing everything into TPLGY

    // The following are the variables for Timing Analysis
    // These are segregated according to interface between different codes
    // Inputs to Timing Analysis from Jiafan
    std::map<std::string, std::vector<float> > gate_sizes;  // The map key is "Gate/Net name" and value is the "Gate Size"
    //---- Each gates adaptive block . vector[0]--gate_size  vector[1]--R  vector[2]--C  vector[3]--Vth

		//double critical_path_delay;
		pair<double,double> power_mu_sigma;

    // Inputs from Timing Analysis to Jiafan
    std::map<std::string, std::vector<double> > pca_params; // Key = "Gate Name" , Value = "a STL Vector of pca parameters"
    std::map<std::string, std::vector<float> >  plac_info;  // Key = "Gate Name" , Value = "a triple containing x, y location and adp block index"
    std::map<int, double> probability;
		std::map<int, string> prob_type;
    std::map<int, bool> adpt_block_exist;//Key = Adaptivity block index, Value = YES/NO there needs an adaptivity circuit

    // Input from Timing Analysis to Edward
    map<std::string, pair<double,double> > mu_sigma;  // Key = "Gate Name" , Value = "a pair containing the mu and sigma"
    map<std::string, pair<double,double> > arrival_time_slack;  // Key = "Gate Name" , Value = "a pair containing the required arrival time and slack"
    float lm_a;
    ITER_JUDGE  trade_off_n;    /* Lagrangian iteration trade off */
    ITER_JUDGE  trade_off_r;    /* JRR iteration trade off */
    FINAL_SOL   final_result;
    std::ofstream	dest_file;		/* write each lagrangian info to the file */

private:
	
	boost::mutex i_mux;
	void rscheduling();
	void scheduling();
	vector<VDES_CONTAINER> NodeSetV;
	vector<VDES_CONTAINER> rNodeSetV;
	//int k_mux;
    void _pp_consistence_relax(BOOLEAN is_refine, int _index, int level);
	void _pp_consistence_restore(int _index, int level);
    bool _is_no_improvement(ITER_JUDGE &, vector<FLT32> &);
    vector<float> _get_inner_lagrangian_obj();

public:
    double critical_path_delay;
		int no_of_adaptive_blocks;
    UINT16D iter_lagrangian;
    //create graph, initialize nodes, initialize tiles, pca, store variations, probability of tuner states
    bool init_lracc_database(CHAR*);
#ifdef IS_ISCAS85_TEST
    bool load_graph_and_tiles(CHAR*);
#else
    bool load_graph_and_tiles(string&);
#endif
	bool load_celllib(CHAR *filename);
    bool load_supplement(CHAR *filename);
#ifdef IS_ISCAS85_TEST
    bool load_at_raq(CHAR *filename);
#else
    bool load_at_raq();
#endif
    bool load_tuners();
    FLT32 update_gates_info();//update gate size info for Rohit
	void load_plca_mark_edge();//order nodes topologically
#ifdef IS_OPTIMAL_CODE
    void opt_gate_solution();
#endif
    void LRACC_primal();
	void LRACC_dual();
	void LRACC_ssta();
	float calculate_leakage();
	//
	void SetAdaptivity(float);
	void sensitivity_analysis(CHANGE_LIST& M, float& gate_area, float& best_leakage, float& best_area, float& best_min_slack);
	void sensitivity_commit(CHANGE_LIST& M, std::vector<bool>& tilefixed, float& best_min_slack, int& result_flag, int& tmp_flag, float alpha);
	//
	void init_timing_analysis(string benchmark, string placement, double output_required_time);
	void check_gate_solution();//check gate_sizes and Vth solution is feasible
    bool record_gate_solution();//record and store the better solutions
    void print_lracc_result();
    void write_solution_to_file();
    FLT32 read_solution_by_file();
    //void print_curve_info();

    void calc_prob();
    //new functions
    void dsta();
		void get_power();
    void get_probability();
		void get_timing_yield();
    int incr_lag_num();

private:
	map<int, vector<float> > block_prob;//block_prob[i][FBB_GRADE_1]--probability of tile i in FBB level
    map<int, vector<float> > for_get_power;//set this probability before call get_power()
	float total_a;

private:
//////////////////////////cutting plane//////////////////////////////////
		vector<float > dual_hist_obj;//primal obj
  	vector< pair<float,map<string,pair<double,double> > > > dual_hist_subg;
		float total_p;
		ocl_opt* project;
///////////////////////////////////////////////////////////////////////////////////////////////

};

#endif

