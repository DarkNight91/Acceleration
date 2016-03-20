#include "timing_analysis.h"
#include "LRACC.h"
#include "Circuit.h"
#include "iostream"
#include "Celllib.h"
#include "parser_helper.h"
#include "werc_102_vlsi.h"

FLT32   Time_Limit = 0;
UINT32D	Area_Limit = 0;
UINT16D  gGrid_num = 0;  


int main(int args, char** argv)
{
	argv[1] = "fft";
	argv[2] = "800";
	argv[3] = "20000";
	string benchfile("E:/Icanmakeit/new/Project2/Project2/benchmarks/");
	string placefile("E:/Icanmakeit/new/Project2/Project2/benchmarks/");
	string timingfile("E:/Icanmakeit/new/Project2/Project2/benchmarks/");

    benchfile += argv[1];
#ifdef IS_ISCAS85_TEST
    benchfile += ".bench";
#else
    benchfile += ".v";
#endif
    placefile += argv[1];
    placefile += ".pl";
    timingfile += argv[1];
    timingfile += ".at";
    //Area_Limit = Area_Limit>MAX_AREA_LIMIT?Area_Limit:MAX_AREA_LIMIT;

    bool    rtn = FALSE;
    LRACC project;
    cout << "Just Compiling Nothing else " << endl;
    rtn = project.init_lracc_database(argv[1]);
#ifdef IS_ISCAS85_TEST
    rtn = project.load_graph_and_tiles(&benchfile.at(0));
#else
    rtn = project.load_graph_and_tiles(benchfile);
#endif

	rtn = project.load_celllib("E:/Icanmakeit/new/Project2/Project2/benchmarks/contest.lib");
    //rtn = project.load_supplement("supplement2.lib");
#ifdef IS_ISCAS85_TEST
    rtn = project.load_at_raq(&timingfile.at(0));//"benchmarks/c432.at");
#else
    rtn = project.load_at_raq();//we may need the timing info for ISPD, add here
#endif

    Area_Limit = atoi(argv[3]);
    Time_Limit = atof(argv[2]);
    FLT32  Area_init = project.update_gates_info();//update the gate size info
    FLT32  Timing_Limit = atof(argv[2]);//project.critical_path_delay - 3;//atof(argv[2]);
    project.init_timing_analysis(&benchfile.at(0), &placefile.at(0), Timing_Limit);//"benchmarks/c432.bench","benchmarks/c432.pl", 4);
    rtn = project.load_tuners();
    project.load_plca_mark_edge();
#ifdef  IS_OPTIMAL_CODE
    project.opt_gate_solution();
#endif

    cout<<"Timing_Limit from input = "<< Timing_Limit <<endl;
    cout<<"Area_Limit = "<< Area_Limit <<"\t Area_init = "<<Area_init<<endl;

    clock_t t1,t2,t3,t4;
    
    while(project.incr_lag_num() <= MAX_LAGRANGIAN_NUM)//or reach the boundary
    { 
		
        project.calc_prob();
        project.LRACC_dual();
        //project.get_probability();
       	//bool best_seen = project.record_gate_solution();
#if !IS_SIMULTANEOUS
        if(project.iter_lagrangian%4==3)
				{
					//t3 = clock();
					project.SetAdaptivity(0.05);
					//t4 = clock();
					//float d = (float)t4 - (float)t3;
					//cout << "Time for one SENSITIVITY!! = " << d / CLOCKS_PER_SEC << endl;

				}
#endif
		t1 = clock();
	    project.LRACC_primal();
		t2 = clock();
		float diff((float)t2 - (float)t1);
		cout << "Time for one iter = " << diff / CLOCKS_PER_SEC << endl;
		getchar();
        project.update_gates_info();//update gate size info
        //project.print_curve_info();

		//t3 = clock();
		project.LRACC_ssta();
		//t4 = clock();
		
	  }; 
		////////////////////calculate statistics using best solution/////////////////////////
  

    //project.calc_prob();
   /// float Timing_Limit = atof(argv[2]);//project.critical_path_delay - 3;//atof(argv[2]);
   /// project.init_timing_analysis(&benchfile.at(0), &placefile.at(0), Timing_Limit);//"benchmarks/c432.bench","benchmarks/c432.pl", 4);
   // project.dsta();
	getchar();
    project.LRACC_ssta();

	//project.get_timing_yield();
    cout << endl << endl << "DONE WITH Timing Yield " << endl;
    //cin >> wait;
	project.get_probability();
#ifndef NDEBUG
   // cout<<"Time Elapsed(s) = "<<diff/CLOCKS_PER_SEC<<endl;
    //project.print_lracc_result();
#endif
	getchar();
    return 0;
}

