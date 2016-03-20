#ifndef CELLLIB_H
#define CELLLIB_H

#include <boost/config.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <iterator>
#include "werc_102_vlsi.h"
#include <math.h>

//#include "Circuit.h"

//////////////////////////////////////////////////////////////////
//Description:
//	cell library includes options like sizes, threshold voltages
//	for each type of gates and corresponding eletric parameters
//	such as input capacitance, output resistance for each option
//	and each tuning grades, area
//Author: Hao He
//Change: Jiafan Wang
//////////////////////////////////////////////////////////////////

/*------------------------------------------------------------------------*/
/*                           Macro Definitions                            */
/*------------------------------------------------------------------------*/
#define INPUT_SLEW_INDEX        (2)
#define MAX_FOOTPRINT_LEN       (4)
#define MAX_GATE_NAME           (8)
#define MAX_GATE_SIZE_TYPE      (PN_SIZE_MAX)
#define MAX_VT_LEVEL_TYPE       (VT_LEVEL_MAX)
#define MAX_SOLUTION_NUM        (MAX_GATE_SIZE_TYPE*MAX_VT_LEVEL_TYPE)
#define MAX_SOL_VALID_NUM       (MAX_SOLUTION_NUM)//MAX_SOLUTION_NUM-3
#define MAX_OBJECTIVE           (3.40E+38)//4503599627370496

#define MAX_DPOWER_ALPHA        (0.5)
#define MAX_VDD_VOLTAGE         (1.2)//0.7v
#define MAX_FREQUENCY           (2)//2GHz
#define MAX_LAGRANGIAN_NUM      (8)//maximum outer Lagrangian iteration loop
#define MAX_LATEST_ITERATION    (3)//maximum number of record
#define EULER_NUMBER            (2.71828)

typedef enum tagIJRR_PN_SIZE
{
    PN_SIZE_1X = 0,
    PN_SIZE_2X,
    PN_SIZE_3X,
    PN_SIZE_4X,
    PN_SIZE_6X,
    PN_SIZE_8X,
    PN_SIZE_10X,
    PN_SIZE_20X,
    PN_SIZE_40X,
    PN_SIZE_80X,
    PN_SIZE_MAX
}IJRR_PN_SIZE;

typedef enum tagIJRR_VT_LEVEL
{
    VT_LEVEL_S = 0,
    VT_LEVEL_M,
    VT_LEVEL_F,
    VT_LEVEL_MAX
}IJRR_VT_LEVEL;

typedef enum{
    INV = 0,
    NAND2,
    NAND3,
    NAND4,
    NAND5,
    NAND6,
    NAND7,
    NAND8,
    NAND9,
    NOR2,
    NOR3,
    NOR4,
    NOR5,
    NOR6,
    NOR7,
    NOR8,
    NOR9,
    XOR2,
    OA12,
    OA22,
    AO12,
    AO22,
    FFMS,
    //AND2,
    //AND3,
    //OR2,
    //OR3,
    IN1,
    OUT1,
    WIRE,
    MAX_CP_TYPE
}COMPONENTTYPE;

typedef struct tagIJRR_NAME_MAP
{
    std::string  netlist_name;//[MAX_GATE_NAME];
    std::string  stdlib_name;//[MAX_GATE_NAME];
}IJRR_NAME_MAP;

typedef std::map<std::string, COMPONENTTYPE>  GMAP_INDEX;

typedef enum{
    ZERO_GRADE = 0,
    FBB_GRADE_1,
    RBB_GRADE_1,
    MAX_BB_GRADE
}TUNE_GRADE;

#define GET_SHIFT_UNIT(_tune_index, _shift_unit)\
    do{\
    if(ZERO_GRADE == (_tune_index))\
    (_shift_unit) = 0;\
    else if(FBB_GRADE_1 == (_tune_index))\
    (_shift_unit) = 1;\
    else if (RBB_GRADE_1 == (_tune_index))\
    (_shift_unit) = -1;\
    else\
    H102_ASSERT(FALSE);\
    }\
    while(0)

/*------------------------------------------------------------------------*/
/*                          Class  Definitions                            */
/*------------------------------------------------------------------------*/
class GRADE;
class COMPONENTOPT
{
public:
	virtual ~COMPONENTOPT() = 0;//abstract class
};
class GATEOPT: public COMPONENTOPT
{
public:
#ifdef  IS_OPTIMAL_CODE
    FLT32   three_sqrt_delta_LW[2];
#endif

	//for gate
	float l;
	IJRR_PN_SIZE    w;
	IJRR_VT_LEVEL   vth;
	float area;//area is not relevant grade
	//adaptability, is not contained in cell lib
	bool is_adaptable;
    std::vector<float>  resist;
    std::vector<float>  off_power;
    float               cap;
    float               offset;

    float leak_power(GRADE g);
    //return r,c of current gate option given a tuning grade
    //may be calculated by a function, or from a LUT
    float r(GRADE g);//nominal value
    float c(GRADE g);
    float off(GRADE g);

    GATEOPT (): l(1), w(PN_SIZE_1X), vth(VT_LEVEL_S), area(0), is_adaptable(false){};
    ~GATEOPT(VOID){};
    VOID set_solution(UINT8D sol_index) 
    {
        w = (IJRR_PN_SIZE)(sol_index%PN_SIZE_MAX);//length is set to be 1
        vth = (IJRR_VT_LEVEL)(sol_index/PN_SIZE_MAX);
    }

    FLT32 get_leakpower_bb(INT8D shift)
    {
        //in 90nm tech, original power/(e^(0.05 or -0.05/nVt)); n = 1.5 what i use and Vt is 0.026
        FLT32   pbb = 0;
        FLT32   n = 1.5;
        FLT32   vT = 0.026;
        pbb = off_power[0]/(pow(EULER_NUMBER, -0.05*shift/(n*vT)));
        return pbb;
    }

    FLT32 get_resist_bb(INT8D shift)
    {
        FLT32   rbb = 0;
        rbb = resist[0]/(pow(1+0.05*shift, 2));
        return rbb;
    }
};

class WIREOPT: public COMPONENTOPT
{
public:
	float width;
	//return r,c of the wire given width
	float r(float l);
	float c(float l);
};

//1st dim: type (NAND,NOT,etc.)
//2nd dim: opts (1st opt, 2nd opt,...)
typedef std::vector< std::vector<COMPONENTOPT *> > CELLLIB;
#define gateopt_cast(x) (dynamic_cast<GATEOPT *>(x))
#define wireopt_cast(x) (dynamic_cast<WIREOPT *>(x))


/*------------------------------------------------------------------------*/
/*                          Global variable                               */
/*------------------------------------------------------------------------*/
extern GMAP_INDEX gGate_Index;
extern IJRR_NAME_MAP gGate_map[MAX_CP_TYPE];
extern UINT32D	Area_Limit;
extern UINT16D gGrid_num;

#endif

