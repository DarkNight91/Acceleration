#include "Circuit.h"
#include "LRACC.h"
#include "Celllib.h"
#include <math.h>

using namespace std;
using namespace boost;

COMPONENT::~COMPONENT()
{
    ;
}

COMPONENTOPT::~COMPONENTOPT()
{
    ;
}

FLT32   GATE::r_all(GATEOPT& o, UINT8D adpt_num)
{
    FLT32   r;
    r = o.resist[adpt_num];
	if (!(BIT_TEST(this->flag, INPUT_DRIVER)) && !(BIT_TEST(this->flag, OUTPUT_LOAD)) && !(BIT_TEST(this->flag, FF_VIRTUAL)))
    {
        GRADE   grad;
        grad.vbb = (ADAPTIVITY_0 == adpt_num)? ZERO_GRADE:FBB_GRADE_1;
        r += 3*r_sd(o, grad);
    }
    
    return r;
}

FLT32   GATE::r_sd(GATEOPT& o,GRADE g)
{
    ;//standard deviation of r is the sum of dr2dl*l_variation+dr2dw*w_variation+...
    return 0;
}


//derivative to process parameters
FLT32   GATE::dr2dl(GATEOPT& o,GRADE g)
{

    return 0;
}

FLT32   GATE::dr2dw(GATEOPT& o,GRADE g)
{
    return 0;
}


/* used to compute the power and area overhead of the adaptivity circuit */
pair<FLT32, FLT32>  NODEPROPERTY::get_gate_overhead(TPLGY_VDES faninv, TPLGY& g, TILES& tune_units)
{
    UINT16D  tile_index;
    TILE    *pTile;
    pair<FLT32,FLT32>   overhead;//First = area overhead, Second = power overhead
    UINT16D  gate_num;

    tile_index = this->info->tile_idx;
    if (MAX_TUNING_FFF == tile_index)
    {
        overhead.first = 0;
        overhead.second = 0;
    } 
    else
    {
        pTile = &(tune_units[tile_index]);//???????????????
        gate_num = pTile->components.size();
        H102_ASSERT(0 != gate_num);
        overhead.first = pTile->tuner.area()/gate_num;
        overhead.second = pTile->tuner.power()/gate_num;
    }
    
    return  overhead;
}


VOID NODEPROPERTY::list_and_merge(TPLGY_VDES fanoutv, TPLGY& g, list<FANIN_SUM> *fanin_set)
{
    pair<TPLGY_IEITER, TPLGY_IEITER> vi;
    TPLGY_VDES  faninv;
    UINT8D       adpt_num;
    COMPONENT_SOLUTION  *pCom_Sol[2];
    FANIN_SUM   fanin_tmp;

    fanin_set[ADAPTIVITY_0].clear();
#if IS_SIMULTANEOUS
    fanin_set[ADAPTIVITY_1].clear();
#endif

    vi = in_edges(fanoutv, g);
    for (; vi.first != vi.second; ++vi.first)
    {
        faninv = source(*vi.first, g);
        if (BIT_TEST(g[faninv].info->flag, INPUT_DRIVER)
            ||(g[faninv].info->tile_idx == g[fanoutv].info->tile_idx))
        {
            pCom_Sol[ADAPTIVITY_0] = &(*g[faninv].optimal_iter[ADAPTIVITY_0]);
#if IS_SIMULTANEOUS
            pCom_Sol[ADAPTIVITY_1] = &(*g[faninv].optimal_iter[ADAPTIVITY_1]);
#endif
        } 
        else
        {
            pCom_Sol[ADAPTIVITY_0] = &(*g[faninv].optimal_iter[ADAPTIVITY_0]);
#if IS_SIMULTANEOUS
            pCom_Sol[ADAPTIVITY_1] = &(*g[faninv].optimal_iter[ADAPTIVITY_1]);

            if (pCom_Sol[ADAPTIVITY_0]->propagates.obj_pru
                >= pCom_Sol[ADAPTIVITY_1]->propagates.obj_pru)
            {
                pCom_Sol[ADAPTIVITY_0] = pCom_Sol[ADAPTIVITY_1];
            } 
            else
            {
                pCom_Sol[ADAPTIVITY_1] = pCom_Sol[ADAPTIVITY_0];
            }
#endif
        }
        
        for (adpt_num = ADAPTIVITY_0; ;++adpt_num)
        {
            fanin_tmp.faninv = faninv;
            fanin_tmp.fanin_opt = pCom_Sol[adpt_num];
            fanin_set[adpt_num].push_back(fanin_tmp);

            if (ADAPTIVITY_1 == adpt_num)
                break;
        }
    }
}


/* used to merge the fanout solutions and pruning inferior combinations */
VOID NODEPROPERTY::merge_and_prune(TPLGY_VDES faninv, TPLGY& g,
    list<FANOUT_SUM*>& fanout_set, BOOLEAN is_refine)
{
    TPLGY_VDES  fanoutv;
    //list<FANOUT_SUM>    fanout_set[2];store the combination of the fanout gates
    list<FANOUT_SUM* >    fanout_tmp;//used as the temporary merged struct
    //FANOUT_SUM  vout_tmp;
    UINT8D       adpt_num;
    list<COMPONENT_SOLUTION>    *pfanout_can[2];//point to the candidates of the fanout NODECOMPONENT
    list<COMPONENT_SOLUTION>    candi_tmp[2];
    FLT32       cap;
    FLT32       obj_opt;
    FLT32       obj_pru;
    BOOLEAN     is_obj_valid;
    BOOLEAN     iter1_del;

    fanout_set.clear();
#if IS_SIMULTANEOUS
    fanout_set[ADAPTIVITY_1].clear();
#endif

    //enumerate and merge the fanout gates of current gate
    pair<TPLGY_OEITER, TPLGY_OEITER> vo = out_edges(faninv, g);

    for (; vo.first != vo.second; ++vo.first)
    {
		//int n = out_degree(faninv, g);
        fanoutv = target(*vo.first, g);
        is_obj_valid = BIT_TEST(g[*vo.first].edge_flag, REDUNDANT_EDGE)? NO:YES;

        //if the downstream gate is in the same block with fanin gate, then merge level by level
        if (BIT_TEST(g[faninv].info->flag, INPUT_DRIVER)
            || BIT_TEST(g[fanoutv].info->flag, OUTPUT_LOAD)
            ||(g[faninv].info->tile_idx == g[fanoutv].info->tile_idx))
        {
            if (is_refine && 1 < in_degree(fanoutv, g))
            {
                candi_tmp[ADAPTIVITY_0].clear();
                candi_tmp[ADAPTIVITY_0].push_back(*g[fanoutv].optimal_iter[ADAPTIVITY_0]);
                pfanout_can[ADAPTIVITY_0] = &(candi_tmp[ADAPTIVITY_0]);
#if IS_SIMULTANEOUS
                candi_tmp[ADAPTIVITY_1].clear();
                candi_tmp[ADAPTIVITY_1].push_back(*g[fanoutv].optimal_iter[ADAPTIVITY_1]);
                pfanout_can[ADAPTIVITY_1] = &(candi_tmp[ADAPTIVITY_1]);
#endif
            } 
            else
            {
                pfanout_can[ADAPTIVITY_0] = &(g[fanoutv].candidates[ADAPTIVITY_0]);
#if IS_SIMULTANEOUS
                pfanout_can[ADAPTIVITY_1] = &(g[fanoutv].candidates[ADAPTIVITY_1]);
#endif
            }
        }
        //if the downstream gate is in different blocks, then merge its adaptivity level first
        else
        {
            if (is_refine && 1 < in_degree(fanoutv, g))
            {
                candi_tmp[ADAPTIVITY_0].clear();
                candi_tmp[ADAPTIVITY_0].push_back(*g[fanoutv].optimal_iter[ADAPTIVITY_0]);
#if IS_SIMULTANEOUS
                candi_tmp[ADAPTIVITY_1].clear();
                candi_tmp[ADAPTIVITY_1].push_back(*g[fanoutv].optimal_iter[ADAPTIVITY_1]);
#endif
            } 
            else
            {
                candi_tmp[ADAPTIVITY_0] = g[fanoutv].candidates[ADAPTIVITY_0];
#if IS_SIMULTANEOUS
                candi_tmp[ADAPTIVITY_1] = g[fanoutv].candidates[ADAPTIVITY_1];
#endif
            }
#if IS_SIMULTANEOUS

            /***************************************************************************/
            /*  merge the solutions with/out adaptivity circuit, keep smaller obj_opt  */
            /***************************************************************************/
            list<COMPONENT_SOLUTION>::iterator  iter_1st = candi_tmp[ADAPTIVITY_0].begin();
            for (; iter_1st != candi_tmp[ADAPTIVITY_0].end();)
            {
                iter1_del = FALSE;
                list<COMPONENT_SOLUTION>::iterator  iter_2nd = candi_tmp[ADAPTIVITY_1].begin();
                for (; iter_2nd != candi_tmp[ADAPTIVITY_1].end();)
                {
                    if ((iter_2nd->propagates.c >= iter_1st->propagates.c)
                        && (iter_2nd->propagates.obj_opt >= iter_1st->propagates.obj_opt))
                    {
                        iter_2nd = candi_tmp[ADAPTIVITY_1].erase(iter_2nd);
                    }
                    else if ((iter_1st->propagates.c >= iter_2nd->propagates.c)
                        && (iter_1st->propagates.obj_opt >= iter_2nd->propagates.obj_opt))
                    {
                        iter1_del = TRUE;
                        break;
                    }
                    else
                    {
                        ++iter_2nd;
                    }
                }

                if (iter1_del)
                {iter_1st = candi_tmp[ADAPTIVITY_0].erase(iter_1st);}
                else
                {++iter_1st;}
            }

            //the results are all stored in ADAPTIVITY_0
            candi_tmp[ADAPTIVITY_0].splice(candi_tmp[ADAPTIVITY_0].end(), candi_tmp[ADAPTIVITY_1]);
            H102_ASSERT(0 == candi_tmp[ADAPTIVITY_1].size());
#endif
            pfanout_can[ADAPTIVITY_0] = &(candi_tmp[ADAPTIVITY_0]);
#if IS_SIMULTANEOUS
            pfanout_can[ADAPTIVITY_1] = pfanout_can[ADAPTIVITY_0];
#endif
        }

        //combine the candidate solutions in the same adaptivity level to the temporary structure
        for (adpt_num = ADAPTIVITY_0;; ++adpt_num)
        {
            fanout_tmp.clear();
			FANOUT_SUM vout_tmp;
            //create a empty fanout so we could standardize the process
            if (fanout_set.empty())
                fanout_set.push_back(&vout_tmp);

            //merge the new fanout with existing fanout_set and enlarge the fanout_tmp
            list<FANOUT_SUM*>::iterator  fanout_iter = fanout_set.begin();
            for (; fanout_iter != fanout_set.end(); ++fanout_iter)
            {
                list<COMPONENT_SOLUTION>::iterator sol_iter = pfanout_can[adpt_num]->begin();
                for (; sol_iter != pfanout_can[adpt_num]->end(); ++sol_iter)
                {
					FANOUT_SUM vout_tmp2;
					cap = (*fanout_iter)->c + sol_iter->propagates.c;
					obj_opt = (*fanout_iter)->obj_opt + is_obj_valid * sol_iter->propagates.obj_opt;
                    //obj_pru(i) = sum(obj_opt(j))+O(i) where j belongs to fanout(i)
                    //O(i)=Area(i)+Delay(i)+Power(i), obj_pru(i)=total cost of fanout cone
					obj_pru = (*fanout_iter)->obj_pru + sol_iter->propagates.obj_opt;
                    vout_tmp2.c = cap;
                    vout_tmp2.obj_opt = obj_opt;
                    vout_tmp2.obj_pru = obj_pru;
					vout_tmp2.fanout = (*fanout_iter)->fanout;
                    vout_tmp2.fanout.push_back(make_pair(sol_iter->opt, sol_iter->adpt_num));
                    fanout_tmp.push_back(&vout_tmp2);
                }
            }

            // prune the fanout combinations with larger value of sum(obj_opt(j)), 
            // which is now equal to obj_pru
            list<FANOUT_SUM*>::iterator  iter_1st = fanout_tmp.begin();
            for (; iter_1st != fanout_tmp.end();)
            {
                iter1_del = FALSE;
                list<FANOUT_SUM*>::iterator  iter_2nd = iter_1st;
                for (++iter_2nd; iter_2nd != fanout_tmp.end();)
                {
					if (((*iter_2nd)->c >= (*iter_1st)->c)
						&& ((*iter_2nd)->obj_pru >= (*iter_1st)->obj_pru))
                    {
                        iter_2nd = fanout_tmp.erase(iter_2nd);//iter_2nd = 
                    }
					else if (((*iter_1st)->c >= (*iter_2nd)->c)
						&& ((*iter_1st)->obj_pru >= (*iter_2nd)->obj_pru))
                    {
                        iter1_del = TRUE;
                        break;
                    }
                    else
                    {
                        ++iter_2nd;
                    }
                }
                if (iter1_del)
                {
                    iter_1st = fanout_tmp.erase(iter_1st);
                }
                else
                {
                    ++iter_1st;
                }
            }

            //update the fanout_set with fanout_tmp
            fanout_set = fanout_tmp;

            if (ADAPTIVITY_1 == adpt_num)
                break;
        }
    }
}


VOID NODEPROPERTY::solution_pruning(TPLGY_VDES faninv, TPLGY& g)
{
    UINT8D   adpt_num;
    BOOLEAN iter1_del;
    list<COMPONENT_SOLUTION>    *pfanout_can;

    for (adpt_num = ADAPTIVITY_0;; ++adpt_num)
    {
        pfanout_can = &(g[faninv].candidates[adpt_num]);

        list<COMPONENT_SOLUTION>::iterator  iter_1st = pfanout_can->begin();
        for (; iter_1st != pfanout_can->end();)
        {
            iter1_del = FALSE;
            list<COMPONENT_SOLUTION>::iterator  iter_2nd = iter_1st;
            for (++iter_2nd; iter_2nd != pfanout_can->end();)
            {
                if ((iter_2nd->propagates.c >= iter_1st->propagates.c)
                    && (iter_2nd->propagates.obj_pru >= iter_1st->propagates.obj_pru))
                {
                    iter_2nd = pfanout_can->erase(iter_2nd);
                }
                else if ((iter_1st->propagates.c >= iter_2nd->propagates.c)
                    && (iter_1st->propagates.obj_pru >= iter_2nd->propagates.obj_pru))
                {
                    iter1_del = TRUE;
                    break;
                }
                else
                {
                    ++iter_2nd;
                }
            }

            if (iter1_del)
            {
                iter_1st = pfanout_can->erase(iter_1st);
            }
            else
            {
                ++iter_1st;
            }
        }

        if (ADAPTIVITY_1 == adpt_num)
            break;
    }
}


/************************************************************************/
/*   To compute the gate delay (Rohit give me the PCA info as a whole   */
/************************************************************************/
float NODEPROPERTY::calc_gate_delay(GATEOPT& o, UINT8D adpt_num, float c, map<string, vector<double> > &pca, map<string, vector<double> >::iterator&  iter_pca1)
{
    FLT32   gate_delay = 0;
    FLT32   r;
 
	if (!BIT_TEST(this->info->flag, INPUT_DRIVER) && !(BIT_TEST(this->info->flag, OUTPUT_LOAD)) && !(BIT_TEST(this->info->flag, FF_VIRTUAL)))
    {
        H102_ASSERT(ADAPTIVITY_0 == adpt_num || (ADAPTIVITY_0+1) == adpt_num);
        r = (ADAPTIVITY_0 == adpt_num)?
            o.resist[ZERO_GRADE]:o.resist[FBB_GRADE_1];
        //o.resist[TUNING_GRADE_L]:o.resist[TUNING_GRADE_H];
        gate_delay = r*c;//mean delay

#ifdef  IS_OPTIMAL_CODE
        /*FLT32   three_sqrt_delta_L_delta_W = this->info->three_sqrt_delta_LW;

        //get the dR/dL()
        dRdL = r*0.15/3;// RESIST_PARAM_P * pow(w,-1);
        delta_L = this->info->deltaL * pow(dRdL, 2);

        //get the dR/dW()
        dRdW = r*0.08/3;//-RESIST_PARAM_P * pow(w,-2);
        delta_W = this->info->deltaW * pow(dRdW, 2);*/

        gate_delay += o.three_sqrt_delta_LW[adpt_num] * c;
#else
        UINT16D  num;
        UINT16D  middle = 0;
        FLT32   delta_L = 0;
        FLT32   delta_W = 0;
        FLT32   dRdL = 0;
        FLT32   dRdW = 0;

		map<string, vector<double> >::iterator  iter_pca = pca.find(this->info->name);
        H102_ASSERT(iter_pca != pca.end());

        FLT32   w;//get width  o.w*0.1 + 0.01
        GET_WIDTH(w, o.w);
        middle = iter_pca->second.size()/2;

        //get the dR/dL()
        dRdL = r*0.15/3;// RESIST_PARAM_P * pow(w,-1);
        for (num = 0; num < middle; ++num)
        {
            delta_L += pow(iter_pca->second[num],2);
        }
        delta_L *= pow(dRdL, 2);

        //get the dR/dW()
        dRdW = r*0.08/3;//-RESIST_PARAM_P * pow(w,-2);
        for (num = middle; num < 2*middle; ++num)
        {
            delta_W += pow(iter_pca->second[num],2);
        }
        delta_W *= pow(dRdW, 2);

        gate_delay += 3*sqrt(delta_L + delta_W)*c;
#endif
    }

    return  gate_delay;
}

// calculate the dynamic and leakage power of the gate
pair<FLT32, FLT32>  NODEPROPERTY::calc_gate_power(GATEOPT& o, UINT8D adpt_num, float cap, map<int, vector<float> > &prob)
{
    //pGate_opt->off_power[0] + 0.5*pGate_opt->cap[0];
    //FLT32   tmp_power;
    vector<FLT64>   local_prob(MAX_BB_GRADE,0);
    FLT32   leak_power = 0;
    FLT32   dyna_power = 0;
    //FLT32   gate_power = 0;

    if (BIT_TEST(this->info->flag, INPUT_DRIVER)
		|| ADAPTIVITY_0 == adpt_num || BIT_TEST(this->info->flag, FF_VIRTUAL))
    {
        dyna_power = MAX_DPOWER_ALPHA*cap*MAX_VDD_VOLTAGE*MAX_VDD_VOLTAGE*MAX_FREQUENCY;
        leak_power = o.off_power[ZERO_GRADE];
    } 
    else
    {
		
        map<int, vector<float> >::iterator  iter_prob = prob.find(this->info->tile_idx);
		if (iter_prob==prob.end())
			cout << "Name: " << this->info->name << " Flag: " << this->info->flag << " Tile idx: " << this->info->tile_idx << endl;
        H102_ASSERT(iter_prob != prob.end());

        /*if (!type[this->info->tile_idx].compare("0"))
        {local_prob[ZERO_GRADE] = iter_prob->second;
          H102_ASSERT(1==iter_prob->second);}
        else
        {
            if (!type[this->info->tile_idx].compare("FBB"))
            {local_prob[FBB_GRADE_1] = iter_prob->second;}
            else if (!type[this->info->tile_idx].compare("RBB"))
            {local_prob[RBB_GRADE_1] = iter_prob->second;}
            else
            {H102_ASSERT(FALSE);}
            //local_prob[TUNING_GRADE_H] = iter_prob->second;//supposing only the prob for changing L is given
            local_prob[ZERO_GRADE] = 1 - iter_prob->second;
        }*/

        UINT8D   grade_num;
        for (grade_num = ZERO_GRADE; grade_num < MAX_BB_GRADE; ++grade_num)
        {
            //leak_power += local_prob[grade_num]*o.off_power[grade_num];
            leak_power += iter_prob->second[grade_num]*o.off_power[grade_num];

            //gate_power += local_prob[grade_num] *tmp_power[grade_num];
            //if (TUNING_GRADE_H == grade_num)
            //    break;
        }
        dyna_power = MAX_DPOWER_ALPHA*cap*MAX_VDD_VOLTAGE*MAX_VDD_VOLTAGE*MAX_FREQUENCY;
    }

    return make_pair(dyna_power,leak_power);
}

