#include "LRACC.h"
#include <algorithm>
#include <numeric>
#include "Circuit.h"
#include <boost/thread.hpp>

int action = -1;
class CompCHANGES
{
public:
	bool operator()(CHANGES a,CHANGES b) {return a.sensi<b.sensi;}
} cmpChanges;

bool compareMapOn2nd(std::map<std::string,std::pair<double,double> >::value_type a,std::map<std::string,std::pair<double,double> >::value_type b)
{
	return a.second.second<b.second.second;
}
float LRACC::calculate_leakage()
{
	float leakage = 0;
	for(TPLGY_VITER viter=boost::vertices(circuit).first;viter!=boost::vertices(circuit).second;viter++)
	{
		if(BIT_TEST(circuit[*viter].info->flag,INPUT_DRIVER)||BIT_TEST(circuit[*viter].info->flag,OUTPUT_LOAD)||(circuit[*viter].info->type==FFMS)) continue;
		GATEOPT *pGate = gateopt_cast(circuit[*viter].optimal_iter[0]->opt);
		int blkidx = (circuit[*viter].info)->tile_idx;
		vector<float>   prob(MAX_BB_GRADE,0);
    if (!tile[blkidx].is_adaptable_dup)
   	{
			prob[ZERO_GRADE] = 1;
			prob[FBB_GRADE_1] = 0;
			prob[RBB_GRADE_1] = 0;
  	} 
   	else
    {
			prob = block_prob[blkidx];
		}
		leakage += std::inner_product(prob.begin(),prob.end(),pGate->off_power.begin(),0.0);//expected leakage power
	}
	for(int i=0;i<tile.size();++i)
	{
		if(tile[i].is_adaptable) leakage+=tile[i].tuner.power_overhead;
	}
	return leakage;
}


void LRACC::SetAdaptivity(float alpha)
{
	if(iter_lagrangian==1)
	{
		for(int i=0;i<tile.size();i++)
		{
			tile[i].is_adaptable = false;
		}
		return;
	}

#if !IS_NO_ADAPTIVITY
	//initialize
	for(int i=0;i<tile.size();i++)
	{
		tile[i].is_adaptable_dup = tile[i].is_adaptable;
	}
		/////////////////////////restore from final_result////////////////////
    list<GATE_FINAL>::iterator  iter_final = final_result.fopt_set.begin();
    for (; iter_final != final_result.fopt_set.end(); ++iter_final)
    {
        if (BIT_TEST(circuit[iter_final->fgate].info->flag, INPUT_DRIVER)
            || BIT_TEST(circuit[iter_final->fgate].info->flag, OUTPUT_LOAD)) continue;
        circuit[iter_final->fgate].optimal_iter[0]->opt = iter_final->pGate;
    }
	//////////////////////////////////////////////////////////////////////
		float gate_area = 0;
		
		for(TPLGY_VITER viter=boost::vertices(circuit).first;viter!=boost::vertices(circuit).second;viter++)
		{
			if(BIT_TEST(circuit[*viter].info->flag,INPUT_DRIVER)||BIT_TEST(circuit[*viter].info->flag,OUTPUT_LOAD)||(circuit[*viter].info->type==FFMS)) continue;
			std::string nm=circuit[*viter].info->name;
			GATEOPT *pGate = gateopt_cast(circuit[*viter].optimal_iter[0]->opt);
			std::vector<float> tmp;
			tmp.push_back(pGate->w);
			tmp.push_back(pGate->resist[0]);
			tmp.push_back(pGate->cap);
			tmp.push_back(pGate->vth);
			gate_sizes[nm] = tmp;
			gate_area += pGate->area;
		}
		
	float best_leakage = 1e100;
	float best_min_slack = -1e100;
	float best_area = 1e100;
	int result_flag = 0;
	std::vector<bool> tilefixed = std::vector<bool>(tile.size(),false);

  //calculate probability
  get_probability();
  for(int i=0;i<tile.size();++i)
  {
     std::cout<<i<<":ZERO GRADE has "<<block_prob[i][0];
     std::cout<<" ,FBB GRADE has "<<block_prob[i][1];
     std::cout<<" ,RBB GRADE has "<<block_prob[i][2]<<endl;
  }

  int tmp_flag = 0;
  while (1)
  {
	  CHANGE_LIST M;
	  bool p = false;
	  if (p){
		  boost::thread analysis(&LRACC::sensitivity_analysis, this, M, gate_area, best_leakage, best_area, best_min_slack);
		  boost::thread commit(&LRACC::sensitivity_commit, this, M, tilefixed, best_min_slack, result_flag, tmp_flag, alpha);

		  std::cout << "Joining threads....." << endl;

		  analysis.join();
		  commit.join();

		  if (action == 0||action==-1) continue;
		  else goto update_final;
	  }
	  else{
		  LRACC::sensitivity_analysis(M, gate_area, best_leakage, best_area, best_min_slack);
		  if (action==1) goto update_final;
		  LRACC::sensitivity_commit(M, tilefixed, best_min_slack, result_flag, tmp_flag, alpha);
		  if (action == 1) goto update_final;
	  }
  }
update_final:
	/////renew lm
	for(TPLGY_EITER eiter = edges(circuit).first;eiter!=edges(circuit).second;++eiter)
	{
		circuit[*eiter].lm*=2;
	}
	final_result.adpt_power = 0;
	for(int i=0;i<tile.size();++i)
	{
		if(tile[i].is_adaptable) final_result.adpt_power+=tile[i].tuner.power_overhead;
		final_result.tile_adpt[i] = tile[i].is_adaptable;
	}
	float dynamic = final_result.power_sum - final_result.leak_sum;
	final_result.leak_sum = best_leakage;
	final_result.power_sum = best_leakage + dynamic;
	final_result.area_sum = best_area;
	final_result.min_slack = best_min_slack;
	final_result.adpt_area = best_area - gate_area;
	final_result.solution_flag = result_flag;
#endif
}


void LRACC::sensitivity_analysis(CHANGE_LIST& M, float& gate_area, float& best_leakage, float& best_area, float& best_min_slack){
	
	std::cout << "Start analysis...." << endl;

	float area = gate_area;
	//sta to get slack and compute area
	for (int num = 0; num < tile.size(); ++num)
	{
		adpt_block_exist[num] = tile[num].is_adaptable_dup;
		cout << "trying to change tile " << num << " to " << adpt_block_exist[num] << "\n";
		if (adpt_block_exist[num])
			area += tile[num].tuner.area();
	}

	LRACC_ssta();
	float leakage = calculate_leakage();
	float min_slack = (*std::min_element(arrival_time_slack.begin(), arrival_time_slack.end(), compareMapOn2nd)).second.second;
	int tmp_flag = 0;
	if (min_slack<0) BIT_SET(tmp_flag, TIME_INFEASIBLE);
	if (area>Area_Limit) BIT_SET(tmp_flag, AREA_INFEASIBLE);
	std::cout << "leakage = " << leakage << ",min slack = " << min_slack << ",area = " << area << std::endl;
	std::cout << "best leakage = " << best_leakage << ",best min slack = " << best_min_slack << ",best area = " << best_area << std::endl;
	bool try_other_tile = false;
	if (best_min_slack<0)
	{
		if (min_slack >= best_min_slack)
		{
			best_min_slack = min_slack;
			best_area = area;
			best_leakage = leakage;
		}
		else// if(min_slack==best_min_slack)
		{
			cout << "wrong tile sensitivity,try others\n";
			try_other_tile = true;
		}
		/*else
		{
		cout<<"slack worse than best\n";
		goto update_final;
		}*/
	}
	else
	{
		if (min_slack<0)
		{
			cout << "slack worse than best\n";
			action = 1; //final
			std::cout << "Finishing analysis...." << endl;
			return;
			//goto update_final;
		}
		if (best_area>Area_Limit)
		{
			if (area<best_area)
			{
				best_min_slack = min_slack;
				best_area = area;
				best_leakage = leakage;
			}
			else
			{
				cout << "slack ok, area worse than best\n";
				action = 1; //final
				std::cout << "Finishing analysis...." << endl;
				return;
				//goto update_final;
			}
		}
		else
		{
			if (area>Area_Limit)
			{
				cout << "slack ok, area worse than best\n";
				action = 1; //final
				std::cout << "Finishing analysis...." << endl;
				return;
				//goto update_final;
			}
			if (best_leakage>leakage)
			{
				best_min_slack = min_slack;
				best_area = area;
				best_leakage = leakage;
			}
			else
			{
				cout << "slack ok, area ok, power worse than best\n";
				action = 1; //final
				std::cout << "Finishing analysis...." << endl;
				return;
				//goto update_final;
			}
		}
	}
	if (try_other_tile)
	{
		if (M.empty()) {
			cout << "no other options\n"; action = 1; //final
			std::cout << "Finishing analysis...." << endl;
			return;
		}
		tile[M.back().tileidx].is_adaptable_dup = M.back().gate;//change the large sensitivity block
		M.pop_back();
		action = 0; //continue
		std::cout << "Finishing analysis...." << endl;
		return;
		//continue;
	}
	std::cout << "Finishing analysis no jump...." << endl;
	return;
}


void LRACC::sensitivity_commit(CHANGE_LIST& M, std::vector<bool>& tilefixed, float& best_min_slack, int& result_flag, int& tmp_flag, float alpha){
	
	std::cout << "Start commit...." << endl;
	//commit changes
	for (int i = 0; i<tile.size(); i++)
	{
		std::cout << "change made:" << i << " goes from " << tile[i].is_adaptable << " to " << tile[i].is_adaptable_dup << std::endl;
		if (tile[i].is_adaptable == tile[i].is_adaptable_dup) continue;
		tile[i].is_adaptable = tile[i].is_adaptable_dup;
		tilefixed[i] = true;
	}
	result_flag = tmp_flag;
	M.clear();

	//update NPaths
	std::map<TPLGY_VDES, float> NPaths;
	for (TPLGY_VITER viter = boost::vertices(circuit).first; viter != boost::vertices(circuit).second; viter++)
	{
		if (BIT_TEST(circuit[*viter].info->flag, INPUT_DRIVER) || BIT_TEST(circuit[*viter].info->flag, OUTPUT_LOAD) || (circuit[*viter].info->type == FFMS)) { NPaths[*viter] = 0; continue; }
		if (best_min_slack>0)
		{
			NPaths[*viter] = arrival_time_slack[circuit[*viter].info->name].second>0.5*best_min_slack;
		}
		else
		{
			NPaths[*viter] = arrival_time_slack[circuit[*viter].info->name].second<0.5*best_min_slack;
		}
		for (TPLGY_IEITER ieiter = boost::in_edges(*viter, circuit).first; ieiter != boost::in_edges(*viter, circuit).second; ieiter++)
		{
			TPLGY_VDES v = source(*ieiter, circuit);
			if (best_min_slack>0)
			{
				NPaths[*viter] += arrival_time_slack[circuit[*viter].info->name].second>0.5*best_min_slack;
			}
			else
			{
				NPaths[*viter] += arrival_time_slack[circuit[*viter].info->name].second<0.5*best_min_slack;
			}
		}
		for (TPLGY_OEITER ieiter = boost::out_edges(*viter, circuit).first; ieiter != boost::out_edges(*viter, circuit).second; ieiter++)
		{
			TPLGY_VDES v = target(*ieiter, circuit);
			if (best_min_slack>0)
			{
				NPaths[*viter] += arrival_time_slack[circuit[*viter].info->name].second>0.5*best_min_slack;
			}
			else
			{
				NPaths[*viter] += arrival_time_slack[circuit[*viter].info->name].second<0.5*best_min_slack;
			}
		}
	}
	//try add adaptivity
	for (int i = 0; i<tile.size(); i++)
	{
		if (tile[i].components.size() == 0) continue;
		float deltaTNS = 0;//new delay - old delay
		float deltaLeak = tile[i].tuner.power();//new leakage - old leakage
		float deltaArea = tile[i].tuner.area();
		vector<float> prob(MAX_BB_GRADE, 0);
		prob = block_prob[i];
		for (int j = 0; j<tile[i].components.size(); j++)
		{
			TPLGY_VDES v = tile[i].components[j];
			GATEOPT *nowOpt = gateopt_cast(circuit[v].optimal_iter[0]->opt);
			deltaLeak += std::inner_product(prob.begin(), prob.end(), nowOpt->off_power.begin(), 0.0);//expected leakage power
			deltaLeak -= nowOpt->off_power[0];
			float deltaR = std::inner_product(prob.begin(), prob.end(), nowOpt->resist.begin(), 0.0) - nowOpt->resist[0];//expected resistance over probability 
			float C = 0;//calculate the downstream C
			for (TPLGY_OEITER oeiter = boost::out_edges(v, circuit).first; oeiter != boost::out_edges(v, circuit).second; oeiter++)
			{
				TPLGY_VDES tgt = boost::target(*oeiter, circuit);
				GATEOPT *pGate = gateopt_cast(circuit[tgt].optimal_iter[0]->opt);
				C += pGate->cap;
			}
			deltaTNS += deltaR*C*sqrt(NPaths[v]);
		}
		if (best_min_slack>0)//we want to make slack smaller
		{
			if (deltaTNS <= 0 && tile[i].is_adaptable&&!tilefixed[i])//if applying adaptivity makes slack no smaller, and it is adaptable
			{
				//we can potentially make it non-adaptive
				float tmp = -deltaTNS / (deltaLeak>alpha ? deltaLeak : alpha);
				M.push_back(CHANGES(0, tmp, i, true));
			}
			else if (deltaTNS>0 && deltaLeak <= 0 && !tile[i].is_adaptable&&!tilefixed[i])//if applying adaptivity makes delay larger but power no greater, and it is not adaptable
			{
				//we can potentially make it adaptive
				float tmp = deltaTNS / ((-deltaLeak)>alpha ? (-deltaLeak) : alpha);
				M.push_back(CHANGES(1, tmp, i, true));
			}
		}
		else//we want to make slack larger, regardless of area and power
		{
			if (deltaTNS<0 && !tile[i].is_adaptable&&!tilefixed[i])//if applying adaptivity makes slack larger, and it is not adaptable
			{
				//we can potentially make it adaptive
				float tmp = -deltaTNS / tile[i].components.size();///(deltaLeak>alpha?deltaLeak:alpha);
				M.push_back(CHANGES(1, tmp, i, true));
			}
			else if (deltaTNS>0 && tile[i].is_adaptable&&!tilefixed[i])//if applying adaptivity makes slack smaller, and it is adaptable
			{
				//we can potentially make it non-adaptive
				float tmp = deltaTNS / tile[i].components.size();///(-deltaLeak)>alpha?(-deltaLeak):alpha;
				M.push_back(CHANGES(0, tmp, i, true));
			}
		}
		tile[i].is_adaptable_dup = tile[i].is_adaptable;
	}
	if (M.size() == 0) {
		std::cout << "no enhancement!\n"; 				
		action = 1; //final
		std::cout << "Finishing commit...." << endl;
		return;
		// goto update_final;
	}
	M.sort(cmpChanges);
	for (int i = 0; i<1; i++)
	{
		if (best_min_slack>0)//we want to decrease slack
		{
			tile[M.front().tileidx].is_adaptable_dup = M.front().gate;//change the small sensitivity block
			M.pop_front();
		}
		else
		{
			tile[M.back().tileidx].is_adaptable_dup = M.back().gate;//change the large sensitivity block
			M.pop_back();
		}
	}
	std::cout << "Finishing commit...." << endl;
}