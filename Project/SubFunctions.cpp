#include "LRACC.h"
#include <algorithm>
#include <numeric>
#include "Circuit.h"
extern UINT32 Area_Limit;
class CompCHANGES
{
public:
	bool operator()(CHANGES a,CHANGES b) {return a.sensi<b.sensi;}
} cmpChanges;

float LRACC::calculate_leakage()
{
	float leakage = 0;
	for(TPLGY_VITER viter=boost::vertices(circuit).first;viter!=boost::vertices(circuit).second;viter++)
	{
		if(BIT_TEST(circuit[*viter].info->flag,INPUT_DRIVER)||BIT_TEST(circuit[*viter].info->flag,OUTPUT_LOAD)) continue;
		GATEOPT *pGate = gateopt_cast(cell_lib[(circuit[*viter].info)->type][circuit[*viter].candidates[0].back()._opt]);
		int blkidx = (circuit[*viter].info)->tile_idx;
		vector<FLT64>   prob(MAX_BB_GRADE,0);
    	if (!tile[blkidx].is_adaptable)
   		{
				prob[ZERO_GRADE] = 1;
				prob[FBB_GRADE_1] = 0;
				prob[RBB_GRADE_1] = 0;
  	  } 
   	 	else
    	{
        map<int, double>::iterator  iter_prob = probability.find(blkidx);

        if (prob_type[blkidx].compare("0"))
        {prob[ZERO_GRADE] = iter_prob->second;}
        else
        {
            if (prob_type[blkidx].compare("FBB"))
            {prob[FBB_GRADE_1] = iter_prob->second;}
            else if (prob_type[blkidx].compare("RBB"))
            {prob[RBB_GRADE_1] = iter_prob->second;}
            else
            {H102_ASSERT(FALSE);}
            prob[ZERO_GRADE] = 1 - iter_prob->second;
        }
			}
			leakage += std::inner_product(prob.begin(),prob.end(),pGate->off_power.begin(),0.0);//expected leakage power
	}
	return leakage;
}

bool compareMapOn2nd(std::map<std::string,std::pair<double,double> >::value_type a,std::map<std::string,std::pair<double,double> >::value_type b)
{
	return a.second.second<b.second.second;
}

void LRACC::GTR(float alpha,float gamma,float beta)
{
	//sort the cell lib in increasing order

	//initialize to min leakage solution
	for(TPLGY_VITER viter = boost::vertices(circuit).first;viter!=boost::vertices(circuit).second;viter++)
	{
		if(BIT_TEST(circuit[*viter].info->flag,INPUT_DRIVER)||BIT_TEST(circuit[*viter].info->flag,OUTPUT_LOAD)) continue;	
		circuit[*viter].candidates[0].clear();
		circuit[*viter].candidates[1].clear();	
		COMPONENT_SOLUTION s;	
		s._opt = 0;
		circuit[*viter].candidates[0].push_back(s);
		circuit[*viter].candidates[1].push_back(s);
	}
	for(int i=0;i<tile.size();i++)
	{
		tile[i].is_adaptable = false;
		tile[i].is_adaptable_dup = false;
	}
	float best_leakage = 1e100;
	float leakage = 1e100;
	while(1)
	{

		//sta to get slack and compute area
		float area = 0;
		for(TPLGY_VITER viter=boost::vertices(circuit).first;viter!=boost::vertices(circuit).second;viter++)
		{
			if(BIT_TEST(circuit[*viter].info->flag,INPUT_DRIVER)||BIT_TEST(circuit[*viter].info->flag,OUTPUT_LOAD)) continue;
			std::string nm=circuit[*viter].info->name;
			GATEOPT *pGate = gateopt_cast(cell_lib[circuit[*viter].info->type][circuit[*viter].candidates[1].back()._opt]);
			std::vector<float> tmp;
			tmp.push_back(pGate->w);
			tmp.push_back(pGate->resist[0]);
			tmp.push_back(pGate->cap);
			tmp.push_back(pGate->vth);
			gate_sizes[nm] = tmp;
			area += pGate->area;
		}
		for (UINT8 num = 0; num < MAX_TUNING_UNIT; ++num)
		{
 		 	adpt_block_exist[num] = tile[num].is_adaptable_dup;
			if(adpt_block_exist[num])
				area += tile[num].tuner.area(); 
		}	
	
		LRACC_ssta();
		//calculate probability
		get_probability();
		/*for(int i=0;i<tile.size();i++)
		{
			float min_slack = 1e100;
			for(int j=0;j<tile[i].components.size();j++)
			{
				string nm = circuit[tile[i].components[j]].info->name;
				if(arrival_time_slack[nm].second<min_slack) min_slack = arrival_time_slack[nm].second;
			}
			probability[i] = 0.5-min_slack/40;
			if(probability[i]<0) probability[i] = 0;
			else if(probability[i]>1) probability[i] = 1;
		}//*/


		//update best_leakage if feasible
		leakage = calculate_leakage();
		float min_slack = (*std::min_element(arrival_time_slack.begin(),arrival_time_slack.end(),compareMapOn2nd)).second.second;	
		if(min_slack>=0)
		{
			if(leakage<best_leakage)
			{
				best_leakage = leakage;
			}
			else
			{
				for(TPLGY_VITER viter=boost::vertices(circuit).first;viter!=boost::vertices(circuit).second;viter++)
				{
					GATEOPT *pGate = gateopt_cast(cell_lib[circuit[*viter].info->type][circuit[*viter].candidates[0].back()._opt]);
					std::cout<<circuit[*viter].info->name<<":"<<pGate->w<<","<<pGate->vth<<std::endl;
					if(arrival_time_slack[circuit[*viter].info->name].second>=0) continue;
					std::cout<<"slack is negative!!!!!!"<<std::endl;
				}
				for(int n = 0;n<MAX_TUNING_UNIT;n++)
					std::cout<<"block "<<n<<":"<<tile[n].is_adaptable<<std::endl;
				std::cout<<"leakage = "<<leakage<<std::endl;
				break;
			}
		}
		//commit changes
		for(TPLGY_VITER viter=boost::vertices(circuit).first;viter!=boost::vertices(circuit).second;viter++)
		{
			if(BIT_TEST(circuit[*viter].info->flag,INPUT_DRIVER)||BIT_TEST(circuit[*viter].info->flag,OUTPUT_LOAD)) continue;
			circuit[*viter].candidates[0] = circuit[*viter].candidates[1];
		}
		for(int i=0;i<tile.size();i++)
		{
			H102_ASSERT(tile[i].is_adaptable_dup==false);
			tile[i].is_adaptable = tile[i].is_adaptable_dup;
		}			
		std::map<TPLGY_VDES,int> NPaths;
		//update NPaths
		for(TPLGY_VITER viter=boost::vertices(circuit).first;viter!=boost::vertices(circuit).second;viter++)
		{
			if(BIT_TEST(circuit[*viter].info->flag,INPUT_DRIVER)||BIT_TEST(circuit[*viter].info->flag,OUTPUT_LOAD)) continue;	
			NPaths[*viter] = 1;//(arrival_time_slack[circuit[*viter].info->name].second<0);
			for(TPLGY_IEITER ieiter = boost::in_edges(*viter,circuit).first;ieiter!=boost::in_edges(*viter,circuit).second;ieiter++)
			{
				TPLGY_VDES v = source(*ieiter,circuit);
				if(arrival_time_slack[circuit[v].info->name].second<0)
					NPaths[*viter]++;
			}
			for(TPLGY_OEITER ieiter = boost::out_edges(*viter,circuit).first;ieiter!=boost::out_edges(*viter,circuit).second;ieiter++)
			{
				TPLGY_VDES v = target(*ieiter,circuit);
				if(arrival_time_slack[circuit[v].info->name].second<0)
					NPaths[*viter]++;
			}
		}
		//record potential changes and their sensitivity in a set
		std::list<CHANGES> M;
		for(TPLGY_VITER viter=boost::vertices(circuit).first;viter!=boost::vertices(circuit).second;viter++)
		{
			if(BIT_TEST(circuit[*viter].info->flag,INPUT_DRIVER)||BIT_TEST(circuit[*viter].info->flag,OUTPUT_LOAD)) continue;		
			int optidx = circuit[*viter].candidates[0].back()._opt;
			int optmax = cell_lib[(circuit[*viter].info)->type].size();
			int blkidx = (circuit[*viter].info)->tile_idx;
		  vector<FLT64>   prob(MAX_BB_GRADE,0);
    	if (!tile[blkidx].is_adaptable)
   		{
				prob[ZERO_GRADE] = 1;
				prob[FBB_GRADE_1] = 0;
				prob[RBB_GRADE_1] = 0;
  	  } 
   	 	else
    	{
        map<int, double>::iterator  iter_prob = probability.find(blkidx);

        if (prob_type[blkidx].compare("0"))
        {prob[ZERO_GRADE] = iter_prob->second;}
        else
        {
            if (prob_type[blkidx].compare("FBB"))
            {prob[FBB_GRADE_1] = iter_prob->second;}
            else if (prob_type[blkidx].compare("RBB"))
            {prob[RBB_GRADE_1] = iter_prob->second;}
            else
            {H102_ASSERT(FALSE);}
            prob[ZERO_GRADE] = 1 - iter_prob->second;
        }
			}
			if((optidx+1)/10==optidx/10)//if upsizable
			{
				GATEOPT *oldOpt = gateopt_cast(cell_lib[(circuit[*viter].info)->type][optidx]);
				GATEOPT *newOpt = gateopt_cast(cell_lib[(circuit[*viter].info)->type][optidx+1]);	
				float deltaTNS = 0;
				float deltaArea = newOpt->area - oldOpt->area;
				float deltaLeak = std::inner_product(prob.begin(),prob.end(),newOpt->off_power.begin(),0.0)- std::inner_product(prob.begin(),prob.end(),oldOpt->off_power.begin(),0.0);//expected leakage power

				float deltaR = std::inner_product(prob.begin(),prob.end(),newOpt->resist.begin(),0.0) - std::inner_product(prob.begin(),prob.end(),oldOpt->resist.begin(),0.0);//expected resistance over probability 
				float C = 0;//calculate the downstream C
				for(TPLGY_OEITER oeiter=boost::out_edges(*viter,circuit).first;oeiter!=boost::out_edges(*viter,circuit).second;oeiter++)
				{
					TPLGY_VDES tgt = boost::target(*oeiter,circuit);
					GATEOPT *pGate = gateopt_cast(cell_lib[(circuit[tgt].info)->type][circuit[tgt].candidates[0].back()._opt]);
					C += pGate->cap;					
				}
				deltaTNS -= deltaR*C*sqrt(NPaths[*viter]);
	
				float deltaC = newOpt->cap-oldOpt->cap;//C;
				for(TPLGY_IEITER ineiter=boost::in_edges(*viter,circuit).first;ineiter!=boost::in_edges(*viter,circuit).second;ineiter++)
				{
					TPLGY_VDES src = boost::source(*ineiter,circuit);
  				if(BIT_TEST(circuit[src].info->flag,INPUT_DRIVER)||BIT_TEST(circuit[src].info->flag,OUTPUT_LOAD))
					{
						prob.clear();
						prob[ZERO_GRADE] = 1;
						prob[FBB_GRADE_1] = 0;
						prob[RBB_GRADE_1] = 0;
					}
					else
					{
						blkidx = (circuit[src].info)->tile_idx;
						prob.clear();
						if(tile[blkidx].is_adaptable)
						{
							map<int, double>::iterator  iter_prob = probability.find(blkidx);

        			if (prob_type[blkidx].compare("0"))
        			{prob[ZERO_GRADE] = iter_prob->second;}
        			else
        			{
            		if (prob_type[blkidx].compare("FBB"))
            		{prob[FBB_GRADE_1] = iter_prob->second;}
            		else if (prob_type[blkidx].compare("RBB"))
            		{prob[RBB_GRADE_1] = iter_prob->second;}
            		else
            		{H102_ASSERT(FALSE);}
            		//local_prob[TUNING_GRADE_H] = iter_prob->second;//supposing only the prob for changing L is given
            		prob[ZERO_GRADE] = 1 - iter_prob->second;
        			}
						}
						else
						{
							prob[ZERO_GRADE] = 1;
							prob[FBB_GRADE_1] = 0;
							prob[RBB_GRADE_1] = 0;
						}
					}
					GATEOPT *pGate = gateopt_cast(cell_lib[(circuit[src].info)->type][circuit[src].candidates[0].back()._opt]); 
					float R = std::inner_product(prob.begin(),prob.end(),pGate->resist.begin(),0.0);//expected resistance over probability
					deltaTNS -= R*deltaC*sqrt(NPaths[src]);
				}
				float tmp = deltaTNS/(deltaLeak>alpha?deltaLeak:alpha);
				//tmp = tmp/exp(deltaArea*beta/((Area_Limit-area)>1?(Area_Limit-area):1));
				if(tmp>0)
					M.push_back(CHANGES(*viter,tmp,optidx+1));
			}
			if(optidx+10<30)//if vth downscalable
			{
				GATEOPT *oldOpt = gateopt_cast(cell_lib[(circuit[*viter].info)->type][optidx]);
				GATEOPT *newOpt = gateopt_cast(cell_lib[(circuit[*viter].info)->type][optidx+10]);	
				float deltaTNS = 0;
				float deltaArea = newOpt->area - oldOpt->area;
				float deltaLeak = std::inner_product(prob.begin(),prob.end(),newOpt->off_power.begin(),0.0);
				deltaLeak -= std::inner_product(prob.begin(),prob.end(),oldOpt->off_power.begin(),0.0);//expected leakage power

				float deltaR = std::inner_product(prob.begin(),prob.end(),newOpt->resist.begin(),0.0) - std::inner_product(prob.begin(),prob.end(),oldOpt->resist.begin(),0.0);//expected resistance over probability 
				float C = 0;//calculate the downstream C
				for(TPLGY_OEITER oeiter=boost::out_edges(*viter,circuit).first;oeiter!=boost::out_edges(*viter,circuit).second;oeiter++)
				{
					TPLGY_VDES tgt = boost::target(*oeiter,circuit);
					GATEOPT *pGate = gateopt_cast(cell_lib[(circuit[tgt].info)->type][circuit[tgt].candidates[0].back()._opt]);
					C += pGate->cap;					
				}
				deltaTNS -= deltaR*C*sqrt(NPaths[*viter]);
	
				float deltaC = newOpt->cap-oldOpt->cap;//C;
				for(TPLGY_IEITER ineiter=boost::in_edges(*viter,circuit).first;ineiter!=boost::in_edges(*viter,circuit).second;ineiter++)
				{
					TPLGY_VDES src = boost::source(*ineiter,circuit);
					if(BIT_TEST(circuit[src].info->flag,INPUT_DRIVER)||BIT_TEST(circuit[src].info->flag,OUTPUT_LOAD))
					{
						prob.clear();
						prob[ZERO_GRADE] = 1;
						prob[FBB_GRADE_1] = 0;
						prob[RBB_GRADE_1] = 0;
					}
					else
					{
						blkidx = (circuit[src].info)->tile_idx;
						prob.clear();
						if(tile[blkidx].is_adaptable)
						{
							map<int, double>::iterator  iter_prob = probability.find(blkidx);

        			if (prob_type[blkidx].compare("0"))
        			{prob[ZERO_GRADE] = iter_prob->second;}
        			else
        			{
            		if (prob_type[blkidx].compare("FBB"))
            		{prob[FBB_GRADE_1] = iter_prob->second;}
            		else if (prob_type[blkidx].compare("RBB"))
            		{prob[RBB_GRADE_1] = iter_prob->second;}
            		else
            		{H102_ASSERT(FALSE);}
            		prob[ZERO_GRADE] = 1 - iter_prob->second;
        			}
						}
						else
						{
							prob[ZERO_GRADE] = 1;
							prob[FBB_GRADE_1] = 0;
							prob[RBB_GRADE_1] = 0;
						}
					}
					GATEOPT *pGate = gateopt_cast(cell_lib[(circuit[src].info)->type][circuit[src].candidates[0].back()._opt]); 
					float R = std::inner_product(prob.begin(),prob.end(),pGate->resist.begin(),0.0);//expected resistance over probability
					deltaTNS -= R*deltaC*sqrt(NPaths[src]);
				}
		
				float tmp = deltaTNS/(deltaLeak>alpha?deltaLeak:alpha);///(deltaArea>alpha?deltaArea:alpha);
				if(tmp>0)
					M.push_back(CHANGES(*viter,tmp,optidx+10));
			}
		}
		
		//try add adaptivity
		for(int i=0;i<tile.size();i++)
		{
			if(tile[i].is_adaptable) continue;
			float deltaTNS = 0;
			float deltaLeak = tile[i].tuner.power();
			float deltaArea = tile[i].tuner.area();
			vector<FLT64>   prob(MAX_BB_GRADE,0);
			prob[ZERO_GRADE] = 1;
			prob[FBB_GRADE_1] = 0;
			prob[RBB_GRADE_1] = 0;
			for(int j=0;j<tile[i].components.size();j++)
			{
				TPLGY_VDES v = tile[i].components[j];
				int optidx = circuit[v].candidates[0].back()._opt;
				GATEOPT *nowOpt = gateopt_cast(cell_lib[(circuit[v].info)->type][optidx]);
				deltaLeak += std::inner_product(prob.begin(),prob.end(),nowOpt->off_power.begin(),0.0);//expected leakage power
				deltaLeak -= nowOpt->off_power[0];
				float deltaR = std::inner_product(prob.begin(),prob.end(),nowOpt->resist.begin(),0.0) - nowOpt->resist[0];//expected resistance over probability 
				float C = 0;//calculate the downstream C
				for(TPLGY_OEITER oeiter=boost::out_edges(v,circuit).first;oeiter!=boost::out_edges(v,circuit).second;oeiter++)
				{
					TPLGY_VDES tgt = boost::target(*oeiter,circuit);
					GATEOPT *pGate = gateopt_cast(cell_lib[(circuit[tgt].info)->type][circuit[tgt].candidates[0].back()._opt]);
					C += pGate->cap;					
				}
				deltaTNS -= deltaR*C*sqrt(NPaths[v]);
			}
			float tmp = deltaTNS/(deltaLeak>alpha?deltaLeak:alpha);
			//tmp = tmp/exp(deltaArea*beta/((Area_Limit-area)>1?(Area_Limit-area):1));
			if(tmp>0)
				M.push_back(CHANGES(*vertices(circuit).first,tmp,i,true));
		}
	
		M.sort(cmpChanges);
		for(int i=0;i<gamma*M.size();i++)
		{
			if(!M.back().block_change)
				circuit[M.back().gate].candidates[1].back()._opt = M.back().newopt;
			else
				tile[M.back().newopt].is_adaptable_dup = true;		
			M.pop_back();
		}
		M.clear();
	}
}



bool LRACC::SGGS(float alpha)
{
	//sort the cell lib in increasing order

	float leakage;
	std::map<TPLGY_VDES,int> NPaths;
	//update NPaths
	for(TPLGY_VITER viter=boost::vertices(circuit).first;viter!=boost::vertices(circuit).second;viter++)
	{
		if(BIT_TEST(circuit[*viter].info->flag,INPUT_DRIVER)||BIT_TEST(circuit[*viter].info->flag,OUTPUT_LOAD)) continue;	
		NPaths[*viter] = 1;//(arrival_time_slack[circuit[*viter].info->name].second<0);
		for(TPLGY_IEITER ieiter = boost::in_edges(*viter,circuit).first;ieiter!=boost::in_edges(*viter,circuit).second;ieiter++)
		{
			TPLGY_VDES v = source(*ieiter,circuit);
			if(arrival_time_slack[circuit[v].info->name].second<0)
				NPaths[*viter]++;
		}
		for(TPLGY_OEITER ieiter = boost::out_edges(*viter,circuit).first;ieiter!=boost::out_edges(*viter,circuit).second;ieiter++)
		{
			TPLGY_VDES v = target(*ieiter,circuit);
			if(arrival_time_slack[circuit[v].info->name].second<0)
			NPaths[*viter]++;
		}
	}
	std::list<CHANGES> M;
	for(TPLGY_VITER viter=boost::vertices(circuit).first;viter!=boost::vertices(circuit).second;viter++)
	{
		if(BIT_TEST(circuit[*viter].info->flag,INPUT_DRIVER)||BIT_TEST(circuit[*viter].info->flag,OUTPUT_LOAD)) continue;		
		int optidx = circuit[*viter].candidates[0].back()._opt;
		int optmax = cell_lib[(circuit[*viter].info)->type].size();
		int blkidx = (circuit[*viter].info)->tile_idx;
		vector<FLT64>   prob(MAX_BB_GRADE,0);
		if(tile[blkidx].is_adaptable)
		{
			map<int, double>::iterator  iter_prob = probability.find(blkidx);
 			if (prob_type[blkidx].compare("0"))
 			{prob[ZERO_GRADE] = iter_prob->second;}
 			else
 			{
     		if (prob_type[blkidx].compare("FBB"))
     		{prob[FBB_GRADE_1] = iter_prob->second;}
     		else if (prob_type[blkidx].compare("RBB"))
     		{prob[RBB_GRADE_1] = iter_prob->second;}
     		else
     		{H102_ASSERT(FALSE);}
     		prob[ZERO_GRADE] = 1 - iter_prob->second;
   		}
		}
		else
		{
			prob[ZERO_GRADE] = 1;
			prob[FBB_GRADE_1] = 0;
			prob[RBB_GRADE_1] = 0;
		}
		if((optidx>0)&&((optidx-1)/10==optidx/10))//if downsizable
		{
			GATEOPT *oldOpt = gateopt_cast(cell_lib[(circuit[*viter].info)->type][optidx]);
			GATEOPT *newOpt = gateopt_cast(cell_lib[(circuit[*viter].info)->type][optidx-1]);	
			float deltaTNS = 0;
			float deltaArea = newOpt->area - oldOpt->area;
			float deltaLeak = std::inner_product(prob.begin(),prob.end(),newOpt->off_power.begin(),0.0);
			deltaLeak -= std::inner_product(prob.begin(),prob.end(),oldOpt->off_power.begin(),0.0);//expected leakage power

			float deltaR = std::inner_product(prob.begin(),prob.end(),newOpt->resist.begin(),0.0) - std::inner_product(prob.begin(),prob.end(),oldOpt->resist.begin(),0.0);//expected resistance over probability 
			float C = 0;//calculate the downstream C
			for(TPLGY_OEITER oeiter=boost::out_edges(*viter,circuit).first;oeiter!=boost::out_edges(*viter,circuit).second;oeiter++)
			{
				TPLGY_VDES tgt = boost::target(*oeiter,circuit);
				GATEOPT *pGate = gateopt_cast(cell_lib[(circuit[tgt].info)->type][circuit[tgt].candidates[0].back()._opt]);
				C += pGate->cap;					
			}
			deltaTNS -= deltaR*C*sqrt(NPaths[*viter]);
	
			float deltaC = newOpt->cap-oldOpt->cap;//C;
			for(TPLGY_IEITER ineiter=boost::in_edges(*viter,circuit).first;ineiter!=boost::in_edges(*viter,circuit).second;ineiter++)
			{
				TPLGY_VDES src = boost::source(*ineiter,circuit);
				if(BIT_TEST(circuit[src].info->flag,INPUT_DRIVER)||BIT_TEST(circuit[src].info->flag,OUTPUT_LOAD))
					{
						prob.clear();
						prob[ZERO_GRADE] = 1;
						prob[FBB_GRADE_1] = 0;
						prob[RBB_GRADE_1] = 0;
					}
					else
					{
						blkidx = (circuit[src].info)->tile_idx;
						prob.clear();
						if(tile[blkidx].is_adaptable)
						{
							map<int, double>::iterator  iter_prob = probability.find(blkidx);

        			if (prob_type[blkidx].compare("0"))
        			{prob[ZERO_GRADE] = iter_prob->second;}
        			else
        			{
            		if (prob_type[blkidx].compare("FBB"))
            		{prob[FBB_GRADE_1] = iter_prob->second;}
            		else if (prob_type[blkidx].compare("RBB"))
            		{prob[RBB_GRADE_1] = iter_prob->second;}
            		else
            		{H102_ASSERT(FALSE);}
            		prob[ZERO_GRADE] = 1 - iter_prob->second;
        			}
						}
						else
						{
							prob[ZERO_GRADE] = 1;
							prob[FBB_GRADE_1] = 0;
							prob[RBB_GRADE_1] = 0;
						}
					}
				GATEOPT *pGate = gateopt_cast(cell_lib[(circuit[src].info)->type][circuit[src].candidates[0].back()._opt]); 
				float R = std::inner_product(prob.begin(),prob.end(),pGate->resist.begin(),0.0);//expected resistance over probability
				deltaTNS -= R*deltaC*sqrt(NPaths[src]);
			}
			float tmp = deltaLeak/(deltaTNS>alpha?deltaTNS:alpha);
				M.push_back(CHANGES(*viter,tmp,optidx-1));
		}
		if(optidx-10>=0)//if vth upscalable
		{
			GATEOPT *oldOpt = gateopt_cast(cell_lib[(circuit[*viter].info)->type][optidx]);
			GATEOPT *newOpt = gateopt_cast(cell_lib[(circuit[*viter].info)->type][optidx-10]);	
			float deltaTNS = 0;
			float deltaArea = newOpt->area - oldOpt->area;
			float deltaLeak = std::inner_product(prob.begin(),prob.end(),newOpt->off_power.begin(),0.0);
			deltaLeak -= std::inner_product(prob.begin(),prob.end(),oldOpt->off_power.begin(),0.0);//expected leakage power

			float deltaR = std::inner_product(prob.begin(),prob.end(),newOpt->resist.begin(),0.0) - std::inner_product(prob.begin(),prob.end(),oldOpt->resist.begin(),0.0);//expected resistance over probability 
			float C = 0;//calculate the downstream C
			for(TPLGY_OEITER oeiter=boost::out_edges(*viter,circuit).first;oeiter!=boost::out_edges(*viter,circuit).second;oeiter++)
			{
				TPLGY_VDES tgt = boost::target(*oeiter,circuit);
				GATEOPT *pGate = gateopt_cast(cell_lib[(circuit[tgt].info)->type][circuit[tgt].candidates[0].back()._opt]);
				C += pGate->cap;					
			}
			deltaTNS -= deltaR*C*sqrt(NPaths[*viter]);

			float deltaC = newOpt->cap-oldOpt->cap;//C;
			for(TPLGY_IEITER ineiter=boost::in_edges(*viter,circuit).first;ineiter!=boost::in_edges(*viter,circuit).second;ineiter++)
			{
				TPLGY_VDES src = boost::source(*ineiter,circuit);
					if(BIT_TEST(circuit[src].info->flag,INPUT_DRIVER)||BIT_TEST(circuit[src].info->flag,OUTPUT_LOAD))
					{
						prob.clear();
						prob[ZERO_GRADE] = 1;
						prob[FBB_GRADE_1] = 0;
						prob[RBB_GRADE_1] = 0;
					}
					else
					{
						blkidx = (circuit[src].info)->tile_idx;
						prob.clear();
						if(tile[blkidx].is_adaptable)
						{
							map<int, double>::iterator  iter_prob = probability.find(blkidx);

        			if (prob_type[blkidx].compare("0"))
        			{prob[ZERO_GRADE] = iter_prob->second;}
        			else
        			{
            		if (prob_type[blkidx].compare("FBB"))
            		{prob[FBB_GRADE_1] = iter_prob->second;}
            		else if (prob_type[blkidx].compare("RBB"))
            		{prob[RBB_GRADE_1] = iter_prob->second;}
            		else
            		{H102_ASSERT(FALSE);}
            		prob[ZERO_GRADE] = 1 - iter_prob->second;
        			}
						}
						else
						{
							prob[ZERO_GRADE] = 1;
							prob[FBB_GRADE_1] = 0;
							prob[RBB_GRADE_1] = 0;
						}
					}
				GATEOPT *pGate = gateopt_cast(cell_lib[(circuit[src].info)->type][circuit[src].candidates[0].back()._opt]); 
				float R = std::inner_product(prob.begin(),prob.end(),pGate->resist.begin(),0.0);//expected resistance over probability
				deltaTNS -= R*deltaC*sqrt(NPaths[src]);
			}	
			float tmp = deltaLeak/(deltaTNS>alpha?deltaTNS:alpha);	
				M.push_back(CHANGES(*viter,tmp,optidx-10));
		}
	}
	M.sort(cmpChanges);
	bool flag = false;
	for(int i=0;i<M.size();i++)
	{
		//sta to get slack
		///*
		for(TPLGY_VITER viter=boost::vertices(circuit).first;viter!=boost::vertices(circuit).second;viter++)
		{
			if(BIT_TEST(circuit[*viter].info->flag,INPUT_DRIVER)||BIT_TEST(circuit[*viter].info->flag,OUTPUT_LOAD)) continue;
			std::string nm=circuit[*viter].info->name;
			GATEOPT *pGate = gateopt_cast(cell_lib[circuit[*viter].info->type][circuit[*viter].candidates[0].back()._opt]);
			std::vector<float> tmp;
			tmp.push_back(pGate->w);
			tmp.push_back(pGate->resist[0]);
			tmp.push_back(pGate->cap);
			tmp.push_back(pGate->vth);
			gate_sizes[nm] = tmp; 
			//std::cout<<"name = "<<nm<<",size = "<<tmp[0]<<",vth = "<<tmp[3]<<std::endl;
		}
		std::string nm = circuit[M.back().gate].info->name;
		GATEOPT *pGate = gateopt_cast(cell_lib[circuit[M.back().gate].info->type][M.back().newopt]);
		std::vector<float> tmp;
		tmp.push_back(pGate->w);
		tmp.push_back(pGate->resist[0]);
		tmp.push_back(pGate->cap);
		tmp.push_back(pGate->vth);
		gate_sizes[nm] = tmp; 
		
 		for (UINT8 num = 0; num < MAX_TUNING_UNIT; ++num)
 		{
  	 	adpt_block_exist[num] = tile[num].is_adaptable;
			//std::cout<<"tile "<<num<<" = "<<adpt_block_exist[num]<<std::endl;
		}	
		LRACC_ssta();
		float min_slack = (*std::min_element(arrival_time_slack.begin(),arrival_time_slack.end(),compareMapOn2nd)).second.second;	
		if(min_slack>=0)	
		//*/
		//if(ISTAfeasible(M.back()))
		{
			std::cout<<"sggs!"<<std::endl;
			flag = true;
			circuit[M.back().gate].candidates[0].back()._opt = M.back().newopt;
			H102_ASSERT(check_feasibility());
		}
		M.pop_back();
	}
	return flag;
}

/*bool LRACC::ISTAfeasible(CHANGES c)
{
	std::map<std::string,std::pair<double,double> > DUP = arrival_time_slack;
	float Cap=0;
	GATEOPT *oldOpt = gateopt_cast(cell_lib[circuit[c.gate].info->type][circuit[c.gate].candidates[0].back()._opt]); 
	GATEOPT *NewOpt = gateopt_cast(cell_lib[circuit[c.gate].info->type][c.newopt]); 
	for(TPLGY_OEITER oeiter=boost::out_edges(c.gate,circuit).first;oeiter!=boost::out_edges(c.gate,circuit).second;oeiter++)
	{
		Cap +=	gateopt_cast(cell_lib[circuit[boost::target(*oeiter,circuit)].info->type][circuit[boost::target(*oeiter,circuit)].candidates[0].back()._opt])->cap;
	}
	int blkidx = (circuit[c.gate].info)->tile_idx;
	std::vector<float> prob;
	if(tile[blkidx].is_adaptable)
	{
		prob.push_back(1-probability[blkidx]);
		prob.push_back(probability[blkidx]);
	}
	else
	{
		prob.push_back(1);
		prob.push_back(0);
	} 
	float deltaR = prob[0]*(NewOpt->resist[0]-oldOpt->resist[0])+prob[1]*(NewOpt->resist[1]-oldOpt->resist[1]);
	float deltaDelay = deltaR*Cap;
	if(arrival_time_slack[circuit[c.gate].info->name].second-deltaDelay<5)//delay changes larger than slack
		return false;
	else
		DUP[circuit[c.gate].info->name].second = deltaDelay+arrival_time_slack[circuit[c.gate].info->name].second;

	float deltaC = NewOpt->cap-oldOpt->cap;
	for(TPLGY_IEITER ieiter = boost::in_edges(c.gate,circuit).first;ieiter!=boost::in_edges(c.gate,circuit).second;ieiter++)
	{
		TPLGY_VDES v = boost::source(*ieiter,circuit);
		if(BIT_TEST(circuit[v].info->flag,INPUT_DRIVER)||BIT_TEST(circuit[v].info->flag,OUTPUT_LOAD)) continue;	
		GATEOPT *Opt = gateopt_cast(cell_lib[circuit[v].info->type][circuit[v].candidates[0].back()._opt]);
		int blkidx = (circuit[v].info)->tile_idx;
		std::vector<float> prob;
		if(tile[blkidx].is_adaptable)
		{
			prob.push_back(1-probability[blkidx]);
			prob.push_back(probability[blkidx]);
		}
		else
		{
			prob.push_back(1);
			prob.push_back(0);
		} 
		float deltaDelay = std::inner_product(prob.begin(),prob.end(),Opt->resist.begin(),0.0)*deltaC;
		if(arrival_time_slack[circuit[v].info->name].second-deltaDelay<5)//delay changes larger than slack
			return false;
		else
			DUP[circuit[c.gate].info->name].second = deltaDelay+arrival_time_slack[circuit[c.gate].info->name].second;
	}
	//update slack
	arrival_time_slack = DUP;
	return true;
}*/

void LRACC::SpeedUpBottleneck()
{
	//sta to get slack
	for(TPLGY_VITER viter=boost::vertices(circuit).first;viter!=boost::vertices(circuit).second;viter++)
	{
		if(BIT_TEST(circuit[*viter].info->flag,INPUT_DRIVER)||BIT_TEST(circuit[*viter].info->flag,OUTPUT_LOAD)) continue;
		std::string nm=circuit[*viter].info->name;
		GATEOPT *pGate = gateopt_cast(cell_lib[circuit[*viter].info->type][circuit[*viter].candidates[0].back()._opt]);
		std::vector<float> tmp;
		tmp.push_back(pGate->w);
		tmp.push_back(pGate->resist[0]);
		tmp.push_back(pGate->cap);
		tmp.push_back(pGate->vth);
		gate_sizes[nm] = tmp;     	
	}
	for (UINT8 num = 0; num < MAX_TUNING_UNIT; ++num)
  {
  	adpt_block_exist[num] = tile[num].is_adaptable;
	}	
	
	LRACC_ssta();

	//float min_slack = (*std::min_element(arrival_time_slack.begin(),arrival_time_slack.end(),compareMapOn2nd)).second.second;
	//if(min_slack>=0) return;
	for(TPLGY_VITER viter=boost::vertices(circuit).first;viter!=boost::vertices(circuit).second;viter++)
	{
		if(BIT_TEST(circuit[*viter].info->flag,INPUT_DRIVER)||BIT_TEST(circuit[*viter].info->flag,OUTPUT_LOAD)) continue;
		if(arrival_time_slack[circuit[*viter].info->name].second>=0) continue;
		std::cout<<"speedup!"<<std::endl;
		int optidx = circuit[*viter].candidates[0].back()._opt;
		if((optidx+1)/10==optidx/10) circuit[*viter].candidates[0].back()._opt = optidx+1;
		if(optidx+10<30) circuit[*viter].candidates[0].back()._opt = optidx+10;
	}
}

bool LRACC::check_feasibility()
{
		//sta to get slack
		for(TPLGY_VITER viter=boost::vertices(circuit).first;viter!=boost::vertices(circuit).second;viter++)
		{
			if(BIT_TEST(circuit[*viter].info->flag,INPUT_DRIVER)||BIT_TEST(circuit[*viter].info->flag,OUTPUT_LOAD)) continue;
			std::string nm=circuit[*viter].info->name;
			GATEOPT *pGate = gateopt_cast(cell_lib[circuit[*viter].info->type][circuit[*viter].candidates[0].back()._opt]);
			std::vector<float> tmp;
			tmp.push_back(pGate->w);
			tmp.push_back(pGate->resist[0]);
			tmp.push_back(pGate->cap);
			tmp.push_back(pGate->vth);
			gate_sizes[nm] = tmp;
#ifdef PRINT_SSTA
			std::cout<<"name = "<<nm<<":r = "<<tmp[1]<<",c = "<<tmp[2]<<std::endl;
#endif
		}
		for (UINT8 num = 0; num < MAX_TUNING_UNIT; ++num)
		{
 		 	adpt_block_exist[num] = tile[num].is_adaptable;
#ifdef PRINT_SSTA
			std::cout<<"tile "<<num<<" = "<<adpt_block_exist[num]<<std::endl;	
#endif
		}	
	
		LRACC_ssta();
#ifdef PRINT_SSTA
		for(TPLGY_VITER viter=boost::vertices(circuit).first;viter!=boost::vertices(circuit).second;viter++)
		{
			std::string nm=circuit[*viter].info->name;	
			std::cout<<"name = "<<nm<<":arrival time = "<<arrival_time_slack[nm].first<<",slack = "<<arrival_time_slack[nm].second<<std::endl;
		}
#endif
		//update best_leakage if feasible
		float min_slack = (*std::min_element(arrival_time_slack.begin(),arrival_time_slack.end(),compareMapOn2nd)).second.second;	
		return min_slack>=0;	
}


