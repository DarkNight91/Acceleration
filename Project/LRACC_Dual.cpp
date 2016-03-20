#include "LRACC.h"
#include "werc_102_vlsi.h"
#include <string>
#include <cstdio>
extern FLT32 Time_Limit;
extern UINT32D Area_Limit;
void LRACC::LRACC_dual()
{
///*
 	//float stepsize = 0.00025*(1-iter_lagrangian+MAX_LAGRANGIAN_NUM);
  float stepsize = 0.001;
  cout<<"#"<<iter_lagrangian<< "\t LRACC_DUAL parameters"<< endl;
	for(TPLGY_EITER eiter = edges(circuit).first;eiter!=edges(circuit).second;eiter++)
	{
		TPLGY_VDES src = source(*eiter,circuit);
		TPLGY_VDES tgt = target(*eiter,circuit);
		float ai,aj,di;
		if ((circuit[src].info->type == IN1) || (circuit[src].info->type == FFMS&&in_degree(src, circuit) == 0)){
			aj = 0;
		}
		else{
			aj = mu_sigma[circuit[src].info->name].first;
		}
		if ((circuit[tgt].info->type == OUT1) || (circuit[tgt].info->type == FFMS&&out_degree(tgt, circuit) == 0))
		{
			ai = Time_Limit;
			di = 0;
		}
		else
		{
			ai = mu_sigma[circuit[tgt].info->name].first;
			di = gate_sizes[circuit[tgt].info->name][1]*gate_sizes[circuit[tgt].info->name][2];
		}
		circuit[*eiter].lm+=stepsize*(aj+di-ai);
		//cout<<"new lm = "<<circuit[*eiter].lm<<std::endl;
		if(circuit[*eiter].lm<0)
			circuit[*eiter].lm = 0;
	}
	//project to KKT condition space
	for (VDES_CONTAINER::iterator ii = order.begin(); ii != order.end(); ++ii)
  {
		if(circuit[*ii].info->type==OUT1) continue;
		if(circuit[*ii].info->type==FFMS&&out_degree(*ii,circuit)==0) continue;
		float out_lm = 0;
		for(TPLGY_OEITER oeiter = out_edges(*ii,circuit).first;oeiter!=out_edges(*ii,circuit).second;oeiter++)
		{
			out_lm+=circuit[*oeiter].lm;
		}
		circuit[*ii].lm = out_lm;
		float in_lm = 0;
		for(TPLGY_IEITER ieiter = in_edges(*ii,circuit).first;ieiter!=in_edges(*ii,circuit).second;ieiter++)
		{
			in_lm+=circuit[*ieiter].lm;
		}
		if(in_lm==0) continue;
		for(TPLGY_IEITER ieiter = in_edges(*ii,circuit).first;ieiter!=in_edges(*ii,circuit).second;ieiter++)
		{
			circuit[*ieiter].lm*=out_lm/in_lm;
		}
	}
	float tmp = (total_a-Area_Limit)/Area_Limit; 
	//if(tmp>0.01) tmp = (1*tmp);
	//else if(tmp<0.01&&tmp>0) tmp = 0.01*1;
	//lm_a += stepsize*tmp;
  //if (lm_a < 0)
  {
  	lm_a = 0;
  }
	cout<<"Area_lm ="<<lm_a<<endl<<endl;
	//*/
	//float stepsize = 0.002*(1-iter_lagrangian+MAX_LAGRANGIAN_NUM);
 	/*float stepsize = 0.02;
  cout<<"#"<<iter_lagrangian<< "\t LRACC_DUAL parameters"<< endl;
	//update lm based on this hyperplane
	for (VDES_CONTAINER::iterator ii = order.begin(); ii != order.end(); ++ii)
  {
		TPLGY_VDES fgate = *ii;
    if (BIT_TEST(circuit[fgate].info->flag, OUTPUT_LOAD))
    {
    	continue;//skip the OUTPUT gate
    }
   	if (BIT_TEST(circuit[fgate].info->flag, INPUT_DRIVER))
    {
    	circuit[fgate].lm = 0;///???????????? just for simplicity
      continue;
    } 
 		string nm = (circuit[fgate].info)->name;
		//cout<<nm<<": "<<circuit[fgate].lm;
		float tmp = -arrival_time_slack[nm].second/Time_Limit;
		if(tmp>=0.01) tmp = (15*tmp);
		else if(tmp<0.01&&tmp>0) tmp = 15*0.01;
		tmp*=gate_sizes[circuit[fgate].info->name][1]*gate_sizes[circuit[fgate].info->name][2];
		circuit[fgate].lm += stepsize*tmp;
    if (circuit[fgate].lm <= 0)
    {
    	circuit[fgate].lm = 0;
    }
		//cout<<" -> lm="<<circuit[fgate].lm<<endl;//the Lagrangian value after calculate
  }
	float tmp = (total_a-Area_Limit)/Area_Limit; 
	//if(tmp>0.01) tmp = (1*tmp);
	//else if(tmp<0.01&&tmp>0) tmp = 0.01*1;
	//lm_a += stepsize*tmp;
  //if (lm_a < 0)
  {
  	lm_a = 0;
  }
	cout<<"Area_lm ="<<lm_a<<endl<<endl;//*/


	/*
	//update history based on last iteration
	dual_hist_obj.push_back(total_p);
	dual_hist_subg.push_back(make_pair(total_a-Area_Limit,arrival_time_slack));	
	//find the supporting hyperplane(this hyperplane will correspond to a history)
	vector<float> min(dual_hist_obj.size(),0);
	for(int i=0;i<dual_hist_obj.size();i++)
	{
		for (VDES_CONTAINER::iterator ii = order.begin(); ii != order.end(); ++ii)
    {
        fgate = *ii;
        if (BIT_TEST(circuit[fgate].info->flag, OUTPUT_LOAD||INPUT_DRIVER))
        {
            continue;//skip the OUTPUT gate
        }
    	  string nm = (circuit[fgate].info)->name;
        min[i] -= ((dual_hist_subg[i].second)[nm].second)*circuit[fgate].lm;
		}
		min[i] += (dual_hist_obj[i]);
		min[i] += (dual_hist_subg[i].first)*lm_a;
	}
	int idx = distance(min.begin(),min_element(min.begin(),min.end()));
	//update lm based on this hyperplane
	for (VDES_CONTAINER::iterator ii = order.begin(); ii != order.end(); ++ii)
  {
		fgate = *ii;
    if (BIT_TEST(circuit[fgate].info->flag, OUTPUT_LOAD))
    {
    	continue;//skip the OUTPUT gate
    }
   	if (BIT_TEST(circuit[fgate].info->flag, INPUT_DRIVER))
    {
    	circuit[fgate].lm = 0;///???????????? just for simplicity
      continue;
    } 
	
 		string nm = (circuit[fgate].info)->name;
		float tmp = -(dual_hist_subg.back().second)[nm].second*(dual_hist_subg.back().second)[nm].first/Time_Limit;
		if(tmp>=0.01) tmp = (30*tmp);
		else if(tmp<0.01&&tmp>0) tmp = 30*0.01;
		circuit[fgate].lm += stepsize*tmp;
    if (circuit[fgate].lm <= 0)
    {
    	circuit[fgate].lm = 0;
    }
		else cout<<nm<<" -> lm="<<circuit[fgate].lm<<endl;//the Lagrangian value after calculate
  }
	float tmp = (dual_hist_subg.back().first)/Area_Limit; 
	if(tmp>0.01) tmp = (20*tmp);
	else if(tmp<0.01&&tmp>0) tmp = 20*0.05;
	lm_a += stepsize*tmp;
  if (lm_a < 0)
  {
  	lm_a = 0;
  }
	cout<<"Area_lm ="<<lm_a<<endl<<endl;//*/
}

