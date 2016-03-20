#include <fstream>
#include <string>
#include <vector>
using namespace std;

template <class DataType>
bool Store(string filename, DataType data)
{
	ofstream f(filename.c_str(),ofstream::app|ofstream::out);
	if(!f.is_open()) {f.close();return false;}
	f<<'#'<<endl;//start data
	f<<1<<' '<<data;//1 element
	f<<endl<<'$'<<endl;//end data
	f.close();
	return true;
}
template bool Store(string, float);
template bool Store(string, unsigned int);
template bool Store(string, string);
template bool Store(string, unsigned long);
template bool Store(string, unsigned char);


template <class DataType>
bool Store(string filename, vector<DataType> data)
{
	ofstream f(filename.c_str(),ofstream::app|ofstream::out);
	if(!f.is_open()) {f.close();return false;}
	f<<'#'<<endl;//start data
	f<<data.size();//m elements
	for(typename vector<DataType>::iterator iter=data.begin();iter!=data.end();iter++)
		f<<' '<<(*iter);
	f<<endl<<'$'<<endl;//end data
	f.close();
	return true;
}
template bool Store(string, vector<float>);
template bool Store(string, vector<unsigned int>);
template bool Store(string, vector<unsigned long>);

template <class DataType>
bool Store(string filename, vector< vector<DataType> > data)
{
	ofstream f(filename.c_str(),ofstream::app|ofstream::out);
	if(!f.is_open()) {f.close();return false;}
	f<<'#'<<endl;//start data
	for(typename vector< vector<DataType> >::iterator iter1=data.begin();iter1!=data.end();iter1++)
	{
		f<<iter1->size();
		for(typename vector<DataType>::iterator iter2=iter1->begin();iter2!=iter1->end();iter2++)
			f<<' '<<*iter2;
		f<<endl;
	}
	f<<'$'<<endl;//end data
	return true;
}
template bool Store(string, vector< vector<float> >);
template bool Store(string, vector< vector<int> >);
template bool Store(string, vector< vector<unsigned int> >);
template bool Store(string, vector< vector<bool> >);

template <class DataType>
bool Read(string filename,unsigned int N,DataType &data)
{
	ifstream f(filename.c_str(),ifstream::in);
	if(!f.is_open()) {f.close();return false;}
	while(N)
	{
		char tmp;
		tmp = f.get();//extract the next character until...
		if(!f.good()) {f.close();return false;}
		if(tmp=='#') N--;//...N-th episode begin indicator
	}
	f.ignore();//ignore the next '\n'
	unsigned int nbr;
	DataType tmp;
	char end;
	f>>nbr>>tmp>>end;
	if(!f.good()||nbr!=1||end!='$') {f.close();return false;}
	data = tmp;
	f.close();
	return true;
}
template bool Read(string, unsigned int, float &);
template bool Read(string, unsigned int, unsigned int &);
template bool Read(string, unsigned int, unsigned long &);
template bool Read(string, unsigned int, unsigned char &);


template <class DataType>
bool Read(string filename,unsigned int N,vector<DataType> &data)
{
	ifstream f(filename.c_str(),ifstream::in);
	if(!f.is_open()) {f.close();return false;}
	while(N)
	{
		char tmp;
		tmp = f.get();//extract the next character until...
		if(!f.good()) {f.close();return false;}
		if(tmp=='#') N--;//...N-th episode begin indicator
	}
	f.ignore();//ignore the next '\n'
	vector<DataType> tmp;
	unsigned int nbr;
	f>>nbr;
	if(!f.good()) {f.close();return false;}
	for(unsigned int i=0;i<nbr;i++)
	{
		DataType cell;
		f>>cell;
		if(f.good()) tmp.push_back(cell);
		else {f.close();return false;}
	}
	char end;
	f>>end;
	if(!f.good()||end!='$') {f.close();return false;}//episode end indicator
	data = tmp;
	f.close();
	return true;
}
template bool Read(string, unsigned int, vector<float> &);
template bool Read(string, unsigned int, vector<unsigned int> &);


template <class DataType>
bool Read(string filename,unsigned int N,vector< vector<DataType> > &data)
{
	ifstream f(filename.c_str(),ifstream::in);
	if(!f.is_open()) {f.close();return false;}
	while(N)
	{
		char tmp;
		tmp = f.get();//extract the next character until...
		if(!f.good()) {f.close();return false;}
		if(tmp=='#') N--;//...N-th episode begin indicator
	}
	f.ignore();//ignore the next '\n'
	data.clear();
	while(1)
	{
		vector<DataType> tmp;
		unsigned int nbr;
		f>>nbr;
		if(!f.good()) {f.close();return false;}
		char next;
		for(unsigned int i=0;i<nbr;i++)
		{
			DataType cell;
			f>>cell;
			if(f.good()) tmp.push_back(cell);
			else {f.close();return false;}
		}
		data.push_back(tmp);
		f.ignore();
		next = f.peek();
		if(next=='$') break;//episode end indicator
	}
	f.close();
	return true;
}
template bool Read(string, unsigned int, vector< vector<float> > &);
template bool Read(string, unsigned int, vector< vector<unsigned int> > &);
template bool Read(string, unsigned int, vector< vector<int> > &);
template bool Read(string, unsigned int, vector< vector<bool> > &);


