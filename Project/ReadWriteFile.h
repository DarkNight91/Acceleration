template <class DataType>
bool Store(string filename, DataType data);

template <class DataType>
bool Store(string filename, vector<DataType> data);

template <class DataType>
bool Store(string filename, vector< vector<DataType> > data);

template <class DataType>
bool Read(string filename,unsigned int N,DataType &data);

template <class DataType>
bool Read(string filename,unsigned int N,vector<DataType> &data);

template <class DataType>
bool Read(string filename,unsigned int N,vector< vector<DataType> > &data);

