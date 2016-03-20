#include <iostream>
#include <string>

using namespace std;

int main()
{
    std::cout<<"Hello C++!" <<endl;
    string  server_url = getenv("SERVER_URL"); 
    std::cout<<"get from environment = "<<server_url<<endl;
    return 1;
}

