#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
using namespace std;
int main () {
  string line;
  ifstream myfile ("../../Model/Test/inp.in");
  vector< vector<int> > data;
  if (myfile.is_open())
    {
      while ( getline (myfile,line) )
	{
	  //data.push_back(atoi(line.c_str()));
	  //cout << line.c_str() << endl;
	  stringstream stream(line);
	  int n;
	  vector<int> tmpdata;
	  while(stream >> n){
	    tmpdata.push_back(n);
	  }
	  data.push_back(tmpdata);
	}
      myfile.close();
    }

  // See data here
  for (vector<vector<int> >::iterator it = data.begin() ; it != data.end(); ++it)
    {
      vector<int> lineVec = it[0];
      for (int i = 0; i < lineVec.size(); i++)
      {
        cout << lineVec[i] << " ";
      }
      cout << endl;
    }
  return 0;
}
