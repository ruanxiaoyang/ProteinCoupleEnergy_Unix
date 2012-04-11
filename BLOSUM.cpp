#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <cstdlib>
#include <unistd.h>
#include <vector>
#include <math.h>
#include <time.h>
#include <getdate.h>
#include <sarray.h>
#include <darray.h>
#include <LogScore.h>
#include <CommonFunctions.h>
#include <BLOSUM.h>

using namespace std;

int main(int argc,char* argv[])
{
	int o;
	char rflag='Y';
	char* fpath=NULL;
	char* opath=NULL;
	//c for using custom block file. s for suppressing relative frequency output. o for scoring matrix file name
	//h for help file
	while((o=getopt(argc,argv,"c:so:h"))!=-1) 
	{
		switch(o)
		{
			case 'c':fpath=optarg;
				cout<<"use custom block file"<<endl;
				break;
			case 's':rflag='N';
				cout<<"surpress relative frequency output"<<endl;
				break;
			case 'o':
				opath=optarg;
				break;
			case 'h':
				system("less InstructionFiles/BLOSUM");
				exit(0);
			case '?':exit(0);
		}
	}
	CreateBlosumScoreMatrix(fpath,rflag,opath);
}
