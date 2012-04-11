#include <iostream>
#include <sstream>

#include <cstring>
#include <unordered_map>
#include <fstream>
#include <vector>
#include <math.h>
#include <time.h>
#include <QuickSort.h>

#include <getdate.h>
#include <sarray.h>
#include <darray.h>
#include <LogScore.h>
#include <CommonFunctions.h>

int CPUNUM=core_count();
#include <display.h>
#include <ImpFile.h>
#include <GloAlign.h>
#include <PAM.h>

using namespace std;

int main()
{
darray<double> PAM1=getDayhoffPAM1();
sarray<double> aafreq=getaafreq();
vector<darray<double> > PAMdb, RelatedOddsdb;
vector<darray<int> > PAMvec,BLOSUMvec;
CRPAMdb(PAM1,PAMdb,400);
CRROdb(PAMdb,aafreq,RelatedOddsdb);
CRScoredb(RelatedOddsdb,10,PAMvec);
ImpBLOSUM("BLOSUM_score.txt",BLOSUMvec);

	sarray<int> lens(4);
	lens[0]=10;lens[1]=50;lens[2]=100;lens[3]=500;
	ofstream ofpam(".PAMTSthrd");
	int T,S;
	for(int p=0;p<400;p+=10)
	{
		for(int s=0;s<4;++s)
		{
			T_Sthld(lens[s],3,99,PAMvec[p],T,S);
			ofpam<<p<<" "<<lens[s]<<" "<<T<<" "<<S<<endl;
		}
	}
	ofstream ofblo(".BLOTSthrd");
	for(int b=0;b<100;++b)
	{
		if(!BLOSUMvec[b].empty())
		{for(int s=0;s<4;++s)
		{
			T_Sthld(lens[s],3,99,BLOSUMvec[b],T,S,'B');
			ofblo<<b<<" "<<lens[s]<<" "<<T<<" "<<S<<endl;
		}
		}
	}
}
