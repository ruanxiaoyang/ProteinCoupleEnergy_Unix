#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <getdate.h>
#include <sarray.h>
#include <darray.h>
#include <LogScore.h>
#include <CommonFunctions.h>
#include <GloAlign.h>
#include <ImpFile.h>

int CPUNUM=core_count();
//#include <EvoDist.h>
#include <PAM.h>
using namespace std;

int main()
{
	darray<double> PAM1=getDayhoffPAM1();
	sarray<double> aafreq=getaafreq();
	vector<darray<double> > PAMdb, RelatedOddsdb;
	vector<darray<int> > Scoredb;
	CRPAMdb(PAM1,PAMdb,400);
	CRROdb(PAMdb,aafreq,RelatedOddsdb);
	CRScoredb(RelatedOddsdb,10,Scoredb);
	proformat prodb;	
	Custom("./proseqdb.txt",prodb);
	
	darray<double> distmatx,varmatx;
	char distype='s';
//	AprxDistMatx(prodb._multiseqid,distype,distmatx,varmatx);
//	cout<<distmatx<<endl;
//	cout<<varmatx<<endl;
	
}
