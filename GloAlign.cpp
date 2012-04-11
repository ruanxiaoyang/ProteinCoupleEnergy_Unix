#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <math.h>
#include <time.h>
#include <QuickSort.h>
#include <pthread.h>
#include <getdate.h>
#include <sarray.h>
#include <darray.h>
#include <LogScore.h>
#include <CommonFunctions.h>

int CPUNUM=core_count();
#include <display.h>
#include <ImpFile.h>
#include <GloAlign.h>

#include <EvoDist.h>
#include <PAM.h>
using namespace std;

int main(int argc,char* argv[])
{
darray<double> PAM1=getDayhoffPAM1();
sarray<double> aafreq=getaafreq();
vector<darray<double> > PAMdb, RelatedOddsdb;
vector<darray<int> > PAMvec,BLOSUMvec;
CRPAMdb(PAM1,PAMdb,400);
CRROdb(PAMdb,aafreq,RelatedOddsdb);
CRScoredb(RelatedOddsdb,10,PAMvec);


	char o,smflag='p',stflag='f',gmaflag='n',dt='s',supflag='n';
		//smflag scoring matrix flag:'p' for PAM and 'b' for BLOSUM
		//stflag scoring type flag:'f' for fixing scoring matrix and 'n' for not
		//gmaflag whether use gamma distance:'n' for not use,'y' for use
		//dt distance type:'s' for simple,'p' for poisson,'g' for gamma
		//supflag:suppress alignment step and only output pairwise distance 
	char* seqfile=NULL;
	const char* dmpath="Pairwise_Evo_Dist.txt";
	string rangestr;
	int id=250,blosumid=62,gp=-1,gma=1,rangeid;
	sarray<int> range;
	ifstream rangefile;
	while((o=getopt(argc,argv,"f:Bd:i:Yp:g:Dr:o:h"))!=-1)
	{
		switch(o)
		{
			case 'f':
				seqfile=optarg;
				break;
			case 'B':smflag='b';
				id=blosumid;
				ImpBLOSUM("BLOSUM_score.txt",BLOSUMvec);
				break;
			case 'd':
				dt=optarg[0];
				if(dt!='s' && dt!='p' && dt!='g' && dt!='S' && dt!='P' && dt!='G')
				{
					cout<<"Unknown distance estimation method specified. Use default \"simple\""<<endl;
					dt='s';
				}
				break;
			case 'i':id=stn<int>(optarg);		
				break;
			case 'Y':
				stflag='n';
				break;
			case 'p':gp=stn<int>(optarg);
				break;
			case 'g':gmaflag='y'; 
				gma=stn<int>(optarg);
				break;
			case 'D':supflag='y';				//suppress alignment step
				break;
			case 'r':rangestr=(string)(optarg);		
				rangefile.open(rangestr.c_str());
				if(!rangefile)	
				{	range=parseid(rangestr);		//parse the index
				}
				else						//if it is a file
				{	while(rangefile>>rangeid)
					{	range.pushback(rangeid);
					}
				}
				QuickSort(range,0,range.size()-1);
				range=range.unique();
				break;
			case 'o':dmpath=optarg;
				break;
			case 'h':system("less InstructionFiles/GloAlign");
				exit(0);
			case '?':exit(0);
		}
	}

	if(seqfile==NULL)			//error check for file existence
	{
		cout<<"No sequence file specified"<<endl;
		exit(0);
	}
	if(smflag=='p' && (id<=0 || id>400))	//error check for PAM scoring scheme
	{
		cout<<"Need a value between 1 and 400 for PAM scoring"<<endl;
		exit(0);
	}
	else if(smflag=='b')			//error check for BLOSUM scoring scheme
	{
		if(id<=0 || id>100)		//vector id validity
		{	cout<<"Need a value between 1 and 100 for BLOSUM scoring"<<endl;
			exit(0);
		}
		if(BLOSUMvec[id-1].getrnum()==0)	//whether the convergence level is available
		{
			cout<<"BLOSUM scoring not available for this level of convergence, check BLOSUM_score.txt for detail"<<endl;
			exit(0);
		}
		if(stflag=='n')			//need to use fixed scoring matrix
		{	cout<<"Blosum does not support dynamic scoring"<<endl;
			exit(0);
		}
	}
	if(dt!='g' && dt!='G' && gmaflag=='y')	//error check for gamma parameter
		cout<<"Not using gamma distance estimation, -g option skipped"<<endl;

	proformat prodb;
	
	if(range._size==0)			//check if analyze specific sequences	
		Custom(seqfile,prodb);
	else
		Custom(seqfile,prodb,range);
	if(prodb._seqnum<=1)
	{	cout<<"Only one sequence selected"<<endl;
		exit(0);
	}
	darray<double> distmatx,varmatx;
	vector<darray<int> > BKMATXvec;
	char distype='s';
	AprxDistMatx(prodb._multiseqid,dt,distmatx,varmatx,BKMATXvec,gma);	//pairwise distance

	ofstream OP(dmpath);
	OP<<distmatx;			//output pairwise distance matrix
	if(supflag=='y')		//stop if only output distance
		exit(0);

	vector<darray<int> > alnmatxvec,Scoredb;
	sarray<int> alnscore;
	if(smflag=='p')
		Scoredb=PAMvec;
	else if (smflag=='b')
		Scoredb=BLOSUMvec;
		
	NWwp(BKMATXvec,prodb._multiseqid,Scoredb,distmatx,gp,alnmatxvec,alnscore,stflag,id);	//pairwise waypoint

	string sm_str="PAM";
	if(smflag=='p')
	{
		if(stflag=='f')
			sm_str="PAM"+ntos(id);
		else if(stflag=='n')
			sm_str="PAM1-400";
	}
	else if(smflag=='b')
		sm_str="BLOSUM"+ntos(id);
	
	int iter=0;
	
	cout<<"Needleman Global Alignment"<<endl
	<<"pairwise align "<<prodb._multiseqid.size()<<" sequences|"<<"Gap Penalty:"<<gp<<"|Scoring matrix:"<<sm_str<<endl<<endl;
	for(int i=0;i<prodb._seqnum-1;++i)
	{
		for(int j=i+1;j<prodb._seqnum;++j)
		{
			cout<<prodb._infvec[i]._gi<<endl
			<<prodb._infvec[j]._gi<<endl
			<<"Estimated Distance:PAM"<<(int)(distmatx(i,j)*100)
			<<"\tDistance type:"<<dt
			<<"\tVariance:"<<(int)(varmatx(i,j))*100
			<<"\tScore:"<<alnscore[iter]<<endl;
			Dsptwo(prodb._multiseqid[i],prodb._multiseqid[j],Scoredb[id-1],alnmatxvec[iter]);
			cout<<endl;
			iter+=1;
		}
	}
	
}
