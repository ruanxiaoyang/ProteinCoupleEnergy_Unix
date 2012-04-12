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
unordered_map<int,darray<int> > PAMTSthrd=readTSthrd('p'),BLOTSthrd=readTSthrd('b');
darray<double> Gonnet=GonnetScorematx();
#include <EvoTree.h>
int CPUNUM=core_count();
#include <display.h>
#include <ImpFile.h>
#include <PAM.h>
sarray<double> aafreq=getaafreq();

#include <MultiAlign.h>
#include <CoupleEnergy.h>
using namespace std;

int main(int argc,char* argv[])
{
	char* fpath=NULL;
	char o,type='c';
	int mask=5,trhd=5;
	while((o=getopt(argc,argv,"f:t:n:m:h"))!=-1)
	{	switch (o)
		{	case 'f':fpath=optarg;
				break;
			case 't':type=optarg[0];
				if(type != 's' && type != 'c' && type !='S' && type !='C')
				{	cout<<"Undefined coupling energe estimation type, autoset to binomial probability"<<endl;
					type='c';
				}
				break;
			case 'n':trhd=stn<int>(optarg);
				break;
			case 'm':mask=stn<int>(optarg);
				if(mask<0 || mask>9)
				{	cout<<"Mask need a value >=0 and <=9, autoset to 5"<<endl;
					mask=5;
				}
				break;
			case 'h':system("less InstructionFiles/CoupleEnergy");
				exit(0);
			case '?':exit(0);
		}
	}
	if(fpath==NULL)
	{	cout<<"Need .msa file"<<endl;
		exit(0);
	}
	vector<string> infvec;
	darray<int> combseq,posaacnt;
	ImpMSA(fpath,infvec,combseq);
	if(trhd<=3 || trhd>combseq.getrnum()/2)
	{	cout<<"Threshold needs a value between 3 and half of sequence number,autoset to 5"<<endl;
		trhd=5;
	}
	
	cout<<"Protein Coupling Sites Estimation"<<endl;
	cout<<"Estimation type:"<<type<<"|Seqeunce number Threshold:"<<trhd<<"|Mask:>="<<mask<<endl<<endl;
	DspMSA(combseq);
	vector<cpedstr> cpedresvec;
	if(type=='s')	
		SED(combseq,trhd,cpedresvec);
	else if(type=='c')
		CPED(combseq,trhd,cpedresvec);
	DspED(cpedresvec,combseq,mask);

}
