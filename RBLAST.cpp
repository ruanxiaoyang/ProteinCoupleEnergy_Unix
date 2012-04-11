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

int CPUNUM=core_count();
#include <display.h>
#include <ImpFile.h>
#include <PAM.h>
#include <GloAlign.h>
#include <EvoDist.h>
#include <RBLAST.h>
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
ImpBLOSUM("BLOSUM_score.txt",BLOSUMvec);
	
	char o,smflag='p',stflag='f',distype='s',supflag='n';
                //smflag scoring matrix flag:'p' for PAM and 'b' for BLOSUM
                //stflag scoring type flag:'f' for fixing scoring matrix and 'n' for not
                //gmaflag whether use gamma distance:'n' for not use,'y' for use
                //dt distance type:'s' for simple,'p' for poisson,'g' for gamma
                //supflag:suppress alignment step and only output pairwise distance
        char* seqfile=NULL;
	const char* dmpath="Pairwise_Block_Evo_Dist.txt";
        string rangestr;
       	int id=249,blosumid=61,gp=-1,rangeid,wlen=2;
        sarray<int> range;
        ifstream rangefile;
        while((o=getopt(argc,argv,"f:Bd:i:w:p:Dr:o:h"))!=-1)
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
                                distype=optarg[0];
                                if(distype!='s' && distype!='p' && distype!='g' && distype!='S' && distype!='P' && distype!='G')
                                {
                                        cout<<"Unknown distance estimation method specified. Use default \"simple\""<<endl;
                                        distype='s';
                                }
                                break;
			case 'i':id=stn<int>(optarg)-1;
                                break;
                        case 'w':wlen=stn<int>(optarg);
				if(wlen<=1 || wlen>=6)
				{	cout<<"Word length should between 2 and 5, auto set to 2\n";
					wlen=2;
				}
                                break;
                        case 'p':gp=stn<int>(optarg);
                                break;
			case 'D':supflag='y';
				break;
                        case 'r':rangestr=(string)(optarg);
                                rangefile.open(rangestr.c_str());
                                if(!rangefile)
                                {       range=parseid(rangestr);                //parse the index
                                }
                                else                                            //if it is a file
                                {       while(rangefile>>rangeid)
                                        {       range.pushback(rangeid);
                                        }
                                }
                                QuickSort(range,0,range.size()-1);
                                range=range.unique();
                                break;
			case 'o':dmpath=optarg;
				break;
			case 'h':system("less InstructionFiles/RBLAST");
				exit(0);
                        case '?':exit(0);
                }
        }
        if(seqfile==NULL)                       //error check for file existence
        {	cout<<"No sequence file specified"<<endl;
             	exit(0);
        }
	if(smflag=='p' && (id<0 || id>=400))    //error check for PAM scoring scheme
        {	cout<<"Need a value between 1 and 400 for PAM scoring"<<endl
			<<"auto set to PAM 250"<<endl;
                id=249;
        }
        else if(smflag=='b')                    //error check for BLOSUM scoring scheme
        {
                if(id<0 || id>=100)             //vector id validity
                {       cout<<"Need a value between 1 and 100 for BLOSUM scoring"<<endl
				<<"auto set to BLOSUM 62"<<endl;
			id=61;
                }
                if(BLOSUMvec[id-1].getrnum()==0)        //whether the convergence level is available
                {	cout<<"BLOSUM scoring not available for this level of convergence, check BLOSUM_score.txt for detail"<<endl
				<<"auto set to BLOSUM 62"<<endl;
			id=61;
                }
        }

        proformat prodb;
        if(range._size==0)                      //check if analyze specific sequences
                Custom(seqfile,prodb);
        else
                Custom(seqfile,prodb,range);

        if(prodb._seqnum<=1)
        {       cout<<"Only one sequence selected"<<endl;
                exit(0);
	}

        vector<darray<int> > alnmatxvec,Scoredb;
        if(smflag=='p')
                Scoredb=PAMvec;
        else if (smflag=='b')
                Scoredb=BLOSUMvec;
	
        sarray<int> alnscore,alnngscore;
	darray<double> distmatx,varmatx,smltmatx;
	RBLAST(prodb._multiseqid,wlen,Scoredb,id,alnmatxvec,alnscore,alnngscore,distmatx,varmatx,smltmatx,gp,smflag,distype);

        ofstream OP(dmpath);
        OP<<distmatx;                   //output pairwise distance matrix
        if(supflag=='y')                //stop if only output distance
                exit(0);

	string sm_str="PAM";
        if(smflag=='p')
		sm_str="PAM"+ntos(id+1);
        else if(smflag=='b')
                sm_str="BLOSUM"+ntos(id+1);

	
	cout<<"RBlast Glocal Alignment"<<endl 
        <<"pairwise align "<<prodb._seqnum<<" sequences\t"<<"Gap Penalty:"<<gp<<"\tScoring matrix:"<<sm_str<<endl<<endl; 
	int iter=0;
        for(int i=0;i<prodb._seqnum-1;++i) 
        { 
                for(int j=i+1;j<prodb._seqnum;++j) 
                {
                        cout<<prodb._infvec[i]._gi<<endl
                        <<prodb._infvec[j]._gi<<endl
                       	<<"Estimated Distance:PAM"<<(int)(distmatx(i,j)*100)<<"|"
                       	<<"Distype:"<<distype<<"|"
                       	<<"Variance:"<<(int)(varmatx(i,j))*100<<"|"
			<<"word length:"<<wlen<<"|"
                        <<"Score:"<<alnscore[iter]<<"|"
			<<"Nogap Score:"<<alnngscore[iter]<<endl;
                        Dsptwo(prodb._multiseqid[i],prodb._multiseqid[j],Scoredb[id],alnmatxvec[iter]);
                        cout<<endl;
                        iter+=1;
                }
        }

}
