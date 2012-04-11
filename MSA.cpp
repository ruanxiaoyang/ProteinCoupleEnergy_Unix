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
#include <GloAlign.h>
#include <EvoDist.h>
#include <MultiAlign.h>
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
	
	char o,smflag='p',distype='s',ordered='n',iteration='y',seqname='y';
                //smflag scoring matrix flag:'p' for PAM and 'b' for BLOSUM
                //dt distance type:'s' for simple,'p' for poisson,'g' for gamma
		//ordered whether output MSA result as its order in original file
        char* seqfile=NULL;
	const char * MSArespath="Res.msa";
        string rangestr;
        int id=249,blosumid=61,gp=-3,rangeid,wlen=2,anchor=80;
        sarray<int> range;
        ifstream rangefile;
        while((o=getopt(argc,argv,"f:Bd:i:w:p:Ia:r:ODo:h"))!=-1)
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
			case 'I':iteration='n';
				break;
			case 'a':anchor=stn<int>(optarg);
				if(anchor<=0 || anchor >100)
				{	cout<<"Anchor point value should be >=0 and <100, auto set to 80\n";
					anchor=80;
				}
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
			case 'O':ordered='y';
				break;
			case 'D':seqname='n';
				break;
			case 'o':MSArespath=optarg;
				break;
			case 'h':system("less InstructionFiles/MSA");
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

        if(prodb._seqnum<=2)
        {       cout<<"At least 3 sequences are required."<<endl;
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

	const char* dmpath="Pairwise_Block_Evo_Dist_MSA.txt";
        ofstream OP(dmpath);
        OP<<distmatx;                   //output pairwise distance matrix

	nbstr evotree;
	NJtree(distmatx,evotree);
	darray<int> combseq;
	_Alignseqs(prodb._multiseqid,alnmatxvec,evotree,combseq,gp,iteration,anchor);
	double score=_MSAsc(combseq,gp);
	
	string sm_str="PAM";
        if(smflag=='p')
                sm_str="PAM"+ntos(id+1);
        else if(smflag=='b')
                sm_str="BLOSUM"+ntos(id+1);

	if(seqname=='y' && ordered=='n')	//order sequence names
		reorderseqname(prodb,evotree);	

	if(ordered=='y')			//order sequences
		reorder(combseq,evotree);

	////////////////////////////////////////////// output to screen
	cout<<"Multiple Sequences Alignment"<<endl
        <<"Sequence number:"<<prodb._seqnum<<"|Gap Penalty:"<<gp<<"|Scoring matrix:"<<sm_str<<"|distance type:"<<distype<<endl
	<<"Word Length:"<<wlen<<"|Score(Gaston):"<<score<<endl<<endl;
	if(seqname=='y')
		DspMSAseqname(prodb._infvec);

	DspMSA(combseq);

	ofstream OM(MSArespath);
	OM<<combseq.getrnum()<<" "<<combseq.getcnum()<<endl;
	DspMSAseqname(prodb._infvec,OM);
	OM<<combseq;	
}
