#include <sstream>
//Get CPU number
int core_count()
{
	system("grep \"physical id\" /proc/cpuinfo | wc -l > ./.cpucount.txt");
	ifstream INPUT(".cpucount.txt");
	string s;
	getline(INPUT,s);
	istringstream ss(s);
	int n;
	ss>>n;
	return n;
}

//Control format.Same as progressbar.
//Take 2 values."_i"usually the start number of iteration."_size"usually the end number of iteration."_acry" round to _acry digits after decimal point
inline void PctMarker(const int & _i,const int & _size,const int & _acry=0)
{
	if(_i>_size)
		return;
	double per;
	if(_size==0)
		per=1;
	else 
		per=(double)_i/_size;
	int tmp=(int)(per*pow(10.0,(double)(2+_acry)))*10;
	cout.flags(ios::fixed);
	cout.width(3+_acry);
	cout.precision(_acry);
	if(tmp%10==0)
	{
			cout<<tmp/pow(10.0,(double)(_acry+1))<<"%";
			if(_i!=_size)
			{
				for(int i=0;i<4+_acry;++i)
					cout<<"\b";
			}
	}
}
//Break string obtain from getline() function into blank ' 'seperated words and store these words in vector<sarray<char>>
//Take 1 value. "_tmpstr" string type
//Return 1 value.vector<sarray<char>> containing word of "_tempstr"
vector<sarray<char> > _Breakstring(string & _tmpstr)
{
	vector<sarray<char> > wordvec;
	sarray<char> word;
	for(int i=0;i<_tmpstr.size();++i)
	{
		if(_tmpstr[i]!=' ' && _tmpstr[i]!='\t')
		{
			word.pushback(_tmpstr[i]);
			continue;
		}
		if(_tmpstr[i]==' ' && word.size()==0)
			continue;
		if((_tmpstr[i]==' ' || _tmpstr[i]=='\t') && word.size()>0)
		{
			wordvec.push_back(word);
			word.clear();
		}
	}
	if(word.size()>0)
		wordvec.push_back(word);
	return wordvec;
}

//compute length of sequence after alignment(gap+original length)
//take 3 values."_matrix" standard 0 or 1 matrix."_rownum"/"_colnum" row and column of "_matrix"
//return 1 value.length of sequence after alignment
int Length(darray<int> & _matrix,const int & _rownum,const int & _colnum)
{
        int previ=-1,prevj=-1,gap=0;
        for(int i=0;i<_rownum;++i)
        {
                for(int j=prevj+1;j<_colnum;++j)
                {
                        if(_matrix(i,j)==1)
                        {
                                gap+=i-previ-1;                                //compute row gap 
                                previ=i;prevj=j;
                                if(prevj==_colnum-1)
                                {
                                        gap+=_rownum-1-i;
                                        return gap+_colnum;
                                }
                                break;
                        }
                }
        }
        return gap+_colnum;
}

//Control format.Output number mark for sequence.Spaced by 10 AA/nucleotide.
//Take 2 values."_size"length of sequence."_os"output file stream
inline void LandMarker(const int & _size,ostream & _os)
{
        int hang=0;
        for(int i=0;i<_size;++i)
        {
                if(i>=10 && i%10==0)
                {
                        _os<<i;
                    hang=(int)(log((double)i)/log(10.0));
                        continue;
                }
                if(hang==0)
            _os<<" ";
                else if(hang>0)
                        hang-=1;
        }
}

inline void MileStone(const int & _size,ostream & _os)
{
        for(int i=0;i<_size;++i)
        {
                if(i%10==0)
                        _os<<"+";
                else
                        _os<<"-";
        }
}

template <class T>
string ntos(const T& num)
{
	std::ostringstream oss;
	oss<<num;
	return oss.str();
}

template <class T>
T stn(const std::string& s, std::ios_base& (*f)(std::ios_base&)=std::dec)	//*f is pointer to function which takes ios_base and return ios_base
{
	T t;
	std::istringstream iss(s);
	if(!(iss >> f >> t).fail())
	{
		return t;
	}
	else
	{
		cout<<"error"<<endl;
		return 0;
	}
}

sarray<int> parseid(string & str)
{
	sarray<int> res;
	sarray<char> mark;
	string stream;
	char tmp;
	if(str[0]==',')
		str=str.substr(1,str.size()-1);
	if(str[str.size()-1]!=',')
		str.push_back(',');
	for(int i=0;i<str.size();++i)
	{ 	
		if(str[i]==',' || str[i]=='-')
		{stream.push_back(' ');
			if(str[i]==',')
				mark.pushback('n');
			else if(str[i]=='-')
				mark.pushback('y');
		}
		else
			stream.push_back(str[i]);
	}
	std::istringstream s(stream);
	int id,iter=0,end;
	while(s>>id)
	{
		if(mark[iter]=='n')
		{
			res.pushback(id);
			iter+=1;
		}
		else if(mark[iter]=='y')
		{
			s>>end;
			for(int i=id;i<=end;++i)
				res.pushback(i);
			iter+=2;
		}
	}
	return res;
}


//generate all combinations that needs to be done for thread function
//take 1 value. "_sqn" the number of sequence
darray<int> generatethreadid(const int _sqn)
{
	int iter=0;
	darray<int> thrdtaskmatx(2,(int)((double)_sqn*(double)(_sqn-1)/2.0));
        for(int i=0;i<_sqn;++i) 
        {
                for(int j=i+1;j<_sqn;++j)
                {
                        thrdtaskmatx(0,iter)=i;
                        thrdtaskmatx(1,iter)=j;
                        iter+=1;
                }
        }
	return thrdtaskmatx;
}

 //Dayhoff/Blosum100 AA accumulate frequency data.Used for creating random AAid sequence
sarray<int> accumAAfreq(const char _choice='D')   
{
        sarray<int> accumfreq(20,0);
        if(_choice=='D' || _choice=='d'||_choice=='G'||_choice=='g')
        {
            accumfreq[0]=87;accumfreq[1]=128;accumfreq[2]=168;accumfreq[3]=215;accumfreq[4]=248;accumfreq[5]=286;accumfreq[6]=336;accumfreq[7]=425;accumfreq[8]=459;accumfreq[9]=496;accumfreq[10]=581;
            accumfreq[11]=662;accumfreq[12]=677;accumfreq[13]=717;accumfreq[14]=768;accumfreq[15]=838;accumfreq[16]=896;accumfreq[17]=906;accumfreq[18]=936;accumfreq[19]=1000;
        }
        else if(_choice=='B' || _choice=='b')
        {
                accumfreq[0]=73;accumfreq[1]=122;accumfreq[2]=164;accumfreq[3]=218;accumfreq[4]=247;accumfreq[5]=281;accumfreq[6]=338;accumfreq[7]=419;accumfreq[8]=444;accumfreq[9]=509;accumfreq[10]=604;
                accumfreq[11]=656;accumfreq[12]=681;accumfreq[13]=727;accumfreq[14]=766;accumfreq[15]=824;accumfreq[16]=877;accumfreq[17]=892;accumfreq[18]=927;accumfreq[19]=1000;
        }
        return accumfreq;
}

//generate random AA id sequence.optimized for 20 kinds of AA.it costs averagely (n^2+n+20)/2n time to generate one AA.when n=4 it has the smallest value (5 computer clock) 
void randAAidseq(const int & _length,sarray<int> & _accumdb,sarray<int> & randseqid)    
{
	int temp;
	randseqid.resize(_length);
    	for(int n=0;n<_length;++n)
	{
		temp=rand()%1001;
		if(temp<=_accumdb[4])
		{
			for(int i=0;i<=4;++i)
	        {
		       if(temp<=_accumdb[i])
		       {
				   randseqid[n]=i;
			      break;
		       }
	        }
		}
		else if(temp<=_accumdb[9])
		{
			for(int i=5;i<=9;++i)
			{
			   if(temp<=_accumdb[i])
		       {
			      randseqid[n]=i;
			      break;
		       }
			}
		}
		else if(temp<=_accumdb[14])
		{
			for(int i=10;i<=14;++i)
			{
			   if(temp<=_accumdb[i])
		       {
			      randseqid[n]=i;
			      break;
		       }
			}
		}
		else
		{
			for(int i=15;i<=19;++i)
			{
			   if(temp<=_accumdb[i])
		       {
			      randseqid[n]=i;
			      break;
		       }
			}
		}
	}
}

unordered_map<int,darray<int> > readTSthrd(const char & type)
{
	unordered_map<int, darray<int> > res;
	int id,l,t,s;
	sarray<int> tmp(3);
	ifstream IN;
	if(type=='p')
		IN.open(".PAMTSthrd");
	else
		IN.open(".BLOTSthrd");
	if(!IN)
	{	cout<<"Threshold file not found, run TS";
		exit(0);
	}

	while(IN>>id && IN>>l && IN>>t && IN>>s)
	{
		tmp[0]=l;tmp[1]=t;tmp[2]=s;
		res[id].push_row(tmp);
	}
	return res;
}

//Menu:[RBLAST]-[BLAST]-[T_Sthld]
//Calculate t and s score distribution for two random sequences.
//Take 5 values."_seqlen" sequence length of the shorter sequence."_wlen" word length."_cutoffper" cutoff percentile."_scmatx" score matrix used."_DorB" Dayhoff or Blosum AA frequency database 
//Compute 2 values."tthld" t score threshold."sthld" s score threshold.
template <class type>  
void T_Sthld(const int & _seqlen,const int & _wlen,const int & _cutoffper,darray<type> & _scmatx,type & tthld,type & sthld,char _DorB='D')
{
	int num=3000,mark=0;					//"num" number of replication
	type tscore,sscore,tempscore;
	sarray<int> seqa,seqb;
	sarray<type> scorevec(num,(type)0),temparr;
	sarray<int> accumdb;
	if(_DorB=='D')
		accumdb=accumAAfreq('D');			//accumulate AA frequency.
	else
		accumdb=accumAAfreq('B');	
	for(int n=0;n<num;++n)
	{
		tscore=(type)0;
		randAAidseq(_wlen,accumdb,seqa);		//generate random AA seq with "_wlen" length
		randAAidseq(_wlen,accumdb,seqb);
		for(int w=0;w<_wlen;++w)
			tscore+=_scmatx(seqa[w],seqb[w]);
		scorevec[n]=tscore;
	}
	QuickSort(scorevec,0,num-1);				//sort "scorevec" from min to max
	tthld=scorevec[(int)((double)num*(double)_cutoffper/100.0)-1];    
	for(int n=0;n<num;++n)
	{
		randAAidseq(_seqlen,accumdb,seqa);
		randAAidseq(_seqlen,accumdb,seqb);
		sscore=(type)INT_MIN;				//initiate "sscore" to a very small value
		for(int i=0;i<_seqlen;++i)
		{
			if(_scmatx(seqa[i],seqb[i])>=0)		//only start pushback when meet a positive score
			{
				mark=0;
				temparr.pushback(_scmatx(seqa[i],seqb[i]));
			}
			else if(_scmatx(seqa[i],seqb[i])<0 && temparr.size()>0)                    
			{
				mark+=1;
				temparr.pushback(_scmatx(seqa[i],seqb[i]));
			}
			if(temparr.size()>0 && (mark==_wlen || i==_seqlen-1))
			{
				tempscore=(type)0;
				temparr.resize(temparr.size()-mark);
				for(int j=0;j<temparr.size();++j)
					tempscore+=temparr[j];
		        if(tempscore>sscore)	//find the highest score in each pair of random alignment and store in "sscore"
				sscore=tempscore;
				temparr.clear();
			}
		}
		scorevec[n]=sscore;              
	}
	QuickSort(scorevec,0,num-1);
	sthld=scorevec[(int)((double)num*(double)_cutoffper/100.0)-1];
}
