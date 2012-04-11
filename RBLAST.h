struct WDlibstr
{
	darray<int> word;
	darray<int> location;
	int size;
	int seqlen;
};

struct getwdlibstr
{
	vector<sarray<int> > * seqdb;
	vector<WDlibstr> * WDlibvec;
	int wlen;
	int thrdid;
};

typedef unordered_map<string,int> wmap;

void CRlib(sarray<int> & _seqid,int & _wlen,WDlibstr & WDlib)
{
	WDlib.word.clear();
	WDlib.location.clear();
        WDlib.size=0;
	sarray<int> word(_wlen,0);
	string str;
	int id=0;
	wmap tmpmap;
        for(int i=0;i<_seqid.size()-_wlen+1;++i)
        {
                for(int w=0;w<_wlen;++w)
                        word[w]=_seqid[i+w];
		str=word.tos();
                if((id=tmpmap[str])==0)
		{
			WDlib.word.push_row(word);
			WDlib.location.push_row(1);
			WDlib.location(WDlib.location.getrnum()-1,0)=i;
			WDlib.size+=1;
			tmpmap[str]=WDlib.size;
		}
                else
                        WDlib.location.push_to_row(id-1,i);
        }
	WDlib.seqlen=_seqid.size();
}

void* THCRlib(void* param)
{
	getwdlibstr* THdt=static_cast<getwdlibstr*>(param);
	int thrdid=THdt->thrdid;
	int start,end,size,quantum;		//"start""end"thread task start and end position."quantum" aliquot of task
        size=THdt->seqdb->size();
        quantum=(int)((double)size/CPUNUM);
        if(quantum<1)			//if CPU number is larger than task number,set "quantum" to 1,otherwise it is 0.
                quantum=1;
        start=thrdid*quantum;
        if(thrdid!=CPUNUM-1)
        end=(thrdid+1)*quantum-1;
        else
                end=size-1;
        if(end>=size)
                return NULL;
	for(int i=start;i<=end;++i)
		CRlib(THdt->seqdb->operator[](i),THdt->wlen,THdt->WDlibvec->operator[](i));
	
}

void CRlib(vector<sarray<int> > & _seqdb,const int & _wlen,vector<WDlibstr> & WDlibvec)
{
	WDlibvec.resize(_seqdb.size());
	vector<getwdlibstr> dtstr(CPUNUM);
	pthread_t threads[CPUNUM];
	for(int i=0;i<CPUNUM;++i)
	{
		dtstr[i]={&_seqdb,&WDlibvec,_wlen,i};
		pthread_create(&threads[i],NULL,THCRlib,(void*)&dtstr[i]);
	}
	for(int i=0;i<CPUNUM;++i)
		pthread_join(threads[i],NULL);
}

void fetchTS(const int & _id,const int & _len,int & T,int & S,const char & _type='p')
{
	int id=((int)(_id/10))*10;
	int rown=0;
	if(_len>50 && _len<=100)
		rown=1;
	else if(_len>100 && _len<=500)
		rown=2;
	else if(_len>500)
		rown=3;
	darray<int> * pd;
	if(_type=='p')
		pd=&(PAMTSthrd[id]);
	else
		pd=&(BLOTSthrd[id]);
	T=pd->operator()(rown,1);
	S=pd->operator()(rown,2);
}

template<class type>
void BLAST(WDlibstr & _qlib,WDlibstr & _tlib,sarray<int> & _qidseq,sarray<int> & _tidseq,vector<darray<type> >& _scmatxvec,const int & _scmatxid,const int & _wlen,darray<int> & alnmatx,const char & _sctype)
{
	type tscore;
	int T,S,rownum=_qlib.seqlen,colnum=_tlib.seqlen,seqlen=(rownum>colnum ? colnum:rownum);
	darray<type> tempmatx(rownum,colnum,(type)0);		//its necessary to initialize to 0 (for _Converge function)
	fetchTS(_scmatxid,seqlen,T,S,_sctype);
//cout<<"T "<<T<<" S "<<S<<endl;
	for(int i=0;i<_qlib.size;++i)
	{	for(int j=0;j<_tlib.size;++j)
		{	tscore=(type)0;
			for(int w=0;w<_wlen;++w)
				tscore+=_scmatxvec[_scmatxid](_qlib.word(i,w),_tlib.word(j,w));	//calculate score of hit
			if(tscore>=T)			//if "tscore" reached threshold then extend hit
			{ 	for(int k=0;k<_qlib.location.getcnum(i);++k)
				{	for(int m=0;m<_tlib.location.getcnum(j);++m)
				_Extend(_wlen,S,_qlib.location(i,k),_tlib.location(j,m),_scmatxvec[_scmatxid],_qidseq,_tidseq,tempmatx);
				}
			}
		}
	}
//cout<<tempmatx<<endl;
	_Converge(tempmatx,rownum,colnum);			//converge way point to get unique alignment pathway
	alnmatx.fast_resize(rownum,colnum);
	for(int i=0;i<rownum;++i)				//retore alignment alnmatx to standard 0 or 1 form
	{	for(int j=0;j<colnum;++j)
		{	if(tempmatx(i,j)==0 || tempmatx(i,j)==INT_MIN)
				alnmatx(i,j)=0;
			else
				alnmatx(i,j)=1;
		}
	}
}

template<class type>
void _Extend(const int & _wlen,const double & _S,const int & _qst,const int & _tst,darray<type> & _scmatx,sarray<int> & _qidseq,sarray<int> & _tidseq,darray<type> & matrix)
{
	type score;
	int checker=0,qery,bseq,tempqst=_qst,temptst=_tst,qst,tst,right=0,left=0;
	//"tempqst"/"temptst" record the start position.
	//"left"/"right"record the left/right terminal where "_wlen" mismatch occur
	while(tempqst<_qidseq.size()-1 && temptst<_tidseq.size()-1 && _qidseq[tempqst]!=_tidseq[temptst])   
			//if start position does not equal, move to right until it is equal
	{	tempqst+=1;
		temptst+=1;
	}
	if(matrix(tempqst,temptst)!=0)	//stop if this position is covered already
		return;                                          
	qst=tempqst;tst=temptst;int mark=tempqst;    
	while(qst>=0 && tst>=0)		//go until "_wlen" mismatch occur or reach the left terminal
	{	if(_qidseq[qst]!=_tidseq[tst])
			checker+=1;
		else
		{	mark=qst;
			checker=0;
		}
		if(checker==_wlen)
			break;
		qst-=1;tst-=1;
	}
	left=mark-tempqst;
	qst=tempqst+1;tst=temptst+1;checker=0;mark=tempqst;		//restore iterator to 1 pass start position
	while(qst<_qidseq.size() && tst<_tidseq.size())	//go forward until "_wlen" mismatch or reach the right terminal
	{	if(_qidseq[qst]!=_tidseq[tst])
			checker+=1;
		else
		{	mark=qst;
			checker=0;
		}
		if(checker==_wlen)
			break;
		qst+=1;tst+=1;
	}
	right=mark-tempqst;
	score=(type)0;
        for(int i=left;i<=right;++i)
		score+=_scmatx(_qidseq[tempqst+i],_tidseq[temptst+i]);
	if(score>=_S)			//if alignment score is larger than threshold then update match matrix information
	{	for(int i=left;i<=right;++i)
			matrix(tempqst+i,temptst+i)=score;
	}
}

template <class type>
void _Converge(darray<type> & _matrix,const int & _rownum,const int & _colnum)
{
	type comp1,comp2;
	int topI,topJ,botI,botJ,topi,topj,boti,botj;
	for(int i=0;i<_rownum;++i)                           
	{	for(int j=_colnum-1;j>=0;--j)		//iterate from right to left                              
		{	if(_matrix(i,j)==0)		//if start point "_matrix(i,j)" equal to 0,move to next column
				continue;
			else if(_matrix(i,j)==INT_MIN)  //if start point is row stop marker INT_MIN,move to next row
				break;
			else
			{
			comp1=_matrix(i,j);
			_Terminal(_matrix,i,j,topi,topj,boti,botj);
	            	for(int ii=i;ii<_rownum;++ii)
	            	{   for(int jj=j;jj>=0;--jj)
		            {	if(ii==i && jj==j)		//avoid compare to itself
					continue;
				if((comp2=_matrix(ii,jj))==INT_MIN)	//stop if stop marker appear   
					break;
				if(comp2!=0)			//erase the part of weaker alignment
				{
					_Terminal(_matrix,ii,jj,topI,topJ,botI,botJ);
					if(comp2>comp1 || (comp2==comp1 && abs(ii-jj)<=abs(i-j)))   //if the newly found segment is stronger than the target
					{	for(int m=0;m<=boti-topi;++m)
				                {	if((topi+m)<=botI && (topj+m)>=topJ)
						                 _matrix(topi+m,topj+m)=0;
				                }
					}
					if(comp2<comp1 ||  (comp2==comp1 && abs(ii-jj)>abs(i-j)))
			            	{	for(int m=0;m<=botI-topI;++m)
						{
					            	if((topI+m)>=topi && (topJ+m)<=botj)
						        	    _matrix(topI+m,topJ+m)=0;
						}
					}
					if(_matrix(i,j)==0)//only when start point has been erased can move to next column
	                			  goto nextj;
				}
			   }
			}
			for(int ii=i;ii<_rownum;++ii)//if start point was preserved (!=0) after convergence,mark lower left area of start point as useless
	            	{	for(int jj=j;jj>=0;--jj)
				{
					if(ii==i && jj==j)
						continue;
					if(_matrix(ii,jj)==INT_MIN)
						break;
					else
						_matrix(ii,jj)=INT_MIN;
				}
			}
		}
		
nextj:;
		}
	}
}

template<class type>
void _Terminal(darray<type> & _matrix,const int & _i,const int & _j,int & topi,int & topj,int & boti,int & botj)
{
	int i=_i-1,j=_j-1;
	while(i>=0 && j>=0 && _matrix(i,j)>0)		//go backward to find top terminal
	{
			i-=1;j-=1;
	}
	topi=i+1;topj=j+1;
	i=_i+1;j=_j+1;
	while(i<_matrix.getrnum() && j<_matrix.getcnum() && _matrix(i,j)>0)                //go forward to find bottom terminal
	{
			i+=1;j+=1;
	}
    boti=i-1,botj=j-1;
}

//Menu:[RBLAST]-[_Blankmark]
//Find out blank areas not covered by word extending strategy
//Take 3 values."alnmatx" alignment matrix (standard 0 or 1 form)."_rownum"/"_colnum" row/column number of "alnmatx"
//Teturn 4-column darry<int> containing coordinates marking the blank area.column 0,1,2,3 are top i,top j,bottom i,bottom j respectively
void _Blankmark(darray<int> & alnmatx,const int & _rownum,const int & _colnum,darray<int> & blank)
{
	bool open=false;
	int i=0,j=0,ii,jj,lr;                                //"lr" last row id
    	while(i<_rownum && j<_colnum)
	{
	    for(jj=j;jj<_colnum;++jj)                        //scan the row for 1
		{
			if(alnmatx(i,jj)!=0)                         //if record found,record the column id and break
			{
				j=jj;
				break;
			}
		}
		if(jj==_colnum)                                  //only start scanning the column when no record found in the row
		{
			for(ii=i;ii<_rownum;++ii)
			{
				if(alnmatx(ii,j)!=0)                     //if record found,record the row id and break
				{
					i=ii;
					break;
				}
			}
		}
		if(jj==_colnum && ii==_rownum && open==false)    //when no record found in row and column and blank has not been openned yet,open a blank
		{
			open=true;
			blank.push_row(4);
			lr=blank.getrnum()-1;
			blank(lr,0)=i,blank(lr,1)=j;
		}
		else if((jj!=_colnum || ii!=_rownum) && open==true) //if a blank is open,and a record was found,close the blank and store the end coordinate
		{
				blank(lr,2)=i-1;
				blank(lr,3)=j-1;
				open=false;
		}
		i+=1;j+=1;
	}
	if(open==true)		//if the lower right terminal of "alnmatx" is 0,filled in the coordinate 
	{
		blank(lr,2)=_rownum-1;blank(lr,3)=_colnum-1;
	}
}

//Menu:[RBLAST]-[_Fillblank]-[_Maxscwp]
//Find score waypoint for blank area after [_Converge]
//Take 2 values."_BKMATX" Needleman global alignment matrix."wp" waypoint darray
//Update 1 value."wp" waypoint darray:column 0 store row id,column 1 store column id.
void _Maxscwp(darray<int> & _BKMATX,darray<int> & wp)
{
	int i=0,j=0,r,s,d,start,row=_BKMATX.getrnum()-1,col=_BKMATX.getcnum()-1;
	wp.fast_resize(0,2);
	bool ropen=false,dopen=false;
	while(i<row && j<col)
	{	start=_BKMATX(i,j);		//scoer of the current cell
		r=_BKMATX(i,j+1);		//score of the right cell
		d=_BKMATX(i+1,j);		//score of the down cell
		s=_BKMATX(i+1,j+1);		//scoer of the diagonal cell
	//when no gap open occur                   //when right gap open                    //when down gap open     
        if(((r!=start && d!=start) || (s>=r && s>=d))||(ropen==true && (r!=start || s>=r))||(dopen==true && (d!=start || s>=d)))
		{  	ropen=false;
			dopen=false;
			wp.push_row(2);
			wp(wp.getrnum()-1,0)=i;
			wp(wp.getrnum()-1,1)=j;
			i+=1;j+=1;
		}
		          //when no gap open occur                   //when right gap open,r only needs to compare with s
		else if((ropen==false && dopen==false && r>=d) || (ropen==true && r==start && r>s))          //when down gap open,no right gap is allowed
		{	ropen=true;
			j+=1;
		}
		else
		{	dopen=true;
			i+=1;
		}
	}
}

//Menu:[RBLAST]-[_Fillblank]
//Alignment the blank area by using Needleman global alignment technique 
//Take 4 values."_blank" blank area coordiates."_qery" quiry sequence."_bseq" base sequence."alnmatx" standard 0 or 1 waypoint matrix
//Update "alnmatx",which should be formatted into standard 0 or 1 form
void _Fillblank(darray<int> & _blank,sarray<int> & _qidseq,sarray<int> & _tidseq,darray<int> & alnmatx)
{
	int brnum;
	if((brnum=_blank.getrnum())==0)                                             //return if no blank area found
		return;
	else
	{
		darray<int> BKMATX,way;
		sarray<int> subseqa,subseqb;
		for(int b=0;b<brnum;++b)
		{	if(_blank(b,2)==_blank(b,0) && _blank(b,3)==_blank(b,1))	//if start point equal to end point,no alignment is needed 
			{	alnmatx(_blank(b,0),_blank(b,1))=1;
				continue;
			}
			else
			{	subseqa=_qidseq.subseq(_blank(b,0),_blank(b,2));	//fetch out subseq
				subseqb=_tidseq.subseq(_blank(b,1),_blank(b,3));
				NWalign(subseqa,subseqb,BKMATX);
				_Maxscwp(BKMATX,way);
				for(int i=0;i<way.getrnum();++i)
					alnmatx(_blank(b,0)+way(i,0),_blank(b,1)+way(i,1))=1;       //fill in alnmatx
			}
		}
	}
}

//Menu:[RBLAST]-[_Getscore]
//Compute score for standard 0 or 1 waypoint matrix
//Take 6 values."_scmatx" score matrix."_alnmatx" standard 0 or 1 alignment matrix."_rownum"/"_colnum" row/column number of "_alnmatx"."_qery"/"_bseq" quiry/base sequence."_gappnt" gap penalty
//Return score.
template<class type>
type _Getscore(darray<type> & _scmatx,darray<int> & _alnmatx,const int & _rownum,const int & _colnum,sarray<int> & _qidseq,sarray<int> & _tidseq,const int & _gappnt)
{
	type score=0;
	int gap=0,previ=-1,prevj=-1;
	for(int i=previ+1;i<_rownum;++i)		//start from previous point
	{	for(int j=prevj+1;j<_colnum;++j)
		{	if(_alnmatx(i,j)==1)
			{	if(i!=0 && j!=0)	//gaps at the end of alignment were not considered
				   gap=i-previ+j-prevj-2;
				score+=_scmatx(_qidseq[i],_tidseq[j])+(type)gap*_gappnt;     
				previ=i;prevj=j;
				break;
			}
		}
	}
	return score;
}
//Compute score for blank unfilled BLAST alignment matrix
//Take 6 values."_scmatx"scoring matrix."_BLASTmatx" blank unfilled BLAST matrix."_rownum"/"_colnum" row/column number."_qidseq"/"_tidseq" qery and base AA id seq
template<class type>
type _Getscore(darray<type> & _scmatx,const int & _rownum,const int & _colnum,sarray<int> & _qidseq,sarray<int> & _tidseq,const char & _distype,darray<int> & _alnmatx,double & dist,double & var,double & smlt)
{
	int match=0;
	type score=0;int i=0,j=0,prevj=0;
	while(i<_rownum && j<_colnum)                                           
	{	while(i<_rownum && j<_colnum && _alnmatx(i,j)==1)
		{	if(_qidseq[i]==_tidseq[j])
				match+=1;
			score+=_scmatx(_qidseq[i],_tidseq[j]);
			i+=1;j+=1;
			prevj=j;	//always start from the column next to the column where 1 was found.
		}
		j+=1;
		if(j==_colnum)
		{	i+=1;
			j=prevj;	//i always increases by 1 when one row has been scanned
		}
	}
	int len=(_rownum>_colnum ? _colnum : _rownum);
	AprxDist(_distype,match,len,dist,var,smlt);
	return score;
}

template<class type>                                                 
void RBLAST(WDlibstr & _qlib,WDlibstr & _tlib,sarray<int> & _qidseq,sarray<int> & _tidseq,const int & _wlen,vector<darray<type> >& _scmatxvec,const int & _scmatxid,darray<int> & alnmatx,type & score,type & ngscore,double & dist,double & var,double & smlt,const int &_gappnt=-3,const char &_sctype='p',const char & _distype='s')
{
	darray<int> blank;
	BLAST(_qlib,_tlib,_qidseq,_tidseq,_scmatxvec,_scmatxid,_wlen,alnmatx,_sctype);
	_Blankmark(alnmatx,_qidseq.size(),_tidseq.size(),blank);	//find blank area not covered by word matching rule
	_Fillblank(blank,_qidseq,_tidseq,alnmatx);		//fill in blank area by Needleman global alignment rule
	score=_Getscore(_scmatxvec[_scmatxid],alnmatx,_qidseq.size(),_tidseq.size(),_qidseq,_tidseq,_gappnt);
	ngscore=_Getscore(_scmatxvec[_scmatxid],_qidseq.size(),_tidseq.size(),_qidseq,_tidseq,_distype,alnmatx,dist,var,smlt);
}

template<class type>
struct getRBLASTstr
{
	vector<sarray<int> > * seqdb;
	vector<WDlibstr> * WDlibvec;
	vector<darray<type> > * scmatxvec;
	vector<darray<int> > * alnmatxvec;
	darray<int> * ttm;
	sarray<type> * alnscore;
	sarray<type> * alnngscore;
	darray<double> * distmatx;
	darray<double> * varmatx;
	darray<double> * smltmatx;
	int wlen;
	int scmatxid;
	int gappnt;
	char sctype;
	char distype;
	int thrdid;
};

template<class type>
void* _THRBLAST(void * param)
{	
	getRBLASTstr<type>* dtstr=static_cast<getRBLASTstr<type>* >(param);
	int thrdid=dtstr->thrdid;
        int start,end,size,quantum,I,J;		//"start""end"thread task start and end position."quantum" aliquot of task
        size=dtstr->ttm->getcnum();
        quantum=(int)((double)size/CPUNUM);
        if(quantum<1)			//if CPU number is larger than task number,set "quantum" to 1,otherwise it is 0.
                quantum=1;
        start=thrdid*quantum;
        if(thrdid!=CPUNUM-1)
        end=(thrdid+1)*quantum-1;
        else
                end=size-1;
        if(end>=size)
                return NULL;	
	for(int i=start;i<=end;++i)
	{
		I=dtstr->ttm->operator()(0,i);
		J=dtstr->ttm->operator()(1,i);
		RBLAST(dtstr->WDlibvec->operator[](I),dtstr->WDlibvec->operator[](J),dtstr->seqdb->operator[](I),dtstr->seqdb->operator[](J),dtstr->wlen,*(dtstr->scmatxvec),dtstr->scmatxid,dtstr->alnmatxvec->operator[](i),dtstr->alnscore->operator[](i),dtstr->alnngscore->operator[](i),dtstr->distmatx->operator()(I,J),dtstr->varmatx->operator()(I,J),dtstr->smltmatx->operator()(I,J),dtstr->gappnt,dtstr->sctype,dtstr->distype);
	}	
}

template<class type>
void RBLAST(vector<sarray<int> > & _seqdb,const int & _wlen,vector<darray<type> > & _scmatxvec,const int & _scmatxid,vector<darray<int> > & alnmatxvec,sarray<type> & alnscore,sarray<type> & alnngscore,darray<double> & distmatx,darray<double> & varmatx,darray<double> & smltmatx,const int & _gappnt=-3,const char & _sctype='p',const char & _distype='s')
{

	darray<int> ttm=generatethreadid(_seqdb.size());
	alnmatxvec.resize(ttm.getcnum());
	alnscore.resize(ttm.getcnum());
	alnngscore.resize(ttm.getcnum());
	int sqn=_seqdb.size();
	distmatx.fast_resize(sqn,sqn);
	varmatx.fast_resize(sqn,sqn);
	smltmatx.fast_resize(sqn,sqn);
	pthread_t threads[CPUNUM];	
	vector<getRBLASTstr<type> > dtstr(CPUNUM);
	vector<WDlibstr> WDlibvec;
        CRlib(_seqdb,_wlen,WDlibvec);

	for(int i=0;i<CPUNUM;++i)
	{	dtstr[i]={&_seqdb,&WDlibvec,&_scmatxvec,&alnmatxvec,&ttm,&alnscore,&alnngscore,&distmatx,&varmatx,&smltmatx,_wlen,_scmatxid,_gappnt,_sctype,_distype,i};
		pthread_create(&threads[i],NULL,_THRBLAST<type>,(void*)&dtstr[i]);
	}
	for(int i=0;i<CPUNUM;++i)
		pthread_join(threads[i],NULL);
}	
 

	
