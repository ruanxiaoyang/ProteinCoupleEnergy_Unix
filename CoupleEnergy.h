void _getpx2pmsa(darray<int> & _combseq,const int & _sta,const int & _end,sarray<double> & _comb,darray<double> & px2pmsa,const char & _scaleup='n')
{
	int row=_combseq.getrnum(),origrow=_comb.size()-1,n=origrow;
	int col=_end-_sta+1,aa;
	for(int i=0;i<row;++i)
	{	for(int j=_sta;j<=_end;++j)
		{	if((aa=_combseq(i,j))>=0)
				px2pmsa(aa,j)+=1.0;
		}
	}
	if(_scaleup=='y')
	{	for(int i=0;i<px2pmsa.getrnum();++i)
		{	for(int j=_sta;j<=_end;++j)
			{	px2pmsa(i,j)=px2pmsa(i,j)*(origrow/(double)row)+0.5;
			}
		}
	}
	double k,p;
	for(int a=0;a<20;++a)
	{	for(int j=_sta;j<=_end;++j)
		{	k=px2pmsa(a,j);
			p=aafreq[a];
			px2pmsa(a,j)=_comb[(int)k]+(k-1)*log(p)+(n-k)*log(1-p);
		}
	}		
}

struct getpx2pmsastr
{
	darray<int> * combseq;
	sarray<double> * comb;
	darray<double> * px2pmsa;
	char scaleup;
	int thrdid;
};

void* _THpx2pmsastr(void* param)
{
	getpx2pmsastr* dtstr=static_cast<getpx2pmsastr*>(param);
	int thrdid=dtstr->thrdid;
        int start,end,size,quantum;             //"start""end"thread task start and end position."quantum" aliquot of task
        size=dtstr->combseq->getcnum();
        quantum=(int)((double)size/CPUNUM);
        if(quantum<1)                   //if CPU number is larger than task number,set "quantum" to 1,otherwise it is 0.
                quantum=1;
        start=thrdid*quantum;
        if(thrdid!=CPUNUM-1)
        end=(thrdid+1)*quantum-1;
        else
                end=size-1;
        if(end>=size || start>=size)
                return NULL;
	_getpx2pmsa(*dtstr->combseq,start,end,*dtstr->comb,*dtstr->px2pmsa,dtstr->scaleup);	
}

void _getcombination(const int & _n,sarray<double> & comb)
{
	comb.resize(_n+1);
	int mid=(_n+1)/2;
	double a=0,b=0;
	comb[0]=comb[_n]=0;
	for(int k=1;k<=mid;++k)
	{	a+=log(_n-(k-1));
		b+=log(k);
		comb[k]=comb[_n-k]=a-b;
	}
}

void _getpx2pmsa(darray<int> & _combseq,sarray<double> & _comb,darray<double> & px2pmsa,const char & _scaleup='n')
{
	px2pmsa.fast_resize(20,_combseq.getcnum(),0);
	vector<getpx2pmsastr> dtstr(CPUNUM);
	pthread_t threads[CPUNUM];
	int code=0;
	for(int i=0;i<CPUNUM;++i)
	{	dtstr[i]={&_combseq,&_comb,&px2pmsa,_scaleup,i};
		code=pthread_create(&threads[i],NULL,_THpx2pmsastr,(void*)(&dtstr[i]));
		if(code!=0)
		{	cout<<"Thread number exceeded system limit. You are running "<<CPUNUM<<" CPUs. Change to a smaller number with -c option"<<endl;
			exit(0);
		}
	}
	for(int i=0;i<CPUNUM;++i)
		pthread_join(threads[i],NULL);
}
	
void _getposaacnt(darray<int> & _combseq,darray<int> & posaacnt)
{	
	posaacnt.fast_resize(20,_combseq.getcnum(),0);
	int id=0;
	for(int i=0;i<_combseq.getrnum();++i)
	{	for(int j=0;j<_combseq.getcnum();++j)
		{	if((id=_combseq(i,j))>=0)
				posaacnt(id,j)+=1;
		}
	}
}

void _getcped(darray<double> & _pp1,darray<double> & _pp2,sarray<double> & cped)
{
	cped.resize(_pp1.getcnum(),0);
	darray<double> diff;
	diff=_pp1-_pp2;
	for(int i=0;i<diff.getrnum();++i)
	{	for(int j=0;j<diff.getcnum();++j)
			cped[j]+=pow(diff(i,j),2.0);
	}
	for(int j=0;j<cped.size();++j)
		cped[j]=pow(cped[j],0.5);
}

void _subcombseq(darray<int> & _combseq,const int & _j,const int & _aaid,const int & _num,darray<int> & subcombseq)
{
	subcombseq.fast_resize(_num,_combseq.getcnum());
	int iter=0;
	for(int i=0;i<_combseq.getrnum();++i)
	{	if(_combseq(i,_j)==_aaid)
		{	for(int j=0;j<_combseq.getcnum();++j)
				subcombseq(iter,j)=_combseq(i,j);
			iter+=1;
		}
	}
}

struct cpedstr
{
	int pos;
	sarray<int> aaidarr;
	vector<sarray<double> > poscpedvec;
};

sarray<int> __lt_trhd(darray<int> & _posaacnt,const int & _j,const int & _trhd,const int & _rownum)
{	
	sarray<int> res;
	for(int i=0;i<_posaacnt.getrnum();++i)
	{	if(_posaacnt(i,_j)>=_trhd && _rownum-_posaacnt(i,_j)>=_trhd)
			res.pushback(i);
	}
	return res;
}

void CPED(darray<int> & _combseq,sarray<double> & _comb,darray<double> & _ppt,const int & _sta,const int & _end,darray<int> & _posaacnt,const int & _trhd,vector<cpedstr> & cpedresvec)
{
	darray<int> subcombseq;
	darray<double> pps;
	for(int j=_sta;j<=_end;++j)
	{	cpedresvec[j].pos=-1;
		cpedresvec[j].aaidarr=__lt_trhd(_posaacnt,j,_trhd,_combseq.getrnum());
		if(cpedresvec[j].aaidarr.size()>0)
		{	cpedresvec[j].pos=j;
			cpedresvec[j].poscpedvec.resize(cpedresvec[j].aaidarr.size());	
			for(int i=0;i<cpedresvec[j].aaidarr.size();++i)
			{	_subcombseq(_combseq,j,cpedresvec[j].aaidarr[i],_posaacnt(cpedresvec[j].aaidarr[i],j),subcombseq);
				_getpx2pmsa(subcombseq,_comb,pps,'y');
				_getcped(_ppt,pps,cpedresvec[j].poscpedvec[i]);
			}
		}
	}
	
}

struct getCPEDstr
{
	darray<int> * combseq;
	sarray<double> * comb;
	darray<int> * posaacnt;
	darray<double> * ppt;	
	vector<cpedstr> * cpedresvec;
	int trhd;
	int thrdid;
};

void* _THCPED(void* param)
{
	getCPEDstr* dtstr=static_cast<getCPEDstr*>(param);
	int thrdid=dtstr->thrdid;
        int start,end,size,quantum;             //"start""end"thread task start and end position."quantum" aliquot of task
        size=dtstr->combseq->getcnum();
        quantum=(int)((double)size/CPUNUM);
        if(quantum<1)                   //if CPU number is larger than task number,set "quantum" to 1,otherwise it is 0.
                quantum=1;
        start=thrdid*quantum;
        if(thrdid!=CPUNUM-1)
        end=(thrdid+1)*quantum-1;
        else
                end=size-1;
        if(end>=size || start>=size)
                return NULL;
	CPED(*dtstr->combseq,*dtstr->comb,*dtstr->ppt,start,end,*dtstr->posaacnt,dtstr->trhd,*dtstr->cpedresvec);
}	

//Coupling Energe Difference
void CPED(darray<int> & _combseq,const int & _trhd,vector<cpedstr> & cpedresvec)
{
	cpedresvec.resize(_combseq.getcnum());
	darray<int> posaacnt;
	_getposaacnt(_combseq,posaacnt);

	sarray<double> comb;		
	_getcombination(_combseq.getrnum(),comb);

	darray<double> ppt;
	_getpx2pmsa(_combseq,comb,ppt);	

	vector<getCPEDstr> dtstr(CPUNUM);
	pthread_t threads[CPUNUM];
	int code;
	for(int i=0;i<CPUNUM;++i)
	{	dtstr[i]={&_combseq,&comb,&posaacnt,&ppt,&cpedresvec,_trhd,i};
		code=pthread_create(&threads[i],NULL,_THCPED,(void*)(&dtstr[i]));
		if(code!=0)
		{	cout<<"Thread number exceeded system limit. You are running "<<CPUNUM<<" CPUs. Change to a smaller number with -c option"<<endl;
			exit(0);
		}
	}
	for(int i=0;i<CPUNUM;++i)
		code=pthread_join(threads[i],NULL);
}







#define ABS(X) ((X) > 0 ? (X) : (-(X)))


void _getsed(darray<int> & _tt,darray<int> & _sub,const int & _ttrow,const int & _subrow,const int & _j,sarray<double> & sed)
{
	sed.resize(_tt.getcnum(),0);
	for(int i=0;i<_tt.getrnum();++i)
	{	for(int j=0;j<_tt.getcnum();++j)
		{	if(j!=_j)
				sed[j]+=ABS(_tt(i,j)/(double)_ttrow-_sub(i,j)/(double)_subrow);
		}
	}
}

void SED(darray<int> & _combseq,const int & _sta,const int & _end,darray<int> & _posaacnt,const int & _thrd,vector<cpedstr> & cpedresvec)
{
	darray<int> subcombseq,subposaacnt;
	for(int j=_sta;j<=_end;++j)
	{	cpedresvec[j].pos=-1;
		cpedresvec[j].aaidarr=__lt_trhd(_posaacnt,j,_thrd,_combseq.getrnum());
		if(cpedresvec[j].aaidarr.size()>0)
		{	cpedresvec[j].pos=j;
			cpedresvec[j].poscpedvec.resize(cpedresvec[j].aaidarr.size());	
			for(int i=0;i<cpedresvec[j].aaidarr.size();++i)
			{	_subcombseq(_combseq,j,cpedresvec[j].aaidarr[i],_posaacnt(cpedresvec[j].aaidarr[i],j),subcombseq);
				_getposaacnt(subcombseq,subposaacnt);		
				_getsed(_posaacnt,subposaacnt,_combseq.getrnum(),subcombseq.getrnum(),j,cpedresvec[j].poscpedvec[i]);
			}
		}
	}
	
}


struct getSEDstr
{
	darray<int> * combseq;
	darray<int> * posaacnt;
	vector<cpedstr> * cpedresvec;
	int trhd;
	int thrdid;
};


void* _THSED(void* param)
{
	getSEDstr* dtstr=static_cast<getSEDstr*>(param);
	int thrdid=dtstr->thrdid;
        int start,end,size,quantum;             //"start""end"thread task start and end position."quantum" aliquot of task
        size=dtstr->combseq->getcnum();
        quantum=(int)((double)size/CPUNUM);
        if(quantum<1)                   //if CPU number is larger than task number,set "quantum" to 1,otherwise it is 0.
                quantum=1;
        start=thrdid*quantum;
        if(thrdid!=CPUNUM-1)
        end=(thrdid+1)*quantum-1;
        else
                end=size-1;
        if(end>=size)
                return NULL;
	SED(*dtstr->combseq,start,end,*dtstr->posaacnt,dtstr->trhd,*dtstr->cpedresvec);
}	

//Simple Estimated Difference
void SED(darray<int> & _combseq,const int & _trhd,vector<cpedstr> & cpedresvec)
{
	cpedresvec.resize(_combseq.getcnum());
	darray<int> posaacnt;
	_getposaacnt(_combseq,posaacnt);		
	vector<getSEDstr> dtstr(CPUNUM);
	pthread_t threads[CPUNUM];
	int code;
	for(int i=0;i<CPUNUM;++i)
	{	dtstr[i]={&_combseq,&posaacnt,&cpedresvec,_trhd,i};
		code=pthread_create(&threads[i],NULL,_THSED,(void*)(&dtstr[i]));
		if(code!=0)
		{	cout<<"Thread number exceeded system limit. You are running "<<CPUNUM<<" CPUs. Change to a smaller number with -c option"<<endl;
			exit(0);
		}
	}
	for(int i=0;i<CPUNUM;++i)
		pthread_join(threads[i],NULL);
}
/****************************/


//calculate coupling energe difference distribution
int _getEDdist(vector<cpedstr> & _cpedresvec,darray<int> & _combseq,sarray<double> & qtl)
{
	int num=0,iter=0;
	for(int i=0;i<_cpedresvec.size();++i)
		num+=_cpedresvec[i].aaidarr.size()*_combseq.getcnum();
	if(num==0)
		return num;
	sarray<double> distribution(num);
	for(int i=0;i<_cpedresvec.size();++i)
	{	if(!_cpedresvec[i].aaidarr.empty())
		{	for(int a=0;a<_cpedresvec[i].poscpedvec.size();++a)
			{	for(int j=0;j<_cpedresvec[i].poscpedvec[a].size();++j)
				{	distribution[iter]=_cpedresvec[i].poscpedvec[a][j];
					iter+=1;
				}
			}
		}
	}
	
	QuickSort(distribution,0,distribution.size()-1);
	iter=0;				//remove all zero
	while(iter<distribution.size() && distribution[iter]==0)
	{	iter+=1;
	}
	distribution=distribution.subseq(iter,distribution.size()-1);
	if(distribution.size()==0)
		return 0;
	int size;
	size=distribution.size()-1;
	qtl.resize(10);
	double qtlevel=0.9;
	for(int i=0;i<10;++i)
	{
		qtl[i]=distribution[(int)((double)size*qtlevel)];
		qtlevel+=0.01;
	}
	return distribution.size();
}

void DspED(vector<cpedstr> & _cpedresvec,darray<int> & _combseq,const int & _mask,ostream & os=cout)
{
	double cped;
	int colnum=_combseq.getcnum();
	int rownum=_combseq.getrnum();
	int width=(int)(log((double)rownum)/log(10.0))+3; 
	sarray<char> offset(width,' ');
 
	int iter=0;
	sarray<double> qtl;

	int flag=0;
	flag=_getEDdist(_cpedresvec,_combseq,qtl);
	if(flag==0)
	{	os<<"No coupling sites found"<<endl;
		return;
	}

	os<<"Quantile 90-99"<<endl<<qtl<<endl<<endl;

        for(int j=0;j<_cpedresvec.size();++j)
	{	if(!_cpedresvec[j].aaidarr.empty())
		{	os<<"Pos "<<j<<endl;
			for(int i=0;i<_cpedresvec[j].aaidarr.size();++i)
			{       os<<"Mutation to "<<idaa(_cpedresvec[j].aaidarr[i])<<endl;
				os<<offset;LandMarker(colnum,os);os<<endl;
				os<<offset;MileStone(colnum,os);os<<endl;
				os<<offset;
				for(int jj=0;jj<_cpedresvec[j].poscpedvec[i].size();++jj)
				{	cped=_cpedresvec[j].poscpedvec[i][jj];
					if(jj==j)
						os<<idaa(_cpedresvec[j].aaidarr[i]);
					else
					{	iter=0;
						while(iter<qtl.size() && cped>=qtl[iter])
						{	iter+=1;
						}
						iter-=1;
						if(iter<_mask)
							os<<"-";
						else
							os<<iter;
					}					
				}
				os<<endl;	
			}
			os<<endl;
		}
        }
}

	
