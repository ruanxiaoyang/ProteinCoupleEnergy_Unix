void _identical(darray<int> & _combseq,sarray<char> & identical)
{
	int rownum=_combseq.getrnum();
	int colnum=_combseq.getcnum();
	int start,i;
	identical.resize(colnum,' ');
	for(int j=0;j<colnum;++j)
	{	start=_combseq(0,j);
		i=1;
		while(i<rownum && _combseq(i,j)==start)
		{	i+=1;
		}
		if(i==rownum)
			identical[j]='*';
	}
}		

void DspMSA(darray<int> & _combseq,ostream & os=cout)
{
	int rownum=_combseq.getrnum();
	int colnum=_combseq.getcnum();
	int width=(int)(log((double)rownum)/log(10.0))+3; 
	sarray<char> offset(width,' ');
  	
	sarray<char> identical;
	_identical(_combseq,identical);
	os<<offset;LandMarker(colnum,os);os<<endl;
	os<<offset;MileStone(colnum,os);os<<endl;
	os<<offset;os<<identical<<endl;
	for(int i=0;i<rownum;++i)
	{	os.setf(ios::right);
		os.width(width-1);
		os<<i<<" ";
		for(int j=0;j<colnum;++j)
		{	os<<idaa(_combseq(i,j));
		}
		os.setf(ios::right);
		os.width(width);
		os<<i;
		os<<endl;
	}
	os<<offset;os<<identical<<endl;
	os<<offset;MileStone(colnum,os);os<<endl;
	os<<offset;LandMarker(colnum,os);os<<endl;

}

void DspMSAseqname(vector<inf> & _infvec,ostream & os=cout)
{
	int rownum=_infvec.size();
	int width=(int)(log((double)rownum)/log(10.0))+3; 
	sarray<char> offset(width,' ');
  	
	for(int i=0;i<rownum;++i)
	{	os.setf(ios::right);
		os.width(width-1);
		os<<i<<" ";
		os<<_infvec[i]._gi;
		os<<endl;
	}
}


//menu:[MSA]-[_Alignseqs]-[_SeqMge]
//merge two sarray<int> sequence according to "_alnmatx" to give two properly aligned seqs stored in darray<int>
//take 3 values."_qery"/"_bseq" query and base sequence."_alnmatx" alignment matrix(row for query,column for target).
//output 1 value."combseq" two aligned seqs stored in darray<int>,first row store target,second row store query,gap was represented by '-' char
void _SeqMge(sarray<int> & _qidseq,sarray<int> & _tidseq,darray<int> & _alnmatx,darray<int> & combseq)
{	
	int rownum=_qidseq.size(),colnum=_tidseq.size();
	int i=0,j=0,previ=0,prevj=0,posi=0,posj=0,gapi,gapj,size;
	size=Length(_alnmatx,rownum,colnum);	//length of aligned seq,use it to initiate "combseq"
	combseq.fast_resize(2,size,-4);
	while(i<rownum && j<colnum)
	{	while(i<rownum && j<colnum && _alnmatx(i,j)==1)
		{
			gapi=i-previ;
			gapj=j-prevj;
			posi+=gapj;
			posj+=gapi;
			for(int I=0;I<=gapi;++I)
				combseq(1,posi+I)=_qidseq[previ+I];
			for(int J=0;J<=gapj;++J)
				combseq(0,posj+J)=_tidseq[prevj+J];
			i+=1;
			j+=1;
			posi+=(gapi+1);
			posj+=(gapj+1);
			previ=i;
			prevj=j;
		
		}
		j+=1;
		if(j==colnum)
		{	i+=1;
			j=prevj;
		}
	}
	for(int i=previ;i<rownum;++i)
		combseq(1,posi+(i-previ))=_qidseq[i];
	for(int j=prevj;j<colnum;++j)
		combseq(0,posj+(j-prevj))=_tidseq[j];
}

void _SeqMge(darray<int> & _qidseq,darray<int> & _tidseq,darray<int> & _alnmatx,darray<int> & combseq)
{	
	int rownum=_qidseq.getcnum(),colnum=_tidseq.getcnum();
	int qn=_qidseq.getrnum(),tn=_tidseq.getrnum();

	int i=0,j=0,previ=0,prevj=0,posi=0,posj=0,gapi,gapj,size;
	size=Length(_alnmatx,rownum,colnum);	//length of aligned seq,use it to initiate "combseq"
	combseq.fast_resize(qn+tn,size,-4);
	while(i<rownum && j<colnum)
	{	while(i<rownum && j<colnum && _alnmatx(i,j)==1)
		{
			gapi=i-previ;
			gapj=j-prevj;
			posi+=gapj;
			posj+=gapi;
			for(int q=0;q<qn;++q)
			{	for(int I=0;I<=gapi;++I)
					combseq(tn+q,posi+I)=_qidseq(q,previ+I);
			}
			for(int t=0;t<tn;++t)
			{	for(int J=0;J<=gapj;++J)
					combseq(t,posj+J)=_tidseq(t,prevj+J);
			}
			i+=1;
			j+=1;
			posi+=(gapi+1);
			posj+=(gapj+1);
			previ=i;
			prevj=j;
		
		}
		j+=1;
		if(j==colnum)
		{	i+=1;
			j=prevj;
		}
	}
	for(int q=0;q<qn;++q)
	{	for(int i=previ;i<rownum;++i)
			combseq(tn+q,posi+(i-previ))=_qidseq(q,i);
	}
	for(int t=0;t<tn;++t)
	{	for(int j=prevj;j<colnum;++j)
			combseq(t,posj+(j-prevj))=_tidseq(t,j);
	}
}




//Menu:[MSA]-[_Alignseqs]-[_SeqMge]
//Merge two darray<char> seqs into a combined one according to Needleman global rule
//Take 3 values."_marr1"/"_marr2" multiple(or sinlge) seqs store as darray<char>."_gappnt" gap penalty
//Output 1 value."combseq" final combination of former seqs 
void _SeqMge(darray<int> & _marrid1,darray<int> & _marrid2,darray<int> & combseq,const int & _gappnt)
{
	int rownum=_marrid1.getcnum(0),colnum=_marrid2.getcnum(0);
	int qn=_marrid1.getrnum(),tn=_marrid2.getrnum();
	int id1,id2,gapi,gapj,size;
	double tmpscore,s,max;
	darray<int> alnmatx(rownum,colnum,0);
	darray<double> BKMATX(rownum+1,colnum+1,0.0);
	for(int j=colnum-1;j>=0;--j)                                                   //calculate BKMATX
	{	for(int i=rownum-1;i>=0;--i)
	    	{	tmpscore=0;
			for(int k=0;k<tn;++k)
			{	for(int m=0;m<qn;++m)
				{	id1=_marrid1(m,i);
					id2=_marrid2(k,j);
					if(id1>=0 && id2>=0)
						tmpscore+=Gonnet(id1,id2);
					else if(id1<0 && id2<0)
						tmpscore+=2;
				}
			}
			s=tmpscore+BKMATX(i+1,j+1);
			max=(BKMATX(i+1,j)>BKMATX(i,j+1)?BKMATX(i+1,j):BKMATX(i,j+1));
			BKMATX(i,j)=(s>max?s:max);
		}
	}
	int I=0,J=0,rgap=0,dgap=0;
	double rsc,dsc,ssc,r,d,start;
	while(I<rownum && J<colnum)		//calculate standard 0 or 1 waypoint matrix
	{	start=BKMATX(I,J);
		rgap+=1;
		dgap+=1;
		r=BKMATX(I,J+1);
		d=BKMATX(I+1,J);
		s=BKMATX(I+1,J+1);
		rsc=r+rgap*_gappnt;
		dsc=d+dgap*_gappnt;
		ssc=s;
		if(((r!=start && d!=start) || (s>=r && s>=d)||(ssc>=rsc && ssc>=dsc))|| (rgap>1 && (r!=start || s>=r || ssc>=rsc)) || (dgap>1 && (d!=start || s>=d|| ssc>=dsc)))
		{	rgap=0;dgap=0;
			alnmatx(I,J)=1;
			I+=1;J+=1;
		}
		else if((rgap==1 && dgap==1 && (r>d ||rsc>=dsc)) || (rgap>1 && r>s && rsc>ssc))
		{	dgap=0;
			J+=1;
		}
		else
		{	rgap=0;
			I+=1;
		}
	}
	size=Length(alnmatx,rownum,colnum);	//calculate the volume of combined sequence
	combseq.fast_resize(qn+tn,size,-4);                                        
	_SeqMge(_marrid1,_marrid2,alnmatx,combseq);
	return;
}

//overloaded functon of _SeqMge(darray<char> & _marr1,darray<char> & _marr2,darray<char> & combseq) to take sarray<char> type in the first argument
void _SeqMge(sarray<int> & _arrid1,darray<int> & _marrid2,darray<int> & combseq,const int & _gappnt)
{
	darray<int> marrid1;
	marrid1.push_row(_arrid1);
	_SeqMge(marrid1,_marrid2,combseq,_gappnt);
}
//overloaded functon of _SeqMge(darray<char> & _marr1,darray<char> & _marr2,darray<char> & combseq) to take sarray<char> type in the second argument
void _SeqMge(darray<int> & _marrid1,sarray<int> & _arrid2,darray<int> & combseq,const int & _gappnt)
{
	darray<int> marrid2;
	marrid2.push_row(_arrid2);
	_SeqMge(_marrid1,marrid2,combseq,_gappnt);
}


//Menu:[MSA]-[_MSAsc]
//Calculate score for MSA.The scoring matrix is Gonnet
//Take 2 value."_combseq"."_gappnt" gap penalty
//Return score;
double _MSAsc(darray<int> & _combseq,const int & _gappnt=-3)
{
        int I,J,row=_combseq.getrnum();
        double sc=0.0;
        for(int j=0;j<_combseq.getcnum();++j)
        {	for(int i=0;i<row;++i)
                {	for(int ii=i+1;ii<row;++ii)
                        {	I=_combseq(i,j);
                                J=_combseq(ii,j);
                                if(I>=0 && J>=0)
					sc+=Gonnet(I,J);
                                else if ((I<0 && J>=0)||(I>=0 && J<0))
                                        sc+=(double)_gappnt;
                        }
                }
        }
        sc/=((double)(row*(row-1))/2.0);
        return sc;
}

//Menu:[MSA]-[_OmitBlank] [MSA]-[_Alignseqs]-[_OmitBlank]
//Omit postions that was purely composed of gaps
//Take 2 values."_combseq" int or char seqs."_posicnt" position count
//Output 1 value."identical" sarray storing positions with identical AA
//Update 1 value."_combseq" pure gap will be omitted
template <class type>
void _OmitBlank(darray<type> & _combseq,darray<int> & _posicnt)
{
	int colnum=0,max,iter=0,index;
	sarray<int> nonemptyposi(_combseq.getcnum());
	darray<type> tmpcombseq;

	for(int j=0;j<_posicnt.getcnum();++j)
	{	if((max=_posicnt.columnmax(j,0,index))!=0)
		{	colnum+=1;
			nonemptyposi[iter]=j;
			iter+=1;
		}
	}
	if(colnum==_combseq.getcnum())                                          //if no pure gap found,no update is needed
		return;
	else
	{	tmpcombseq.fast_resize(_combseq.getrnum(),colnum);
		for(int j=0;j<colnum;++j)
		{	for(int i=0;i<_combseq.getrnum();++i)
				tmpcombseq(i,j)=_combseq(i,nonemptyposi[j]);
		}
		_combseq=tmpcombseq;
	}
}

//Menu:[MSA]-[_Iterate] [MSA]-[_Alignseqs]-[_Iterate]
//This module refines the MSA result to correct obvious errors
//Take 3 values."_combseq" MSA resulted."_ancthld" anchoring point threshold."_amp" gap amplification factor
//update 2 values."posicnt"position AA count of the rectified matrix. "_combseq"
void _Iterate(darray<int> & _combseq,const int & _ancthld,const int & _amp,darray<int> & posicnt)
{
	
	int row=_combseq.getrnum(),col=_combseq.getcnum(),_index_,J,origsc,origdist,markJ,aaid;
	posicnt.fast_resize(20,col,0);
	sarray<int> anchor,ancdist(col);
	for(int j=0;j<col;++j)
	{	for(int i=0;i<row;++i)
		{	if(_combseq(i,j)>=0)
				posicnt(_combseq(i,j),j)+=1;
		}
		if(posicnt.columnmax(j,0,_index_)>=(int)((double)(row*_ancthld)/100.0))
			anchor.pushback(j);
	}
	if(anchor.size()==0)
	{	anchor.resize(col);
		for(int j=0;j<col;++j)
			anchor[j]=j;
	}

	int a=0,pnt=_amp*((int)(row/10.0)==0?1:(int)(row/10.0));//calculate the min distance of each position to anchoring point.
	for(int j=0;j<col;++j)					//multiply the distance with amplification factor
	{	if(a!=anchor.size() && j<anchor[a])
		{	if(a==0)
			   ancdist[j]=pnt*(anchor[a]-j);
			else if(a>0 && a<=anchor.size()-1)
				ancdist[j]=((j-anchor[a-1])>(anchor[a]-j)?pnt*(anchor[a]-j):pnt*(j-anchor[a-1]));
		}
		else if(a!=anchor.size() && j==anchor[a])
		{	ancdist[j]=0;
			a+=1;
		}
		else if(a==anchor.size())
			ancdist[j]=pnt*(j-anchor[a-1]);
	}
	bool shifted=false;	//"shifted" controls whether a new round of iteration is needed
	for(int i=0;i<row;++i)
	{	for(int j=0;j<col;++j)
		{	if(_combseq(i,j)<0)
				continue;
			else
			{
				aaid=_combseq(i,j);
				_combseq(i,j)=-4;
				origdist=ancdist[j];
				origsc=posicnt(aaid,j)-origdist;//score of a position is the count number of the AA minus gap penalty  
				J=j;
				while(J>=0 && _combseq(i,J)<0)	//shift to the left gap
					J-=1;
				J+=1;
				while(J<col && _combseq(i,J)<0)	//search the blank space for a better position
				{	if((posicnt(aaid,J)-ancdist[J])>origsc || ((posicnt(aaid,J)-ancdist[J])==origsc && ancdist[J]<=origdist))
					{     //J has a better score  //J has same score but lower distance                               
						markJ=J;
						origdist=ancdist[J];
						origsc=posicnt(aaid,J)-origdist;
					}
					J+=1;
				}
				_combseq(i,markJ)=aaid;
				if(markJ!=j)
				{
					shifted=true;
					posicnt(aaid,markJ)+=1;	//update the position count
					posicnt(aaid,j)-=1;
				}
				j=J-1;
			}
		}
	}
	if(shifted==true)	//run another round of iteration until convergence
		_Iterate(_combseq,_ancthld,_amp,posicnt);
	return;
}

//Menu:[MSA]-[_Alignseqs]
//This module align multiple seqs.It is the main component of [MSA] 
//"_nbtree" structure hold evolutionary tree(see EvoTree for detail)
//"_mwit"midway iteration switcher."_mancthld"midway anchoring point threshold."_mamp"midway gap amplification factor
//Other parameters see [MSA]
void _Alignseqs(vector<sarray<int> > & _seqdb,vector<darray<int> > & _alnmatxvec,nbstr & _nbtree,darray<int> & seqsmatx,const int & _gappnt=-3,const char & _mwit='y',const int & _mancthld=95,const int & _mamp=3)
{
	int dist,qseqid,tseqid,id,size;
	size=_seqdb.size();
	vector<darray<int> > combseqvec(_nbtree._nbsize);
	darray<int> combseq,alnmatx,posicnt;
	cout<<"Aligning sequence...";
	for(int i=0;i<_nbtree._nbsize;++i)
	{	PctMarker(i,_nbtree._nbsize-1);
		if(_nbtree._lkbk(i,0)<0 && _nbtree._lkbk(i,1)<0)//when two seqs are both single seq,use RBLAST alignment
		{
			qseqid=_nbtree._marr1(i,0);
			tseqid=_nbtree._marr2(i,0);
			id=(int)((double)qseqid*((double)size-(double)(qseqid+1)/2.0))+tseqid-qseqid;
			alnmatx=_alnmatxvec[id-1];
			_SeqMge(_seqdb[qseqid],_seqdb[tseqid],alnmatx,combseq);	//create aligned seqs
		}
		else if(_nbtree._lkbk(i,0)>=0 && _nbtree._lkbk(i,1)<0)
			_SeqMge(combseqvec[_nbtree._lkbk(i,0)],_seqdb[_nbtree._marr2(i,0)],combseq,_gappnt);
		else if(_nbtree._lkbk(i,0)<0 && _nbtree._lkbk(i,1)>=0)
			_SeqMge(_seqdb[_nbtree._marr1(i,0)],combseqvec[_nbtree._lkbk(i,1)],combseq,_gappnt);
		else if(_nbtree._lkbk(i,0)>=0 && _nbtree._lkbk(i,1)>=0)
			_SeqMge(combseqvec[_nbtree._lkbk(i,0)],combseqvec[_nbtree._lkbk(i,1)],combseq,_gappnt);
		if(combseq.getrnum()>=3 && (_mwit=='Y'||_mwit=='y'))    //do midway iteration if it is enabled,
		{
			_Iterate(combseq,_mancthld,_mamp,posicnt);
            		_OmitBlank(combseq,posicnt);
		}
		combseqvec[i]=combseq;
		combseq.clear();
	}
	seqsmatx=combseqvec[_nbtree._nbsize-1];
	cout<<endl;
}

//reorder combseq to the order appeared in original file
//take 2 values."_combseq" the darray<int> holding all aligned seqs."_evotree" a nbstr struct object
//update "_combseq" to new order.
void reorder(darray<int> & _combseq,nbstr & _evotree)
{
	sarray<int> oa=_evotree._marr1[_evotree._marr1.getrnum()-1];
	sarray<int> ob=_evotree._marr2[_evotree._marr2.getrnum()-1];
	sarray<int> order=oa.c(ob);
	darray<int> tmp=_combseq;
	int id,iter=0;
	for(int i=order.size()-1;i>=0;--i)
	{
		id=order[i];
		for(int j=0;j<tmp.getcnum();++j)
			_combseq(id,j)=tmp(iter,j);
		iter+=1;
	}	
}

void reorderseqname(proformat & _prodb,nbstr & _evotree)
{
	sarray<int> oa=_evotree._marr1[_evotree._marr1.getrnum()-1];
        sarray<int> ob=_evotree._marr2[_evotree._marr2.getrnum()-1];
        sarray<int> order=oa.c(ob);
	vector<inf> seqinf(_prodb._seqnum);
	int iter=0;
	for(int i=order.size()-1;i>=0;--i)
	{	seqinf[iter]._gi=_prodb._infvec[i]._gi;
		iter+=1;
	}
	_prodb._infvec=seqinf;
}
		
	


