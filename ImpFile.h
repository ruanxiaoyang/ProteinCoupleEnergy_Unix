struct inf                                                                  //Structure for all format.
{
	sarray<char> _gi;
	sarray<char> _ginum;
	sarray<char> _format;
	sarray<char> _accession;
	sarray<char> _locus;
};
//////////////////////////////////////////////////////////////////////////////standard fasta format
struct proformat                                                          //structure hold all information of sequence file
{
	vector<sarray<int> > _multiseqid;
	vector<inf> _infvec;
	int _seqnum;
	char _format;
};

int Fasta(char * _path,proformat & seqsandinf)
{
	clock_t start,end;
	start=clock();
	char aa;
	seqsandinf._format='F';
	inf proinf;
	seqsandinf._seqnum=0;
	sarray<int> aaseqid;
	ifstream input(_path);
	if(!input)
	{
		cout<<"Didnt find "<<_path<<endl;
		return 0;                                                             //return 0 if the profile was not found
	}
	while(input.get(aa) && aa=='>')
	{
        	while(input.get(aa))
		{
			if(aa!='|')
			proinf._gi.pushback(aa);
			else
				break;
		}
		while(input.get(aa))
		{
			if(aa!='|')
			proinf._ginum.pushback(aa);
			else
				break;
		}
		while(input.get(aa))
		{
			if(aa!='|')
			proinf._format.pushback(aa);
			else
				break;
		}
		while(input.get(aa))
		{
			if(aa!='|')
			proinf._accession.pushback(aa);
			else
				break;
		}
		while(input.get(aa))
		{
			if(aa!='\n')
			proinf._locus.pushback(aa);
			else
				break;
		}
		seqsandinf._infvec.push_back(proinf);
		proinf._gi.clear();proinf._ginum.clear();proinf._format.clear();proinf._accession.clear();proinf._locus.clear();
		while(input.get(aa))
		{
			if(aa=='>')
				break;
			if((int)aa>=65 && (int)aa<=122)
				aaseqid.pushback(aaid(aa));
		}
		seqsandinf._multiseqid.push_back(aaseqid);
		seqsandinf._seqnum+=1;
		aaseqid.clear();
		input.unget();
	}
	end=clock();
	ofstream outlog("./log",ios_base::app);
	outlog<<getdate()<<" Process time "<<double(end-start)/1000000<<" sec"<<endl;
	outlog<<seqsandinf._multiseqid.size()<<" sequences successfully imported"<<endl<<endl; 
	input.close();
	return 1;
}

/////////////////////////////////////////////////////////////////////////////custom format,same structure but different name

int Custom(const char *_path,proformat & seqsandinf)
{
	clock_t start,end;
	start=clock();
	seqsandinf._format='C';
	char aa;
	seqsandinf._seqnum=0;
	sarray<int> aaseqid;
	inf proinf;
	ifstream input(_path);
	if(!input)
	{
		cout<<_path<<" not found!"<<endl;
		return 0;                                                             //return 0 if the profile was not found
	}    
	
	while(input.get(aa) && aa=='>')
	{
		while(input.get(aa))
		{
			if(aa!='\n')
				proinf._gi.pushback(aa);
			else
				break;
		}
		seqsandinf._infvec.push_back(proinf);
		proinf._gi.clear();
		while(input.get(aa))
		{
			if(aa=='>')
				break;
			if((int)aa>=65 && (int)aa<=122)
				aaseqid.pushback(aaid(aa));
		}
		seqsandinf._multiseqid.push_back(aaseqid);
		seqsandinf._seqnum+=1;
		aaseqid.clear();
		input.unget();
	}
	end=clock();
	ofstream outlog("./log",ios_base::app);
    	outlog<<getdate()<<" Process time "<<double(end-start)/1000000<<" sec"<<endl;
	outlog<<seqsandinf._multiseqid.size()<<" sequences successfully imported"<<endl<<endl; 
	input.close();
	return 1 ;
}

int Custom(const char *_path,proformat & seqsandinf,sarray<int> & range)
{
	clock_t start,end;
	start=clock();
	seqsandinf._format='C';
	char aa;
	seqsandinf._seqnum=0;
	sarray<int> aaseqid;
	inf proinf;
	ifstream input(_path);
	if(!input)
	{
		cout<<_path<<" not found!"<<endl;
		return 0;                                                             //return 0 if the profile was not found
	}    
	int fit=1,rit=0;
	while(input.get(aa))
	{
	   if(aa=='>')
	   {
		if(fit!=range[rit])
		{
			fit+=1;
			continue;
		}
		else 
		{
			fit+=1;
			rit+=1;
		}
		while(input.get(aa))
		{
			if(aa!='\n')
				proinf._gi.pushback(aa);
			else
				break;
		}
		seqsandinf._infvec.push_back(proinf);
		proinf._gi.clear();
		while(input.get(aa))
		{
			if(aa=='>')
				break;
			if((int)aa>=65 && (int)aa<=122)
				aaseqid.pushback(aaid(aa));
		}
		seqsandinf._multiseqid.push_back(aaseqid);
		seqsandinf._seqnum+=1;
		aaseqid.clear();
		input.unget();
		if(rit==range._size)
			break;
		}
	}
	end=clock();
	ofstream outlog("./log",ios_base::app);
    	outlog<<getdate()<<" Process time "<<double(end-start)/1000000<<" sec"<<endl;
	outlog<<seqsandinf._multiseqid.size()<<" sequences successfully imported"<<endl<<endl; 
	input.close();
	return 1 ;
}

void ImpBLOSUM(const char* _fpath,vector<darray<int> > & BLOSUMvec)
{
	BLOSUMvec.resize(100);
	ifstream INPUT(_fpath);
	if(!INPUT)
	{
		cout<<_fpath<<" not found, try to run BLOSUM or ./BLOSUM first to get BLOSUM scoring matrix"<<endl;
		exit(0);
	}
	char a;
	int id,sc;
	string str;
	darray<int> tmpvec;
	while(INPUT.get(a))
	{
		if(a=='>')
		{
			INPUT>>id;
			getline(INPUT,str);
			INPUT>>str;
			if(str=="Failed")
			{
				BLOSUMvec[id-1]=tmpvec;
				continue;
			}
			else
			{
				getline(INPUT,str);
				tmpvec.fast_resize(20,20);
				for(int i=0;i<20;++i)
				{
					for(int j=0;j<20;++j)
					{
						INPUT>>sc;
						tmpvec(i,j)=sc;
					}
				}
				BLOSUMvec[id-1]=tmpvec;
				tmpvec.fast_resize(0,0);
			}		
		}	
	}
}

void ImpMSA(const char* _fpath, vector<string> & infvec, darray<int> & combseq)
{
	ifstream IN(_fpath);
	int row,col,aa;
	IN>>row>>col;
	infvec.resize(row);
	combseq.fast_resize(row,col);	

	string str;
	getline(IN,str);	//remove newline character in first line
	for(int i=0;i<row;++i)
	{	getline(IN,str);
		infvec[i]=str;
	}
	for(int i=0;i<row;++i)
	{	for(int j=0;j<col;++j)
		{	IN>>aa;
			combseq(i,j)=aa;
		}
	}
}
		

