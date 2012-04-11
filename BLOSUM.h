//////////////////////////////////////////////////////////////
//                         BLOSUM                           //
//            STEVEN HENIKOFF AND JORJA G.HENIKOFF          //
//               Proc. Natl. Acad. Sci. USA                 //
//                Vol. 89, pp. 10915-10919                  //
//           Copyright Belongs to Ruan Xiaoyang             //
//                    ruansun@163.com                       //
//////////////////////////////////////////////////////////////
//This module creates blosum scoring matrix from standard block file obtained from http://blocks.fhcrc.org/blocks/uploads/blosum/.
//For each run,a log file will be saved in ./log indicating the block version and number of blocks contained

//block file format
/*********************************/

//ID   GLU_CARBOXYLATION; BLOCK
//AC   BL00011; distance from previous block=(1,64)
//DE   Vitamin K-dependent carboxylation domain proteins.
//BL   ECA motif; width=40; 99.5%=703; strength=2331

//the sequence should be formatted like 

//ACP2_SPIOL  (    82)  EIGADSLDTVEIVMKLEEEFGVTVEEENAQTITTIQ
//ACP_ECOLI  (     31)  DLGADSLDTVELVMALEEEFDTEIPDEEAEKITTVQ
//ACP_NEUCR  (     86)  DLGLDSLDTVEVVMAIEEEFSIEIPDKDADQIHSVD
//and ended with "//"
//**note:DATA FILE IS "newline sensitive,case sensitive" AND "blank insensitive"

//Protein information was stored in _infvec.
//Block information was stored in _blkfile


struct inf_blk                           //structure used to hold information of whole block and title information
{
	string  _blockinfo;
	vector<darray<char> > _blkfile;
};

//Menu:[BLOSUM]-[_BLOSUMsingleblock]-[_CprsMatx]
//             -[_Breakstring]
//     [_Impblkfile]-[_Breakstring]


//Menu:[_Impblkfile]
//Import block file for BLOSUM use
//Update 1 value."_blkandinf" block and information structure
//"_defau" whether use the default block file (N/Y)
void _Impblkfile(inf_blk & _blkandinf,char* _fpath)               //import block file from harddisk
{
	clock_t start,end;
	start=clock();

	char answer,a;
	ifstream ifblkdb;
	bool entseq=false;
	sarray<char> seq;
	darray<char> block;
	string tempstr,ch_a="Blocks";
	vector<sarray<char> > checkarr;
	const char* fpath="blocks.dat";
	if(_fpath)
		fpath=_fpath;
	ifblkdb.open(fpath);
	if(!ifblkdb.is_open())							//check the existence of blocks.data file
	{
		cout<<"Can not find "<<fpath<<" in execution directory, type \"BLOSUM -h\" for help"<<endl;
		exit(0);
	}
/////////////////////check database validity
	
	ifblkdb.seekg(0,ios::beg);			//set to the start position
	getline(ifblkdb,tempstr);			//obtain the title line of block file
	checkarr=_Breakstring(tempstr);			//format the string line into vector<sarray<char>> type to make the following comparison easier
	if(!(checkarr[0]==ch_a))				//check if the block file has correct title.if not,remind the user
	{
		cout<<"This block file has wrong title,this may(or may not) cause problem.Continue anyway?(Y/N):";
        	cin>>answer;	 
        	if(answer=='N'||answer=='n')
		{
			ifblkdb.close();
			exit(0);
		}
	}
	
	cout<<"Import block file..."<<endl;
	_blkandinf._blockinfo=tempstr;			//record title information
	while(ifblkdb.get(a))				//search for "ID"--the start marker of a block
	{
		if(a!='I')
			continue;
		else
		{
			ifblkdb.get(a);
			if(a=='D')
				break;
			else
				continue;
		}
	}
	sarray<int> sizeloader;		//a temp array to hold the column size of each block
	int size,n=0;

//////////////////import file to _blkandinf._blkfile
	while(ifblkdb.get(a))
	{
		if(a=='/' || a=='\n')
			getline(ifblkdb,tempstr);
		for(int i=0;i<4;++i)						//skip four lines 
			getline(ifblkdb,tempstr);  
		for(int i=0;i<tempstr.size();++i)
		{
			if(tempstr[i]=='=')
			{
				i+=1;
				while((int)tempstr[i]>=48 && (int)tempstr[i]<=57)
				{
					sizeloader.pushback((int)tempstr[i]-48);
					i+=1;
				}
				break;
			}
		}
		size=0;
		for(int i=0;i<sizeloader.size();++i)
		{
			if(sizeloader[i]!=0)
			   size+=(int)(sizeloader[i]*pow(10.0,sizeloader.size()-i-1));
		}
		sizeloader.clear();
	  	entseq=false;
		seq.resize(size);
		while(ifblkdb.get(a))
		{
			if(a=='/')                               
			{
				_blkandinf._blkfile.push_back(block);		//all blocks were stored in _blkandinf._blkfile
				block.clear();
				break;
			}
			if(entseq==false)
			{
				if(a==' '|| a=='\n')
					continue;
		    		if(a==')')  
					entseq=true;
			}
			else
			{
				if(a!=' ' && a!='\n')
				{
					seq[n]=a;		//store AA in seq until met newline 
				    n+=1;
				}
				else if(a=='\n')
				{
					n=0;
					block.push_row(seq);	//store sequence in block
					entseq=false;
				}
			}
		}
	}
	ifblkdb.close();
	end=clock();
	ofstream outlog("./log",ios_base::app);
	outlog<<getdate()<<" Process time "<<double(end-start)/1000000<<" sec"<<endl;
	outlog<<"Import block file from "<<fpath<<endl;
	for(int i=0;i<checkarr.size();++i)
		outlog<<checkarr[i]<<" ";
	outlog<<"Contain "<<_blkandinf._blkfile.size()<<" Blocks"<<endl<<endl;
}

//menu:[BLOSUM]-[_BLOSUMsingleblock]-[_CprsMatx]
//compress grouping matrix to provide nonredundant grouping information
//take 2 values."_grouparr" upper right nonsymmetrical pairwise correlation matrix of seqs in one block."_rownum" row number of previous matrix
//update "_grouparr" to provide grouping information
void _CprsMatx(darray<int> & _grouparr,const int & _rownum)
{
	for(int i=_rownum-1;i>=1;--i)	//start from lower right corner,search above row for occurance of correlation 
	{
		for(int m=i-1;m>=0;--m)
		{
	        	if(_grouparr(m,i)==1)
			{
				for(int j=i;j<_rownum;++j)
				{
					if(_grouparr(i,j)==1)
					{
						_grouparr(i,j)=0;
						_grouparr(m,j)=1;
					}
				}
				goto nexti;
			}
		}
		for(int j=i+1;j<_rownum;++j)
		{
			if(_grouparr(i,j)==1)
			{
			   for(int m=i-1;m>=0;--m)
			   {
				   if(_grouparr(m,j)==1)
				   {
					  for(int n=i;n<_rownum;++n)
				      {
					      if(_grouparr(i,n)==1)
					      {
						     _grouparr(i,n)=0;
						     _grouparr(m,n)=1;
					     }
				      }
					  goto nexti;
					}
			   }
			}
		}
nexti:;
	}
}

//menu:[BLOSUM]-[_BLOSUMsingleblock]
//calculate mutation frequency count for each block
//take 3 values."_cvgstd" converge standard. "_sngblk" single block. "_simil" prepared pairwise similarity matrix
//update 1 value."Freqcnt" a 20x20 AA mutation frequency count matrix 
int _BLOSUMsingleblock(const int & _cvgstd,darray<char> & _sngblk,darray<double> & _simil,darray<double> & Freqcnt)
{
	int rownum=_sngblk.getrnum(),colnum=_sngblk.getcnum(),mark=0,lastrow;
	darray<int> grouparr(rownum,rownum,0),groupinf;
	darray<double> aagrp(rownum,2,0.0);
	grouparr(rownum-1,rownum-1)=1;
	for(int i=0;i<rownum-1;++i)
	{
		grouparr(i,i)=1;
		for(int j=i+1;j<rownum;++j)
		{
			if(100*_simil(i,j) >= _cvgstd)
				grouparr(i,j)=1;
		}
	}
	_CprsMatx(grouparr,rownum);	//see [_CprsMatx]
	for(int i=0;i<rownum;++i)	//find sequences grouped together,store in darray<int> groupinf
	{
		if(grouparr(i,i)==1) 
		{
			groupinf.push_row(1);
			groupinf(groupinf.getrnum()-1,0)=i;
			for(int j=i+1;j<rownum;++j)
			{
			    	if(grouparr(i,j)>0)
					groupinf.push_to_row(groupinf.getrnum()-1,j);
			}
		}
	}
	if(groupinf.getrnum()<=1)	//if all seqs converged,return 0.which means this block contributed no information
		return 0;
	
	lastrow=0;			//reform groupinf to aagrp format for each of operation
	for(int i=0;i<groupinf.getrnum();++i)
	{
		for(int j=0;j<groupinf.getcnum(i);++j)
		{
			aagrp(lastrow,0)=1.0/groupinf.getcnum(i);	//1st col, contribution of the AA
			aagrp(lastrow,1)=i;				//2nd col, group of AA
			lastrow+=1;
		}
	}

	sarray<int> aacol(rownum);
	for(int k=0;k<colnum;++k)
	{
	   int iter=0;
	   for(int i=0;i<groupinf.getrnum();++i)
	   {
		for(int j=0;j<groupinf.getcnum(i);++j)
		   {                      
			   aacol[iter]=aaid(_sngblk(groupinf(i,j),k));	//amino acid integer ID array for a column
			   iter+=1;
		   }
	   }

	   for(int i=0;i<rownum-1;++i)					//pairwise mutation count
	   {
		   for(int j=i+1;j<rownum;++j)
		   {
			  if(aagrp(i,1)==aagrp(j,1))			//two AAs in same group do not provide information
				   continue;
			  else
			       Freqcnt((int)aacol[i],(int)aacol[j])+=aagrp(i,0)*aagrp(j,0);	//update Freqcnt matrix
		   }
	   }
	}
	return 1;
}

//calculate pairwise sequence similarity for each block
//take 1 value. "_blocks" a vector of darray<char> containing sequences in a block
//update 1 value. "simil" a vector of darray<double> containing sequence pairwise similarity
void _calsimil(vector<darray<char> > & _blocks,vector<darray<double> > & simil)
{
	int seqnum,colnum;
	for(int i=0;i<_blocks.size();++i)
	{
		darray<char> & tmpref=_blocks[i];
		seqnum=tmpref.getrnum();
		colnum=tmpref.getcnum();
		darray<double> matx(seqnum,seqnum,0);
		for(int sa=0;sa<seqnum-1;++sa)
		{
			for(int sb=sa+1;sb<seqnum;++sb)
			{
				int iden=0;
				for(int pos=0;pos<colnum;++pos)
				{
					if(tmpref(sa,pos)==tmpref(sb,pos))
						iden+=1;
				}
				matx(sa,sb)=double(iden)/double(colnum);
			}
		}
		simil[i]=matx;
	}
}

//calcuate mutation frequency for each convergence standard (1-100)
//take 2 values."_blocks" the vector for sequence blocks."_similmatxvec" prepared similarity matrix for all blocks
//update 1 value."freqcntvec" vector for frequency count on each convergence level
//return a sarray containing the number of blocks used on each convergence level
sarray<int> _batchcvgstd(vector<darray<char> > & _blocks,vector<darray<double> > & _similmatxvec,vector<darray<double> > & freqcntvec) {
	freqcntvec.resize(100);
	int iter=0,blockused;
	sarray<int> blockn(100);	//an array storing the number of blocks used on each convergence level
	for(int cvg=1;cvg<=100;++cvg)
	{
		darray<double> freqcnt(20,20,0);
		blockused=0;
		for(int i=0;i<_blocks.size();++i)
		{
			blockused+=_BLOSUMsingleblock(cvg,_blocks[i],_similmatxvec[i],freqcnt);
		}
		cout<<"converge >"<<cvg<<"\tblock used "<<blockused<<endl;
		blockn[iter]=blockused;
		freqcntvec[iter]=freqcnt;
		iter+=1;
	}
	return blockn;
}

//calculate relative frequency and blosum score for certain frequency count table
//take 1 value."_freqcnt" the frequency count at certain convergence level
//update 2 values."RltvFreq" relative frequency,"score" blosum score
//return 0 if total frequency count is zero, or any AA has zero frequency count
int _singlerltvfreq(darray<double> & _freqcnt,darray<double> & RltvFreq,darray<int> & score)
{
	RltvFreq.fast_resize(20,20,0);
	score.fast_resize(20,20,0);
	darray<double> Freq(20,20,0.0),AA_AAProb(20,20,0.0); //"Freq"="freqcnt"/total frequency count."AA_AAProb" random frequency
	sarray<double> AAFreq(20,0.0);
	double Totalfreqcnt=0,temp=0;
	int zerofreq=0;
	for(int i=0;i<20;++i)
	{
		for(int j=i;j<20;++j)
		{
			if(j!=i)
			{
				_freqcnt(i,j)+=_freqcnt(j,i);                        //add A->B to B->A mutation 
			    	_freqcnt(j,i)=_freqcnt(i,j);    
				if(_freqcnt(i,j)==0)
					zerofreq+=1;
				Totalfreqcnt+=_freqcnt(i,j);                        //add up total frequency count
			}
			if(j==i)
				Totalfreqcnt+=_freqcnt(i,i);
		}
	}
	if((int)Totalfreqcnt==0)
		return 0;
	for(int i=0;i<20;++i)
	{
		for(int j=0;j<20;++j)
		{
			if(i<=j)
				Freq(i,j)=_freqcnt(i,j)/Totalfreqcnt;
			if(i==j)
                		temp+=_freqcnt(i,j);                                //temp was used to calculate AA frequency
			else
				temp+=_freqcnt(i,j)/2.0;
		}
		AAFreq[i]=temp/Totalfreqcnt;	//calculate AA frequency
		if(temp==0)
			return 0;
		temp=0;
	}
	for(int i=0;i<20;++i)
	{
		for(int j=i;j<20;++j)
		{
			if(i==j)
   				AA_AAProb(i,j)=AAFreq[i]*AAFreq[j];
			else
				AA_AAProb(i,j)=2*AAFreq[i]*AAFreq[j];
			RltvFreq(i,j)=Freq(i,j)/AA_AAProb(i,j);
			double sc=-4;			//value with 0 frequency will have its score set to -4
			if(RltvFreq(i,j)>0.00001)
				sc=2*(log(RltvFreq(i,j))/log(2));//log2 based and multiply by 2 (according to standard algorithm)
			if(sc>0)		//round to the nearest integer
				sc+=0.5;
			else
				sc-=0.5;	
			score(j,i)=score(i,j)=int(sc);
		}
	}
	return 1;
}

//calculate relative frequency and blosum score for a vector of frequency count
//take 1 value."_freqcntvec" a vector of frequency count
//update 3 values."rltvfreqvec" a vector of relative frequency,"scorevec" a vector of score,"succflag" an array indicating whether relative frequency was successfully obtained
void _batchrltvfreq(vector<darray<double> > & _freqcntvec,vector<darray<double> > & rltvfreqvec,vector<darray<int> > & scorevec,sarray<int> & succflag)
{
	rltvfreqvec.resize(_freqcntvec.size());
	scorevec.resize(_freqcntvec.size());
	succflag.resize(_freqcntvec.size());
	darray<double> rltvfreq;
	darray<int> score;
	for(int i=0;i<_freqcntvec.size();++i)
	{
		succflag[i]=_singlerltvfreq(_freqcntvec[i],rltvfreq,score);
		rltvfreqvec[i]=rltvfreq;
		scorevec[i]=score;
	}
}

//the main function
//take 1 value,"_defau" [Y/N] whether to use the default block data
//output two files to disk."./relative_frequency.txt" the relative frequency used to calculate blosum score."./BLOSUM_score.txt" the transformed blosum score.
 
void CreateBlosumScoreMatrix(char* _fpath=NULL,const char _oprltv='Y',char * _opath=NULL)
{
	inf_blk _blkandinf;
	_Impblkfile(_blkandinf,_fpath);
	int blocknum=_blkandinf._blkfile.size();
	int seqnum=0,colnum=0,AAnum=0;
	sarray<char> aaorder=aa20();
	for(int i=0;i<blocknum;++i)
	{
		seqnum+=_blkandinf._blkfile[i].getrnum();	//record number of sequence in each block
		colnum+=_blkandinf._blkfile[i].getcnum(0);	//record number of column in each block
		AAnum+=_blkandinf._blkfile[i].getrnum()*_blkandinf._blkfile[i].getcnum(0);	//record total AA number in each block
	}

	vector<darray<double> > similmatxvec;		//compute similarity matrix for each block and store in vector
	similmatxvec.resize(blocknum);
	_calsimil(_blkandinf._blkfile,similmatxvec);
	
	vector<darray<double> > freqcntvec,rltvfreqvec;	//compute frequency count for each convergence standard
	vector<darray<int> > scorevec;
	sarray<int> blockn,succflag;
	blockn=_batchcvgstd(_blkandinf._blkfile,similmatxvec,freqcntvec);

	_batchrltvfreq(freqcntvec,rltvfreqvec,scorevec,succflag);
	
	ofstream log("./log",ios_base::app);
	log<<getdate()<<" Create custom BLOSUM relative frequency matrix"<<endl;
	log<<"block number:"<<blocknum<<" sequence number:"<<seqnum<<" column number:"<<colnum<<" AAnum:"<<AAnum<<endl;
	
	const char* outpath="./BLOSUM_score.txt";	//the default scoring matrix output file
	if(_opath)					//use custom output file if it exists
		outpath=_opath;	

	ofstream recscore(outpath);//record converge standard,block file information and relative frequency matrix
	for(int i=0;i<succflag.size();++i)
	{
		recscore<<endl<<">"<<i+1<<"\t"<<blockn[i]<<"\t"<<_blkandinf._blockinfo<<endl;
		if(succflag[i]==0)
		{
			recscore<<"Failed because zero frequency count"<<endl;		
			continue;
		}
		recscore<<aaorder<<endl;	
		scorevec[i].record(recscore);
	}
	if(_oprltv=='Y')				//output relative frequency file by default
	{
		ofstream recrltv("./RltvFreq.txt");
		for(int i=0;i<succflag.size();++i)
		{
			recrltv<<endl<<">"<<i+1<<"\t"<<blockn[i]<<"\t"<<_blkandinf._blockinfo<<endl;
			if(succflag[i]==0)
			{
				recrltv<<"Failed because zero frequency count"<<endl;		
				continue;
			}
			recrltv<<aaorder<<endl;
			rltvfreqvec[i].record(recrltv);
		}
	}
}
