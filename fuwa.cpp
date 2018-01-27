/*
 * Programme: Fuwa
 * Description: A decision-tree-based variant caller
 * Version: v1.0
 * Update: 2018.01.21
 *
 */

#include <algorithm>
#include <iostream>
#include <stdint.h>
#include <cstdlib>
#include <fstream>
#include <cfloat>
#include <vector>
#include <cmath>
#include <ctime>
#include <sam.h>
#include <zlib.h>
#include <map>
#include <unistd.h>
using	namespace	std;

struct	Site{
	uint64_t	k;
	unsigned	rs;
	bool	operator()(Site	X,	Site	Y){	return	X.k<Y.k;	}
}__attribute__ ((packed));

struct	Variant{
	uint64_t	k;
	float	qual;
	unsigned	char	dbsnp,	rdep,	adep;
	string	indel;
	static	bool	by_k(Variant	X,	Variant	Y){	return	X.k<Y.k;	}
	static	bool	by_q(Variant	X,	Variant	Y){	return	X.qual>Y.qual;	}
};

uint8_t	table[16];
uint8_t	code[256];
double	phred[256];
string	ref_dir = "";
int	cur_chr=-1;
uint64_t	pileupn=0;
string	ref;
samfile_t	*bam;
vector<Site>	site;
bool male = false;
float qualThreshold = 0.6;
float snpQualThreshold = 0.8;
vector<Variant>	variant;
vector<float>	attr;
string	attr_name[]={
	"Depth",	"Ratio",	"DeltaL",
	"SumBQ",	"AveBQ",
	"AveMQ",	"WorMQ",	"PoorMQ",
	"VarPos",	"VarStr",	"AveAS",
	"BiasStr",
};
const	size_t	attrn=sizeof(attr_name)/sizeof(string);
const	double	l01=log(0.01),	l99=log(0.99),	l50=log(0.5);
double	tree_like=0,	min_split;
size_t	max_bin;
FILE	*ftree;

void	make_code(void){
	for(size_t	i=0;	i<256;	i++){	phred[i]=1-exp10f(-0.1*i);	code[i]=4;	}
	for(size_t	i=0;	i<16;	i++)	table[i]=4;
	table[1]=0;	table[2]=1;	table[4]=2;	table[8]=3;
	code['A']=code['a']=0;	code['C']=code['c']=1;	code['G']=code['g']=2;	code['T']=code['t']=3;
}

bool	load_site(const	char	*F){
	gzFile	in=gzopen(F,	"rb");
	if(in==Z_NULL)	return	false;
	unsigned	n;
	gzread(in,	&n,	4);
	site.resize(n);
	gzread(in,	&site[0],	site.size()*sizeof(Site));
	gzclose(in);
	return	true;
}

bool	load_ref(int	I){
	if(I==cur_chr)	return	true;
	cur_chr=I;
	string	fn=ref_dir+bam->header->target_name[I]+".fa";
	FILE	*f=fopen(fn.c_str(),	"rt");
	if(f==NULL){	cerr<<"fail to open "<<fn<<'\n';	return	false;	}
	ref.clear();	ref.reserve(bam->header->target_len[I]);
	char	buff[1024];
	while(fgets(buff,	1024,	f)!=NULL)	if(buff[0]!='>')	for(char	*p=buff;	*p>=33;	p++)	ref.push_back(*(code+*p));
	fclose(f);
	if(ref.size()!=bam->header->target_len[I]){	
		cerr<<fn<<'\t'<<ref.size()<<"\t!=\t"<<bam->header->target_name[I]<<'\t'<<bam->header->target_len[I]<<'\n';	return	false;	
	}
	return	true;
}

double	LRT(double	A,	double	B,	double	C,	double	D){
	if(A+B<=0||A+C<=0||B+D<=0||C+D<=0)	return	0;
	double	n=A+B+C+D,	p=(A+B)/n,	q=(A+C)/n,	l=0;
	if(A>0)	l+=A*logf(A/(p*q*n));
	if(B>0)	l+=B*logf(B/(p*(1.0-q)*n));
	if(C>0)	l+=C*logf(C/((1.0-p)*q*n));
	if(D>0)	l+=D*logf(D/((1.0-p)*(1.0-q)*n));
	return	l;
}

int	pileup_func(uint32_t tid,	uint32_t pos,	int n,	const	bam_pileup1_t	*pl,	void *data){
	if(++pileupn%100000000==0)	cerr<<'=';
	char	*chr=bam->header->target_name[tid];
	if(data!=NULL||strlen(chr)>2||(!male&&chr[0]=='Y'))	return	0;
	load_ref(tid);	
	unsigned	ra=ref[pos];
	if(ra==4)	return	0;
//	bool	Haploid=strcmp(chr,	"Y")==0||(male&&strcmp(chr,	"X")==0)||strcmp(chr,	"MT")==0;
	unsigned	cnt[6]={0,0,0,0,0,0};
	double	ebd[6]={0,0,0,0,0,0};
	map<string,	unsigned>	ins,	del;
	{
		string	s;
		const	bam_pileup1_t	*p=pl;
		for(int	i=0;	i<n;	i++,	p++){
			float	wmap=phred[p->b->core.qual];
			if(p->indel>0){	
				cnt[4]++;	ebd[4]+=wmap;
				s.clear();
				for(int	j=0;	j<=p->indel;	j++)	s.push_back("ACGTN"[table[bam1_seqi(bam1_seq(p->b), p->qpos+j)]]);
				ins[s]++;
			}
			else	if(p->indel<0){	
				cnt[5]++;	ebd[5]+=wmap;	
				s.clear();
				for(int	j=0;	j<=-p->indel;	j++)	s.push_back("ACGTN"[(int)ref[pos+j]]);
				del[s]++;
			}
			else	if(!p->is_del){
				uint8_t	c=table[bam1_seqi(bam1_seq(p->b), p->qpos)];
				if(c<4){	cnt[c]++;	ebd[c]+=wmap*phred[bam1_qual(p->b)[p->qpos]];	}
			}
		}
	}
	double	sebd=ebd[0]+ebd[1]+ebd[2]+ebd[3]+ebd[4]+ebd[5];
	for(size_t	aa=0;	aa<6;	aa++)	if(ebd[aa]>sebd/11&&aa!=ra){
		double	Depth=ebd[aa],	Ratio=ebd[aa]/sebd,	SumBQ=0;
		double	SumMQ=0,	WorMQ=100,	PoorMQ=0,	AvePos=0,	VarPos=0,	AveAS=0;
		char	tagAS[2]={'A','S'};
		double	Str[2][2]={{0,0},{0,0}};
		const	bam_pileup1_t	*p=pl;
		for(int	i=0;	i<n;	i++,	p++){
			if(p->b->core.qual<WorMQ)	WorMQ=p->b->core.qual;
			if(p->b->core.qual<10)	PoorMQ+=1;
			double	as=bam_aux2i(bam_aux_get(p->b,	tagAS));
			if(p->indel>0){
				if(aa==4){
					SumBQ+=30;
					SumMQ+=p->b->core.qual;
					AvePos+=p->qpos;	VarPos+=p->qpos*p->qpos;
					AveAS+=as;
					Str[1][bam1_strand(p->b)]+=1;
				}
				else{
					Str[0][bam1_strand(p->b)]+=1;
				}
			}
			else	if(p->indel<0){	
				if(aa==5){
					SumBQ+=30;
					SumMQ+=p->b->core.qual;
					AvePos+=p->qpos;	VarPos+=p->qpos*p->qpos;
					AveAS+=as;
					Str[1][bam1_strand(p->b)]+=1;
				}
				else{
					Str[0][bam1_strand(p->b)]+=1;
				}
			}
			else	if(!p->is_del){
				uint8_t	c=table[bam1_seqi(bam1_seq(p->b), p->qpos)];
				if(c==aa){	
					SumBQ+=bam1_qual(p->b)[p->qpos];
					SumMQ+=p->b->core.qual;
					AvePos+=p->qpos;	VarPos+=p->qpos*p->qpos;
					AveAS+=as;
					Str[1][bam1_strand(p->b)]+=1;
				}
				else{
					Str[0][bam1_strand(p->b)]+=1;
				}
			}
		}
		double	DeltaL=max(ebd[ra]*l01+ebd[aa]*l99,	ebd[ra]*l50+ebd[aa]*l50);
		if(ebd[ra])	DeltaL-=ebd[ra]*logf(ebd[ra]/(ebd[ra]+ebd[aa]));
		if(ebd[aa])	DeltaL-=ebd[aa]*logf(ebd[aa]/(ebd[ra]+ebd[aa]));
		DeltaL/=ebd[aa]+ebd[ra];
		Variant	v;
		if(aa==4){
			size_t	m=0;
			for(map<string,	unsigned>::iterator	mi=ins.begin();	mi!=ins.end();	++mi)	if(mi->second>m){	m=mi->second;	v.indel=mi->first;	}
		}
		else	if(aa==5){
			size_t	m=0;
			for(map<string,	unsigned>::iterator	mi=del.begin();	mi!=del.end();	++mi)	if(mi->second>m){	m=mi->second;	v.indel=mi->first;	}
		}
		v.k=((uint64_t)chr[0]<<56)|((uint64_t)chr[1]<<48)|((uint64_t)(pos+1)<<16)|((uint64_t)("ACGTID"[ra])<<8)|"ACGTID"[aa];
		{
			Site	s;
			s.k=((uint64_t)chr[0]<<56)|((uint64_t)chr[1]<<48)|((uint64_t)(pos+1)<<16)|((uint64_t)'N'<<8)|"ACGTID"[aa];
			size_t	k=lower_bound(site.begin(),	site.end(),	s,	Site())-site.begin();
			v.dbsnp=k<site.size()&&site[k].k==s.k;
			if(ebd[ra]<255&&ebd[aa]<255){	v.rdep=ebd[ra]+0.5;	v.adep=ebd[aa]+0.5;	}
			else{	v.rdep=ebd[ra]/max(ebd[ra],	ebd[aa])*255+0.5;	v.adep=ebd[aa]/max(ebd[ra],	ebd[aa])*255+0.5;	}
		}
//		if(Haploid&&ebd[aa]>sebd/11&&ebd[aa]<(sebd*10)/11)	continue;
		variant.push_back(v);
		attr.push_back(Depth);	attr.push_back(Ratio);	attr.push_back(DeltaL);
		attr.push_back(SumBQ);	attr.push_back(SumBQ/cnt[aa]);
		attr.push_back(SumMQ/cnt[aa]);	attr.push_back(WorMQ);	attr.push_back(PoorMQ/n);
		attr.push_back(VarPos-(AvePos/cnt[aa])*AvePos);	attr.push_back(Str[1][0]*Str[1][1]/(Str[1][0]+Str[1][1]));	attr.push_back(AveAS/cnt[aa]);
		attr.push_back(LRT(Str[0][0],	Str[0][1],	Str[1][0],	Str[1][1]));
	}
	return	0;
}

struct	Sort{
	float	x,	y;
	bool	operator()(Sort	X,	Sort	Y){	return	X.x<Y.x;	}
};

struct Node{
	float	mean,	threshold;
	uint32_t	variable;
	Node	*left;
	Node	*right;
};

struct	Decision{
	uint32_t	variable;
	float	threshold;
	bool	greater;
};

double	choose(vector<uint32_t>	&I,	uint32_t	&X,	float	&T,	float	&M){
	size_t	n=I.size();
	double	sy=0;
	for(size_t	i=0;	i<n;	i++)	sy+=variant[I[i]].dbsnp;
	M=(sy+0.5)/(n+1);
	vector<float>	threshold(attrn);	vector<Sort>	vs(n);
	double	maxg=-FLT_MAX;	X=0;
	for(size_t	f=0;	f<attrn;	f++){
		for(size_t	i=0;	i<n;	i++){	vs[i].x=attr[I[i]*attrn+f];	vs[i].y=variant[I[i]].dbsnp;	}
		sort(vs.begin(),	vs.end(),	Sort());
		double	xsy=0,	best=-FLT_MAX;
		for(size_t	i=0;	i<n;	i++){
			xsy+=vs[i].y;
			double	g=LRT(i+1-xsy,	n-1-i-(sy-xsy),	xsy,	sy-xsy);
			if(g>best&&(i==n-1||vs[i].x!=vs[i+1].x)){	best=g;	threshold[f]=vs[i].x;	}
		}
		if(best>maxg){	maxg=best;	X=f;	}
	}
	T=threshold[X];
	return	maxg;
}

void	add_node(Node	**N,	vector<uint32_t>	&I,	vector<Decision>	&D){
	*N=new	Node;
	double	gain=choose(I,	(*N)->variable,	(*N)->threshold,	(*N)->mean);
	if(gain>min_split||(gain>0&&I.size()>max_bin)){
		vector<uint32_t>	left,	right;
		vector<Decision>	leftd=D,	rightd=D;
		Decision	d;	d.variable=(*N)->variable;	d.threshold=(*N)->threshold;
		d.greater=false;	leftd.push_back(d);	d.greater=true;	rightd.push_back(d);
		for(size_t	i=0;	i<I.size();	i++){
			if(attr[I[i]*attrn+(*N)->variable]<=(*N)->threshold)	left.push_back(I[i]);
			else	right.push_back(I[i]);
		}
		vector<uint32_t>().swap(I);	vector<Decision>().swap(D);
		if(left.size())	add_node(&((*N)->left),	left,	leftd);	else	(*N)->left=NULL;
		if(right.size())	add_node(&((*N)->right),	right,	rightd);	else	(*N)->right=NULL;
	}
	else{
		(*N)->left=(*N)->right=NULL;
		fprintf(ftree,	"%u\t%g",	(unsigned)I.size(),	(*N)->mean);
		for(size_t	i=0;	i<D.size();	i++)	fprintf(ftree,	"\t%s%s%g",	attr_name[D[i].variable].c_str(),	D[i].greater?">":"<=",	D[i].threshold);
		fprintf(ftree,	"\n");
		for(size_t	i=0;	i<I.size();	i++)	variant[I[i]].qual=(*N)->mean;
		tree_like+=I.size()*(*N)->mean*log((*N)->mean)+I.size()*(1-(*N)->mean)*log(1-(*N)->mean);
	}	
}

void free_node(Node	*N){
	if(N!=NULL){free_node(N->left);	free_node(N->right); delete	N;}
}

void one_tree(void){
	cerr<<"\ncandidates\t"<<variant.size()<<'\n';
	min_split=log(variant.size())*2.25;
	//max_bin=variant.size()/log(variant.size());
        max_bin=100000;
	Node *root=NULL;
	vector<uint32_t> id(variant.size());
	for(size_t	i=0; i<variant.size(); i++) id[i]=i;
	vector<Decision> d;
	add_node(&root,	id,	d);
	free_node(root);
}

void	vcf(const	char	*F){
	sort(variant.begin(),	variant.end(),	Variant::by_k);
	string	fn=F;	fn+=".vcf";
	FILE	*fo=fopen(fn.c_str(),	"wt");
	fprintf(fo,	"##fileformat=VCFv4.0\n");
	fprintf(fo,	"##source=Fudan:wgstools:fuwa\n");
	fprintf(fo,	"##reference=1000Genomes-NCBI37\n");
	fprintf(fo,	"##INFO=<ID=CONF,Number=1,Type=String,Description=\"Confident or suspicious\">\n");
	fprintf(fo,	"##INFO=<ID=TYPE,Number=1,Type=String,Description=\"SNP or indel\">\n");
	fprintf(fo,	"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	fprintf(fo,	"##FORMAT=<ID=ED,Number=1,Type=String,Description=\"EffectiveDepth, Ref to Alt\">\n");
	fprintf(fo,	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n",	F);
	for(size_t	i=0;	i<variant.size();	i++){
		bool isIndel = (variant[i].k&255)=='I'||(variant[i].k&255)=='D';
		if(isIndel && variant[i].qual < qualThreshold) {continue;}
		if(!isIndel && variant[i].qual < snpQualThreshold) {continue;}
		unsigned	rs;
		{
			Site	s;
			s.k=((variant[i].k>>16)<<16)|((uint64_t)'N'<<8)|(variant[i].k&255);
			size_t	k=lower_bound(site.begin(),	site.end(),	s,	Site())-site.begin();
			rs=k<site.size()&&site[k].k==s.k?site[k].rs:0;
		}
		fputc(variant[i].k>>56,	fo);
		if((variant[i].k>>48)&255)	fputc((variant[i].k>>48)&255,	fo);
		fprintf(fo,	"\t%u\t",	(unsigned)((variant[i].k>>16)&0xffffffff));
		if(rs){	fprintf(fo,	"rs%u",	rs);	variant[i].dbsnp=1;	}
		else{	fprintf(fo,	".");	variant[i].dbsnp=0;	}
		if((variant[i].k&255)=='I')	fprintf(fo,	"\t%c\t%s\t",	(char)((variant[i].k>>8)&255),	variant[i].indel.c_str());
		else	if((variant[i].k&255)=='D')	fprintf(fo,	"\t%s\t%c\t",	variant[i].indel.c_str(),	(char)((variant[i].k>>8)&255));
		else	fprintf(fo,	"\t%c\t%c\t",	(char)((variant[i].k>>8)&255),	(char)(variant[i].k&255));
		fprintf(fo,	"%g\tPASS\t",	variant[i].qual);
		fprintf(fo,	"CONF=%s",	variant[i].qual<0.9?"suspicious;":"confident;");
		// fprintf(fo,	"TYPE=%s",	(variant[i].k&255)=='I'||(variant[i].k&255)=='D'?"INDEL;":"SNP;");
		fprintf(fo,	"TYPE=%s",	isIndel?"INDEL;":"SNP;");
		fprintf(fo,	"\tGT:ED\t%s:%u/%u\n",	variant[i].adep>variant[i].rdep*10?"1/1":"0/1",	variant[i].rdep,	variant[i].adep);
	}
	fclose(fo);
	// char	cmd[256];
	// sprintf(cmd,	"bgzip -f %s.vcf",	F);
	// if(system(cmd)==-1)	cerr<<"bgzip failed\n";
	// sprintf(cmd,	"tabix -p vcf %s.vcf.gz",	F);
	// if(system(cmd)==-1)	cerr<<"tabix failed\n";
}

void printUsage() {
	printf("Usage: fuwa [options]\n\n"
			"  -i input     input bam file\n"
			"  -d dbSNP     dbSNP gz file\n"
			"  -r ref_dir   reference directory\n"
			"  -o output    [optional] output file name without extension (default: input file name)\n"
			"  -s SNP_qual  [optional] Filtering threshold of quality score for SNPs (range: [0, 1], default: 0.8)\n"
			"  -q qual      [optional] Filtering threshold of quality score for indels (range: [0, 1], default: 0.6)\n"
			"  -m           [optional] the sample is male. By default the sample is considered female\n\n"
			);
}

int	main(int	ac,	char	**av){
	cerr<<"*******************\n";
	cerr<<"*       Fuwa      *\n";
	cerr<<"* Author: Yi Wang *\n";
	cerr<<"* Version: 1.0    *\n";
	cerr<<"*******************\n";
	int ch;
	opterr = 0;
	char *inputFile = NULL;
	char *outputFileName = NULL;
	char *dbSNP = NULL;
	while((ch = getopt(ac, av, "i:d:r:o:q:m")) != -1) {
		switch(ch) {
			case 'i':
				printf("Input file: %s\n", optarg);
				inputFile = optarg;
				break;
			case 'o':
				printf("Output file: %s.vcf %s.tree\n", optarg, optarg);
				outputFileName = optarg;
				break;
			case 'd':
				printf("dbSNP: %s\n", optarg);
				dbSNP = optarg;
				break;
			case 'r':
				printf("Ref_dir: %s\n", optarg);
				ref_dir = optarg;
				if(*ref_dir.rbegin()!='/')	ref_dir += '/';
				break;
			case 'q':
				printf("Indel qual threshold: %s\n", optarg);
				qualThreshold = atof(optarg);
				break;
			case 's':
				printf("SNP qual threshold: %s\n", optarg);
				snpQualThreshold = atof(optarg);
				break;
			case 'm':
				printf("Male\n");
				male = true;
				break;
			default:
				printUsage();
				return 1;
		}
	}
	// if(ac<6){	cerr<<"fuwacall input.bam male(1/0) dbsnp.gz ref_dir/ qual output\n";	return	0;	}
	// male=atoi(av[2])>0;	ref_dir=av[4];	if(*ref_dir.rbegin()!='/')	ref_dir+='/'; qualThreshold=atof(av[5]); 
	if(inputFile == NULL || dbSNP == NULL || ref_dir.length() == 0) {
		printUsage();
		return 1;
	}
	if(!load_site(dbSNP)){	cerr<<"fail to load "<< dbSNP <<'\n';	return	0;	}
	make_code();
	bam=samopen(inputFile,	"rb",	0);
	if(bam==0){	cerr<<"fail to open "<< inputFile <<'\n';	return	0;	}
	sampileup(bam,	-1,	pileup_func,	NULL);
	samclose(bam);
	// cerr<< outputFileName <<endl;
	string	fn = outputFileName ? outputFileName : inputFile;
	// fn+=".tree";
	ftree=fopen((fn + ".tree").c_str(),	"wt");
	one_tree();
	fclose(ftree);
	vcf(outputFileName?outputFileName:inputFile);
	return	0;
}
