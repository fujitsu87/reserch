#include <math.h>
#include "ldpc_structure.c"


double Lanhpro(struct LIST *p,unsigned int n,double tc[n],struct LIST *La){
	struct LIST	*q,*a;
	double	product=1.0;
		q=p;
		a=La;
		while(q!=NULL && a!=NULL){
			product*=tanh((tc[q->INDEX]+a->VALUE)*0.5);
			q=q->NEXT;
			a=a->NEXT;
		}
		return product;
}
double Tanhpro(unsigned int n,double Tc[n],struct LIST *Ta,int dc,double Tac[dc]){
	int i=0;
	double	value;
	struct LIST	*a;
	double	product=1.0;

		a=Ta;

		while(a!=NULL){
			value=(a->VALUE+Tc[a->INDEX])/(1+a->VALUE*Tc[a->INDEX]);
			if(Tac[i]<=-1.0){
				Tac[i]=-(1.0-1e-10);
			}
			else if(value>=1.0){
				Tac[i]=(1.0-1e-10);
			}

			else {
				Tac[i]=value;
			}
			a=a->NEXT;
			product*=Tac[i];
			i++;
		}
		return product;
}


double	tanhaplusb(double a,double b,int Thr,double apparea){
	double tanhab;
	double	t;
	//差の総和が0になり，両端で一致する近似曲線
	double a1=0.34111691664032814349660727125094;
	double a2=-0.84111691664032814349660727125094;
//差の二乗総和が最小になる近似曲線
//	double a1=0.37393861172567665174583033747805;
//	double a2=-0.85989550047409341706106911748307;
	//力業で解いた近似
//		double a1=0.3715990;
//		double a2=-0.8619350;

	t=a*b;

	 if(t==-1.0){
	 	tanhab=1.0-1e-10;
 	}
	 else{
		 	if(Thr==0){
				//tanhab=(a+b)*(1-t+pow(t,2)-pow(t,3)+pow(t,4)-pow(t,5)+pow(t,6)-pow(t,7)+pow(t,8));
				//
	 		tanhab=(a+b)/(1+t);
			//tanhab=(a+b)*(0.3715990*t*t-0.8619350*t+1.0);

	   }
	 	 else{
	// 		 	tanhab=(a+b)*(1-t+pow(t,2));
					if(apparea<-0.21){				//-0.3<x1.0
							 a1=3716182287.0/7141373213.0;
							 a2=-537742559.0/549336401.0;
					 }
					else if(apparea<-0.11){			//-0.2<x1.0
							a1=427441.0/1025280.0;
							a2=-22912351.0/25632000.0;
					}
					else if(apparea<-0.01){			//-0.1<x1.0
						a1=3211208519.0/8493133418.0;
						a2=-333240090.0/386051519.0;
					}
					else if(apparea<0.01){			//-0.0<x1.0
						a1=0.37393861172567665174583033747805;
						a2=-0.85989550047409341706106911748307;
					}
					else if(apparea<0.11){			//0.1<x1.0
						a1=5953.0/16914.0;
						a2=-32975047.0/39113625.0;
					}
					else if(apparea<0.21){			//0.2<x1.0
						a1=1389240853.0/3812077312.0;
						a2=-25379481.0/29781854.0;
					}

					else{
						a1=	
						a2=-0.85989550047409341706106911748307;
					}


				tanhab=(a+b)*(a1*t*t+a2*t+1.0);
	//			tanhab=(a+b)*(0.3715990*t*t-0.8619350*t+1.0);
//			tanhab=(a+b)*(0.5644*t*t-1.062125*t+0.99806);
			//	tanhab=(a+b)*(0.6963*t*t-1.1797*t+0.9838);
//					tanhab=(a+b)*(0.898*t*t-1.3395*t+0.9422);
		//tanhab=(a+b)*((8.0/27.0)*t*t-(20.0/27.0)*t+(17.0/18.0));
			}

		 }
//		printf("a=%f b=%f t=%f tanhab=%f  Thr=%f\n",a,b,t,tanhab,Thr);
		if(fabs(tanhab)>=1.0){
			tanhab=tanhab/fabs(tanhab)*(1.0-1e-10);
		}

	return tanhab;

}

void UpdateTe(double Tpro,struct LIST *c_Te,unsigned int dc,double Tac[dc]){
struct LIST *p;
int		i=0;
	p=c_Te;
	while(p!=NULL){
		p->VALUE=Tpro/Tac[i];
		if(fabs(p->VALUE)>=1.0){
			p->VALUE=(p->VALUE/fabs(p->VALUE))*(1.0-1e-10);
		}
		p=p->NEXT;
		i++;
	}

}
/*************************************************/
//関数名:Tanhpro2
//説明：各Tacの値を計算して格納
/*************************************************/
void Tanhpro2(unsigned int n,double Tc[n],struct LIST *Ta,int dc,double Tac[dc],double apparea,int THRESHOLD){
	int i=0;
	double	value;
	struct LIST	*a;

		a=Ta;

		while(a!=NULL){
			value=(a->VALUE+Tc[a->INDEX])/(1+a->VALUE*Tc[a->INDEX]);
			value=tanhaplusb(a->VALUE,Tc[a->INDEX],THRESHOLD,apparea);
			if(Tac[i]<=-1.0){
				Tac[i]=-(1.0-1e-10);
			}
			else if(value>=1.0){
				Tac[i]=(1.0-1e-10);
			}

			else {
				Tac[i]=value;
			}
			a=a->NEXT;
			i++;
		}
}
/*************************************************/
//関数名:UpdateTe2
//説明：外部値Teの更新
/*************************************************/
void UpdateTe2(struct LIST *c_Te,unsigned int dc,double Tac[dc]){
struct LIST *p;
double pro=1.0;
int		i=0,j=0;
	p=c_Te;

	for(i=0;i<dc;i++){
		for(j=0;j<dc;j++){
			if(i!=j){
				pro*=Tac[j];
			}
		}
		p->VALUE=pro;
		p=p->NEXT;
		pro=1.0;
	}

}
void TanhSumList(struct LIST *Te,double TeSum[2]){

	struct LIST *p;
	p=Te;
	TeSum[0]=1.0;
	TeSum[1]=1.0;

	while(p!=NULL){
		TeSum[0]*=1+p->VALUE;
		TeSum[1]*=1-p->VALUE;
		p=p->NEXT;
	}

}
// double TanhSumList2(struct LIST *Te, double Thr){//,FILE *file_pt){
//
// 	struct LIST *p;
// 	double	value;
// 	p=Te;
//
// 	value=p->VALUE;
// 	p=p->NEXT;
// 	//printf("sumlist2\n");
// 	while(p!=NULL){
// 		value=tanhaplusb(value,p->VALUE,Thr);//,file_pt);
// 		p=p->NEXT;
// 	}
// 	return value;
// }
double SumList(struct LIST *p){
	struct LIST	*q;
	double	sum=0.0;
		q=p;
		while(q!=NULL){
			sum+=q->VALUE;
			q=q->NEXT;

		}
		return sum;
}
double UpdateTa2(unsigned int dv ,struct LIST *Ta,struct LIST *Te,int Thr,double apparea){
	struct LIST *p,*q;
	double a,b,*tTe,*tTa,Tpro;
	int		i,j;
	p=Ta;
	q=Te;
	//printf("UpdateTa2\n");
	tTe=(double *)malloc(sizeof (double)*dv);
	if(tTa==NULL){
		printf("Memory cannot be allocated!!\n");	exit(1);
	}
	tTa=(double *)malloc(sizeof (double)*dv);
	if(tTa==NULL){
		printf("Memory cannot be allocated!!\n");	exit(1);
	}


	for(i=0;i<dv;i++){
		tTa[i]=0.0;
	}
	i=0;
	Tpro=0.0;
	while(q!=NULL){
		tTe[i]=q->VALUE;
		i++;
		q=q->NEXT;
	}
	for(i=0;i<dv;i++){
		Tpro=tanhaplusb(Tpro,tTe[i],Thr,apparea);
		for(j=0;j<dv;j++){
			if(i!=j){
				tTa[j]=tanhaplusb(tTa[j],tTe[i],Thr,apparea);//,file_pt);
			}
		}
	}
	i=0;
	q=Te;
	while(p!=NULL){
	//	printf("	Ta[%d]=%e",i,tTa[i]);
		p->VALUE=tTa[i];
		p->INDEX=q->INDEX;
		i++;
		q=q->NEXT;
		p=p->NEXT;
	}
//printf("\n");

	free(tTa);
	free(tTe);

	return Tpro;
}
void UpdateTa(struct LIST *Ta,double TeSum[2],struct LIST *Te){
	struct LIST *p,*q;
	double a,b;
	p=Ta;
	q=Te;

	while(p!=NULL && q!=NULL){
		a=TeSum[0]/(1+q->VALUE);
		b=TeSum[1]/(1-q->VALUE);
		p->VALUE=(a-b)/(a+b);
		p->INDEX=q->INDEX;
		#ifdef debug
			// printf("a=%e b=%e	a+b=%e,a-b=%e	",a,b,a+b,a-b);
			// printf("Ta=%e\n",p->VALUE);
		#endif
		p=p->NEXT;
		q=q->NEXT;

	}

}

void UpdateLa(unsigned int n,double sum_Le[n],struct LIST *v_La[n],struct LIST *v_Le[n]){
	struct LIST	*p,*q;
	int	 i;
		for(i=0;i<n;i++){
			p=v_La[i];
			q=v_Le[i];
			while(p!=NULL && q!=NULL){
				p->VALUE=sum_Le[i]-q->VALUE;
				p->INDEX=q->INDEX;
				p=p->NEXT;
				q=q->NEXT;
			}
		}
}
/****************************************************************************/
//関数名: UpdateTa3
//説明:Sumproduct復号のTaの更新関数
/****************************************************************************/
double UpdateTa3(struct LIST *v_Ta,struct LIST *v_Te,unsigned int dv,double apparea,int THRESHOLD){
	double *A,*B,*C; //Teの値の一次記憶変数
	struct LIST *q;
	int i;
	double Tpro;
	q=v_Te;

	A=(double *)malloc(sizeof (double)*(dv));
	if(A==NULL){
		printf("Memory cannot be allocated!!\n");	exit(1);
	}
	B=(double *)malloc(sizeof (double)*(dv-1));
	if(B==NULL){
		printf("Memory cannot be allocated!!\n");	exit(1);
	}
	C=(double *)malloc(sizeof (double)*(dv-1));
	if(C==NULL){
		printf("Memory cannot be allocated!!\n");	exit(1);
	}

	//Taの更新のための事前計算
	A[0]=q->VALUE;
	Tpro=q->VALUE;
//	printf("%e\n",q->VALUE);
	//B[dv-2]=q->VALUE;
	q=q->NEXT;
	i=1;
	while(q!=NULL){
//		printf("%d",i);
		A[i]=tanhaplusb(A[i-1],q->VALUE,THRESHOLD,apparea);
		//B[dv-2-i]=tanhaplusb(B[dv-1-i],q->VALUE,-1.0);
		//printf("%d %e\n",i,q->VALUE);
		C[i-1]=q->VALUE;
		i++;
		q=q->NEXT;
	}
	B[dv-2]=C[dv-2];
	for(i=0;i<=(dv-3);i++){
		B[dv-3-i]=tanhaplusb(B[dv-2-i],C[dv-3-i],THRESHOLD,apparea);
	}
	Tpro=A[dv-1];//tanhaplusb(A[0],B[0],-1.0);
	//Taの更新
	q=v_Ta;
	q->VALUE=B[0];
	//printf("	Ta[0]%e",q->VALUE);

	q=q->NEXT;
	for(i=1;i<(dv-1);i++){
		q->VALUE=tanhaplusb(A[i-1],B[i],THRESHOLD,apparea);
	//	printf("	Ta[1]%e",q->VALUE);
		q=q->NEXT;
	}
	q->VALUE=A[dv-2];
//	printf("	Ta[2]%e",q->VALUE);

//printf("\n");
// 	for(i=0;i<(dv-1);i++){
 //		printf("A[%d]=%e,B[%d]=%e,C[%d]=%e\n",i,A[i],i,B[i],i,C[i]);
 	//}
// exit(1);

	free(A);	free(B); free(C);

	return Tpro;
}




void Sumproduct(unsigned int m,unsigned int n,double RecieveDATA[n],double sigma2,unsigned int dc,unsigned int dv,struct LIST *c_neighbor[m],unsigned int Loop,unsigned int BE[Loop]){
	
	int i,j,k;
	double Tpro,divide,Lpro;
	struct LIST	**c_La,**v_La; //事前値
	struct LIST	**c_Le,**v_Le; //外部値
	double	*Lc; //チャネル地
	struct LIST *p,*q,*r;
	double	value,*LeSum;
	double	test;
	
	
		c_Le=LDPCInitialization(m,c_Le);
		c_La=LDPCInitialization(m,c_La);
		v_Le=LDPCInitialization(n,v_Le);
		v_La=LDPCInitialization(n,v_La);
	
	
	
		Lc=(double *)malloc(sizeof (double)*n);
		if(Lc==NULL){
			printf("Memory cannot be allocated!!\n");	exit(1);
		}
		LeSum=(double *)malloc(sizeof (double)*n);
		if(LeSum==NULL){
			printf("Memory cannot be allocated!!\n");	exit(1);
		}
		for(i=0;i<n;i++){
			Lc[i]=2.0*RecieveDATA[i]/sigma2;
			}
	
		for(i=0;i<m;i++){
				p=c_neighbor[i];
			while(p!=NULL){
				if(c_La[i]==NULL){
					c_La[i]=AddListEnd(c_La[i],p->INDEX,0.0);
					c_Le[i]=AddListEnd(c_Le[i],p->INDEX,0.0);
					q=c_La[i];
					r=c_Le[i];
				}
				else{
					q=AddListEnd(q,p->INDEX,0.0);
					r=AddListEnd(r,p->INDEX,0.0);
				}
				p=p->NEXT;
			}
	}
	
	
		for(k=0;k<Loop;k++){
			for(i=0;i<m;i++){
				q=c_La[i];
				r=c_Le[i];
				Lpro=Lanhpro(c_neighbor[i],n,Lc,c_La[i]);
				while(q!=NULL){
					test=tanh((Lc[q->INDEX]+q->VALUE)*0.5);
					divide=Lpro/(tanh(0.5*(Lc[q->INDEX]+q->VALUE)));
					if(abs(divide)>=(1.0-1e-10)){
						value=1e10*divide/abs(divide);
					}else{
					value=2.0*atanh(divide);
				}
					r->VALUE=value;
					r->INDEX=q->INDEX;
					r=r->NEXT;
					q=q->NEXT;
				}
			}
	
			ChangeList(m,c_Le,n,v_Le);
	
	
			for(i=0;i<n;i++){
				LeSum[i]=SumList(v_Le[i]);
			//	printf("%e\n",LeSum[i]+Lc[i]);
				if((LeSum[i]+Lc[i])<0){
					BE[k]++;
				//	printf("%d\n",BE[k]);
				}
			}
		//	printf("\n");
			ChangeList(m,c_La,n,v_La);
			UpdateLa(n,LeSum,v_La,v_Le);
			ChangeList(n,v_La,m,c_La);
			// FreeList(m,c_Le);
			// FreeList(n,v_La);
	
	
		}
		FreeListAll(m,c_Le);
		FreeListAll(m,c_La);
		FreeListAll(n,v_Le);
		FreeListAll(n,v_La);
	
		free(Lc);
		free(LeSum);
	
	}
	

void SumproductT(unsigned int m,unsigned int n,double RecieveDATA[n],double sigma2,unsigned int dc,unsigned int dv,struct LIST *c_neighbor[m],unsigned int Loop,unsigned int BE[Loop]){

int i,j,k;
double Tpro,divide;
struct LIST	**c_Ta,**v_Ta; //事前値
struct LIST	**c_Te,**v_Te; //外部値
double	*Tac,*Tc; //チャネル地
struct LIST *p,*q,*r,*p1,*q1;
double	value,TeSum[2];
double	test;

#ifdef debug
	FILE	*file_pt_d;
	char	file[50]="./debug/TLeLalog.txt";

	file_pt_d=fopen(file,"w");
		if(file_pt_d==NULL){
			printf("file open failed!!\n");
			exit(1);
		}

		for(i=0;i<n;i++){
					fprintf(file_pt_d,"%e	",RecieveDATA[i]);
				}
				fprintf(file_pt_d,"\n");

#endif



	c_Te=LDPCInitialization(m,c_Te);
	c_Ta=LDPCInitialization(m,c_Ta);
	v_Te=LDPCInitialization(n,v_Te);
	v_Ta=LDPCInitialization(n,v_Ta);

	Tac=(double *)malloc(sizeof (double)*dc);
	if(Tac==NULL){
		printf("Memory cannot be allocated!!\n");	exit(1);
	}
	Tc=(double *)malloc(sizeof (double)*n);
	if(Tc==NULL){
		printf("Memory cannot be allocated!!\n");	exit(1);
	}


	for(i=0;i<n;i++){
		Tc[i]=tanh(RecieveDATA[i]/sigma2);
		}
		#ifdef debug
		printf("Tc\n");
		for(i=0;i<n;i++){
				printf("%e  ",Tc[i]);
        fprintf(file_pt_d,"%e ",Tc[i]);
			}
			printf("\n");
			#endif

			for(i=0;i<m;i++){
					p=c_neighbor[i];
				while(p!=NULL){
					if(c_Ta[i]==NULL){
						c_Ta[i]=AddListEnd(c_Ta[i],p->INDEX,0.0);
						c_Te[i]=AddListEnd(c_Te[i],p->INDEX,0.0);
						q=c_Ta[i];
						r=c_Te[i];
					}
					else{
						q=AddListEnd(q,p->INDEX,0.0);
						r=AddListEnd(r,p->INDEX,0.0);
					}
					p=p->NEXT;
				}
				//printf("%d\n",i);
				//ShowList(c_Ta[i]);
		}
ChangeList(m,c_Ta,n,v_Ta);
	for(k=0;k<Loop;k++){
		for(i=0;i<m;i++){
			Tpro=Tanhpro(n,Tc,c_Ta[i],dc,Tac);
			UpdateTe(Tpro,c_Te[i],dc,Tac);
	}

		ChangeList(m,c_Te,n,v_Te);
		for(i=0;i<n;i++){
			TanhSumList(v_Te[i],TeSum);
			UpdateTa(v_Ta[i],TeSum,v_Te[i]);
			 if(TeSum[0]*(1+Tc[i])<TeSum[1]*(1-Tc[i])){
			 	BE[k]++;
			 }


		}

		ChangeList(n,v_Ta,m,c_Ta);
		// FreeList(m,c_Le);
		// FreeList(n,v_La);

	}
	FreeListAll(m,c_Te);
	FreeListAll(m,c_Ta);
	FreeListAll(n,v_Te);
	FreeListAll(n,v_Ta);

	free(Tc);
	free(Tac);
#ifdef debug
	fclose(file_pt_d);
#endif
}
void SumproductT2(unsigned int m,unsigned int n,double RecieveDATA[n],double sigma2,unsigned int dc,unsigned int dv,struct LIST *c_neighbor[m],unsigned int Loop,unsigned int BE[Loop]){

int i,j,k;
double Tpro,divide;
struct LIST	**c_Ta,**v_Ta; //事前値
struct LIST	**c_Te,**v_Te; //外部値
double	*Tac,*Tc; //チャネル地
struct LIST *p,*q,*r,*p1,*q1;


	c_Te=LDPCInitialization(m,c_Te);
	c_Ta=LDPCInitialization(m,c_Ta);
	v_Te=LDPCInitialization(n,v_Te);
	v_Ta=LDPCInitialization(n,v_Ta);

	Tac=(double *)malloc(sizeof (double)*dc);
	if(Tac==NULL){
		printf("Memory cannot be allocated!!\n");	exit(1);
	}
	Tc=(double *)malloc(sizeof (double)*n);
	if(Tc==NULL){
		printf("Memory cannot be allocated!!\n");	exit(1);
	}


	for(i=0;i<n;i++){
		Tc[i]=tanh(RecieveDATA[i]/sigma2);
		}
			for(i=0;i<m;i++){
					p=c_neighbor[i];
				while(p!=NULL){
					if(c_Ta[i]==NULL){
						c_Ta[i]=AddListEnd(c_Ta[i],p->INDEX,0.0);
						c_Te[i]=AddListEnd(c_Te[i],p->INDEX,0.0);
						q=c_Ta[i];
						r=c_Te[i];
					}
					else{
						q=AddListEnd(q,p->INDEX,0.0);
						r=AddListEnd(r,p->INDEX,0.0);
					}
					p=p->NEXT;
				}
		}
ChangeList(m,c_Ta,n,v_Ta);

	for(k=0;k<Loop;k++){
			for(i=0;i<m;i++){
			Tpro=Tanhpro(n,Tc,c_Ta[i],dc,Tac);
			UpdateTe(Tpro,c_Te[i],dc,Tac);
	}

		ChangeList(m,c_Te,n,v_Te);
		for(i=0;i<n;i++){
//			Tpro=UpdateTa2(dv,v_Ta[i],v_Te[i],Thr,file_pt);
			Tpro=UpdateTa3(v_Ta[i],v_Te[i],dv,0.0,0);
				if(tanhaplusb(Tpro,Tc[i],0,0.0)<0){
				BE[k]++;
			}

		}

		ChangeList(n,v_Ta,m,c_Ta);

	}
	FreeListAll(m,c_Te);
	FreeListAll(m,c_Ta);
	FreeListAll(n,v_Te);
	FreeListAll(n,v_Ta);

	free(Tc);
	free(Tac);
}
/****************************************************************************/
//関数名:SumproductApp
//説明: Sumproduct復号をTanh領域に移した上で全ての除算を乗算と加算で置き換えた関数
/****************************************************************************/
void SumproductApp(unsigned int m,unsigned int n,double RecieveDATA[n],double sigma2,unsigned int dc,unsigned int dv,struct LIST *c_neighbor[m],unsigned int Loop,unsigned int BE[Loop],double apparea,int THRESHOLD){

int i,j,k;
double Tpro,divide,test;
struct LIST	**c_Ta,**v_Ta; //事前値
struct LIST	**c_Te,**v_Te; //外部値
double	*Tac,*Tc; //チャネル地
struct LIST *p,*q,*r,*p1,*q1;


	c_Te=LDPCInitialization(m,c_Te);
	c_Ta=LDPCInitialization(m,c_Ta);
	v_Te=LDPCInitialization(n,v_Te);
	v_Ta=LDPCInitialization(n,v_Ta);

	Tac=(double *)malloc(sizeof (double)*dc);
	if(Tac==NULL){
		printf("Memory cannot be allocated!!\n");	exit(1);
	}
	Tc=(double *)malloc(sizeof (double)*n);
	if(Tc==NULL){
		printf("Memory cannot be allocated!!\n");	exit(1);
	}



	for(i=0;i<n;i++){
		Tc[i]=tanh(RecieveDATA[i]/sigma2);	//チャネル値の計算
			}
	for(i=0;i<m;i++){
				p=c_neighbor[i];
				while(p!=NULL){
					if(c_Ta[i]==NULL){
						c_Ta[i]=AddListEnd(c_Ta[i],p->INDEX,0.0);
						c_Te[i]=AddListEnd(c_Te[i],p->INDEX,0.0);
						q=c_Ta[i];
						r=c_Te[i];
					}
					else{
						q=AddListEnd(q,p->INDEX,0.0);
						r=AddListEnd(r,p->INDEX,0.0);
					}
				p=p->NEXT;
			}
		}
		ChangeList(m,c_Ta,n,v_Ta);
		for(k=0;k<Loop;k++){
			for(i=0;i<m;i++){
				Tanhpro2(n,Tc,c_Ta[i],dc,Tac,apparea,THRESHOLD);
				UpdateTe2(c_Te[i],dc,Tac);
	}

		ChangeList(m,c_Te,n,v_Te);
		for(i=0;i<n;i++){
			//Tpro=UpdateTa2(dv,v_Ta[i],v_Te[i],-1.0);
			Tpro=UpdateTa3(v_Ta[i],v_Te[i],dv,apparea,THRESHOLD);
			// if(test!=Tpro){
			// 	printf("test %e Tpro %e \n",test,Tpro);
			// 	printf("calc error at loop %d\n",k);
			// 	exit(1);
			// }


			if(tanhaplusb(Tpro,Tc[i],THRESHOLD,apparea)<0){
				BE[k]++;
			}

		}

		// FreeList(m,c_Le);
		ChangeList(n,v_Ta,m,c_Ta);
		// FreeList(n,v_La);

	}
		//exit(1);
	FreeListAll(m,c_Te);
	FreeListAll(m,c_Ta);
	FreeListAll(n,v_Te);
	FreeListAll(n,v_Ta);

	free(Tc);
	free(Tac);
}
