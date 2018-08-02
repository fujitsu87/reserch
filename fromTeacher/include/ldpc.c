#ifndef LDPC
#define LDPC	1

#include "./random_number.c"
#include "./ldpc_structure.c"

void GetSparseVector(unsigned int n,unsigned char d,unsigned int index[d]){

  int i,j,k;

  for(i=0;i<d;i++){
    index[i]=(unsigned int)gsl_rng_uniform_int(RAN,n);
    k=0;
    do{
      if(k>0){k=0;}
      for(j=0;j<i;j++){
          if(index[i]==index[j]){
            k++;
          }
        }
        if(k>0){
          index[i]=(unsigned int)gsl_rng_uniform_int(RAN,n);
        }
      }
    while(k!=0);
    }
}
void MoveColumns(unsigned int n,struct LIST *v_neighbor[n],unsigned char weight[n],unsigned char dv){
  int i,j,k,l;
  struct LIST *p;
  unsigned char  tw;

  k=0;
  i=0;
  l=n-1;
while(k==0){
    if(weight[i]>=dv){
      p=v_neighbor[i];
      tw=weight[i];
      v_neighbor[i]=v_neighbor[l];
      weight[i]=weight[l];
      v_neighbor[l]=p;
      weight[l]=tw;
      l--;

    }
    else{
      i++;
      }
    if(i==l){k++;}
}
}

unsigned int CheckIndex (unsigned int j,unsigned int index[],unsigned char dv,unsigned char weight[]){
  unsigned int k,l;
  k=0;
  for(l=0;l<j;l++){
      if(index[j]==index[l]  ){   //index[j]は他の値と異なっていて且つ列重みがdvになったところには追加しない
        k++;
      }
    }
    if(weight[index[j]]>=dv){
      k++;
    }
    return k;
}


void MakeLDPC(unsigned int n,unsigned int m,unsigned char dv,unsigned char dc, struct LIST *v_neighbor[n]){
  unsigned int i,j,k,l,dvless=0;
  unsigned int *index;
  unsigned char *weight;
  unsigned int dover;
  unsigned int dn=0;
  struct LIST **p;
  int threshold;

  if(n*dv!=m*dc){
    printf("Failed to make LDPC\n");
    exit(1);
  }

  p=LDPCInitialization(n,p);
  index=(unsigned int *)calloc(dc,sizeof(unsigned int));
  if(index==NULL){
    exit(1);
  }
  weight=(unsigned char *)calloc(n,sizeof(unsigned char));
  if(weight==NULL){
    exit(1);
  }

  for(i=0;i<m;i++){
    if(i>=(m-dv)){       //最後のdv行分の処理
      dvless=0;
      for(j=0;j<n;j++){           //必要な列重みが残り行数以上であればそこに1を立てる
        if(dv-weight[j]>=(m-i)){
          index[dvless]=j;
          dvless++;
        }
      }
      for(j=dvless;j<dc;j++){   //列重みの辻褄合わせのところを覗いてランダムにindexを埋める．
        index[j]=(unsigned int)gsl_rng_uniform_int(RAN,n);
        do{
          k=CheckIndex(j,index,dv,weight);
            if(k>0){
              index[j]=(unsigned int)gsl_rng_uniform_int(RAN,n);

            }
          }
        while(k!=0);
        }
          }
    else{
    for(j=0;j<dc;j++){      //ランダムに行重み分だけindexを埋める．
      index[j]=(unsigned int)gsl_rng_uniform_int(RAN,n);
      do{
          k=CheckIndex(j,index,dv,weight);
          if(k>0){
            index[j]=(unsigned int)gsl_rng_uniform_int(RAN,n);
          }
        }
      while(k!=0);
      }
    }
      for(j=0;j<dc;j++){        //得られたindexを保存
        if(v_neighbor[index[j]]==NULL){
          v_neighbor[index[j]]=AddListEnd(NULL,i,1.0);
          p[index[j]]=v_neighbor[index[j]];
          weight[index[j]]=1;
        }
        else{
          p[index[j]]=AddListEnd(p[index[j]],i,1.0);
          weight[index[j]]++;
        }
      }
      dover=0;
      for(j=0;j<dc;j++){
        if(weight[index[j]]>=dv){
          dover++;        //列重みが基準に達している列の数
        }
      }
    }

  free(index);
  free(weight);
}


#endif