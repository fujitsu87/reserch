#ifndef LDPC_STRUCTURE
#define LDPC_STRUCTURE 1

#include  <stdlib.h>
#include  <stdio.h>
#include "./matrix.c"

struct LIST{
  unsigned int INDEX;
  double VALUE;
  struct LIST *NEXT;
};

struct LIST** LDPCInitialization(unsigned int m,struct LIST **neighbor){
  int i;

  neighbor=(struct LIST**)malloc (sizeof(struct LIST)*m);
  if(neighbor==NULL){
    printf("Memory cannot be allocated!!\n");
    exit(1);
  }

  for(i=0;i<m;i++){
    neighbor[i]=NULL;
  }
  return neighbor;
}

struct LIST *AddListEnd(struct LIST *p,unsigned int index,double value){

  struct LIST *q,*tn;

  if((q=(struct LIST *)malloc (sizeof(struct LIST)))==NULL){
      printf("malloc error");
      exit (EXIT_FAILURE);
  }
  q->INDEX=index;   //INDEX値を代入
  q->VALUE=value;   //VALUE値を代入
  q->NEXT=NULL;     //最後に付け足すから次のポインタはNULL

  tn=p;

  if(tn==NULL){   //付け足す前のポインタがNULLなら作ったものを代入
    return q;
  }
  else {
    if(tn->NEXT!=NULL){
      printf("AddlistEnd Error!!\n\n");
      exit(1);
    }
    else{
      tn->NEXT=q;
    }
    return q;
  }
}

void FreeList(unsigned int m,struct LIST *neighbor[m]){
 struct LIST *p,*tn;
 int i;


 for(i=0;i<m;i++){
   tn=neighbor[i];
    while(tn!=NULL){   //各配列の要素を解放
      p=tn->NEXT;
      free(tn);
      tn=p;
    }
    neighbor[i]=NULL;
  }
}
void FreeListAll(unsigned int m,struct LIST **neighbor){

   FreeList(m,neighbor);
   free(neighbor);
}
void ChangeList(unsigned int m,struct LIST *neighbor_src[m],unsigned int n,struct LIST *neighbor_des[n]){
  int i,j;
  struct LIST *p,**tn;

  tn=(struct LIST**)malloc (sizeof(struct LIST)*n);
  if(tn==NULL){
    printf("Memory cannot be allocated!!\n");
    exit(1);
  }

  FreeList(n,neighbor_des);


  for(i=0;i<m;i++){
    p=neighbor_src[i];
    while(p!=NULL){
      if(neighbor_des[p->INDEX]==NULL){
          neighbor_des[p->INDEX]=AddListEnd(NULL,i,p->VALUE);
          tn[p->INDEX]=neighbor_des[p->INDEX];
        }
      else {
          tn[p->INDEX]=AddListEnd(tn[p->INDEX],i,p->VALUE);
      }
      p=p->NEXT;
    }
  }
}


void PrintMatrixV(unsigned int m,unsigned int n,struct LIST *v_neighbor[n]){
    int i,j;
    int **ldm;
  struct LIST *p;

    ldm=IMatrix(m,n,ldm);

    for(i=0;i<n;i++){
      p=v_neighbor[i];
      while(p!=NULL){
            ldm[p->INDEX][i]=1;
            p=p->NEXT;
          }
    }
    for(i=0;i<m;i++){
      for(j=0;j<n;j++){
        printf("%d",ldm[i][j]);
        }
        printf("\n");
      }
      IMatrixFree(ldm);
}

void PrintMatrixC(unsigned int m,unsigned int n,struct LIST *c_neighbor[m]){
    int i,j;
    int **ldm;
  struct LIST *p;

      ldm=IMatrix(m,n,ldm);

      for(i=0;i<m;i++){
        p=c_neighbor[i];
        while(p!=NULL){
            ldm[i][p->INDEX]=1;
            p=p->NEXT;
        }
      }
    for(i=0;i<m;i++){
      for(j=0;j<n;j++){
        printf("%d",ldm[i][j]);
        }
        printf("\n");
      }

      IMatrixFree(ldm);
}

#endif
