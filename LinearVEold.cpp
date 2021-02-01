//usage: ./LinearVEold <input_graph> <output_file>

#include <stdio.h>
#include <stdlib.h>
#ifdef _WIN32
	#include <chrono>
   using namespace std::chrono;
#elif _linux
	#include "rfw_timer.h"
#endif

int n, m;
int* edges;
int* G; int* firstOut;

void create_adj();
void compute_count();
int* count;

int *dfs, *idfs, *parent, *T, *firstChild, *tempChild;
int *l, *low, *high, *highp, *bcount, *up, *down, *L, *R, *M, *Mp;
int *invHighp, *firstInvHighp, *invHigh, *firstInvHigh;
int *invM, *invMfirst, *invMp, *invMpfirst;
int *backEdges, *backEdgesFirst, *backEdgesTemp;

void construct_tree();
void compute_bcount();
void sort_backEdgesHigh();
void compute_M();
void computeInvM();
void compute_high();
void compute_highp();
void sort_backEdgesLow();
void sort_ChildrenHighp();

int *ufparent, *ufrank, *representative;
void ufInit();
int find(int);
void unite(int,int,int);

int main(int n_args, char** args)
{
   FILE* fp = fopen(args[1],"r");
   fscanf(fp,"%d %d",&n,&m);
   edges = (int*)malloc(sizeof(int)*2*m);
   for(int i=0;i<m;i++)
   {
      int x,y;
      fscanf(fp,"%d %d",&x,&y);
      edges[2*i]=x-1; edges[2*i+1]=y-1;
   }
   fclose(fp);
   
   create_adj();

   #ifdef _WIN32
		std::chrono::duration<double> totalTime;
		high_resolution_clock::time_point t1 = high_resolution_clock::now();
	#elif _linux
		RFWTimer timer(true);
		double t;
	#endif
    
   compute_count();
   #ifdef _WIN32
		high_resolution_clock::time_point t2 = high_resolution_clock::now();
      duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
		 printf("%f\n",time_span.count());
	#elif _linux
		t = timer.getTime();
		printf("Total time= %g\n", t);
	#endif

   fp = fopen(args[2],"w");
   for(int i=0;i<n;i++)
   {
      fprintf(fp,"%d\n",count[i]);
   }
   #ifdef _WIN32
      fprintf(fp,"%lf\n",(double)time_span.count());
   #elif _linux
      fprintf(fp,"%lf\n", t);
   #endif
   fclose(fp);

   return 0;
}

void compute_count()
{
   count = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){count[i]=0;}

   dfs = (int*)malloc(sizeof(int)*n);
   idfs = (int*)malloc(sizeof(int)*n);
   parent = (int*)malloc(sizeof(int)*n);
   tempChild = (int*)malloc(sizeof(int)*n);
   l = (int*)malloc(sizeof(int)*n);
   low = (int*)malloc(sizeof(int)*n);
   backEdges = (int*)malloc(sizeof(int)*2*(m-n+1));
   int backEdgesCounter=0;

   //perform DFS
   for(int i=0;i<n;i++){dfs[i]=-1; l[i]=i; low[i]=i;}
   int* temp_vertex = (int*)malloc(sizeof(int)*n);
   int* temp_out = (int*)malloc(sizeof(int)*n);
   int Nr=0;
   dfs[0]=Nr; idfs[Nr++]=0; parent[0]=-1;
   temp_vertex[0]=0; temp_out[0]=firstOut[0];
   int SP=0;
   while(SP!=-1)
   {
      int v=temp_vertex[SP];
      char descend=0;
      for(int i=temp_out[SP];i<firstOut[v+1];i++)
      {
         int u=G[i];
         if(dfs[u]==-1)
         {
            dfs[u]=Nr; idfs[Nr++]=u; parent[u]=v;
            temp_vertex[SP+1]=u; temp_out[SP+1]=firstOut[u]; temp_out[SP]=i;
            descend=1; break;
         }
         if(dfs[u]<dfs[v] && u!=parent[v])
         {
            if(dfs[u]<dfs[l[v]])
            {
               l[v]=u;
               if(dfs[u]<dfs[low[v]])
               {
                  low[v]=u;
               }  
            }
            backEdges[backEdgesCounter++]=v;
            backEdges[backEdgesCounter++]=u;
         }
         else if(v==parent[u])
         {
            if(dfs[low[u]]<dfs[low[v]])
            {
               low[v]=low[u];
            }
         }
      }
      if(descend){SP++;continue;}
      SP--;
   }

   //gather parameters
   construct_tree();
   compute_bcount();
   compute_M();
   computeInvM();
   compute_high();
   compute_highp();

   //e is a back-edge
   for(int i=2;i<n;i++)
   {
      int v=idfs[i];
      count[parent[v]]+=bcount[v]==1;
   }  

   //high(u)=v
   for(int v=0;v<n;v++)
   {
      int u=firstInvHigh[v];
      int c=firstChild[v];
      int endH=firstInvHigh[v+1];
      int endT=firstChild[v+1];
      while(u!=endH)
      {    
         while(!(c==endT-1 || dfs[T[c+1]]>dfs[invHigh[u]]))
         {
            c++;
         }
         if((T[c]!=invHigh[u])&&(low[invHigh[u]]==v || dfs[invHigh[u]]<=dfs[Mp[T[c]]]))
         {
            count[v]++;
         }
         u++;
      } 
   }

   //high(u)<v
   for(int x=0;x<n;x++)
   {
      int u=invMfirst[x];
      int c=invMpfirst[x];
      int endM=invMfirst[x+1];
      int endMp=invMpfirst[x+1];
      while(u!=endM && c!=endMp)
      {
         while(c!=endMp && dfs[invMp[c]]>=dfs[invM[u]]){c++;}
         if(c==endMp){break;}
         if(dfs[high[invM[u]]]<dfs[parent[invMp[c]]])
         {
            int n_edges=0;
            int h=dfs[high[invM[u]]];
            while(c!=endMp && h<dfs[parent[invMp[c]]])
            {
               while(u!=endM && dfs[invMp[c]]<dfs[invM[u]])
               {
                  n_edges++;
                  u++;
               }
               count[parent[invMp[c]]]+=n_edges;
               c++;
            }
         }
         else
         {
            u++;
         }
      }
   }

   //M(u)=v
   sort_ChildrenHighp(); //sort the lists of children in decreasing order w.r.t. the highp of their elements

   for(int v=1;v<n;v++)
   {
      if(invMfirst[v]==invMfirst[v+1]){continue;}
      int u=invMfirst[v]+1;
      int c=firstChild[v];
      int endM=invMfirst[v+1];
      int endT=firstChild[v+1];
      int min=v;
      while(u!=endM && c!=endT)
      {
         min=highp[T[c]];
         while(u!=endM && dfs[invM[u]]>dfs[min])
         {
            count[v]++;
            u++;
         }
         min=low[T[c]];
         c++;
         while(c!=endT && dfs[highp[T[c]]]>=dfs[min])
         {
            if(dfs[low[T[c]]]<dfs[min]) 
            {
               min=low[T[c]]; 
            }
            c++;
         }
         while(u!=endM && dfs[invM[u]]>dfs[min]){u++;}
      }
      while(u!=endM)
      {
         if(dfs[invM[u]]<=dfs[min])
         {
            count[v]++;
         }
         u++;
      }
   }
   
   //M(u)>v
   for(int x=0;x<n;x++)
   {
      int c = invMpfirst[x];
      int u = invMfirst[x];
      int endMp = invMpfirst[x+1];
      int endM = invMfirst[x+1];  
      while(c!=endMp && u!=endM)
      {
         while(u!=endM && dfs[invM[u]]>=dfs[parent[invMp[c]]]){u++;}
         if(u==endM){break;}
         if(dfs[highp[invMp[c]]]<dfs[invM[u]])
         {
            int n_edges=0;
            int first = u;
            while(u!=endM && dfs[highp[invMp[c]]]<dfs[invM[u]])
            {
               n_edges++;
               u++;
            }
            int last = u-1;
            count[parent[invMp[c]]]+=n_edges;
            c++;
            while(c!=endMp && dfs[invM[last]]<dfs[parent[invMp[c]]])
            {
               while(dfs[invM[first]]>=dfs[parent[invMp[c]]])
               {
                  n_edges--;
                  first++;
               }
               count[parent[invMp[c]]]+=n_edges;
               c++;
            }
         }
         else
         {
            c++;
         }
      }
   }
}

void construct_tree()
{
   //sort the list of the children of every vertex in increasing order
   T = (int*)malloc(sizeof(int)*n);
   firstChild = (int*)malloc(sizeof(int)*(n+1));
   for(int i=0;i<=n;i++){firstChild[i]=0;}
   for(int i=1;i<n;i++){firstChild[parent[idfs[i]]+1]++;}
   for(int i=0;i<n;i++){tempChild[i]=0;}
   int tc_r=0;
   tempChild[0]=firstChild[1];
   for(int i=1;i<n;i++){firstChild[i+1]+=firstChild[i]; tempChild[i]=firstChild[i+1];}
   for(int i=1;i<n;i++)
   {
      int v=idfs[i];
      if(parent[v]==0){T[tc_r++]=v;}
      else{T[tempChild[parent[v]-1]++]=v;}
   }
}

void compute_bcount()
{
   up = (int*)malloc(sizeof(int)*n);
   down = (int*)malloc(sizeof(int)*n);
   bcount = (int*)malloc(sizeof(int)*n);
   for(int v=0;v<n;v++){up[v]=0; down[v]=0; bcount[v]=0;}

   sort_backEdgesHigh(); //sort the back-edges in increasing order w.r.t. their higher end

   for(int v=0;v<n;v++){tempChild[v]=firstChild[v];}
   
   int N=m-n+1;
   for(int i=0;i<N;i++)
   {
      int u=backEdges[2*i];
      int v=backEdges[2*i+1];
      up[u]++;
      while(tempChild[v]!=firstChild[v+1]-1 && dfs[T[tempChild[v]+1]]<dfs[u]){tempChild[v]++;}
      down[T[tempChild[v]]]++;     
   } 

   for(int i=n-1;i>1;i--)
   {
      int v=idfs[i];
      bcount[v]=up[v]-down[v];
      for(int j=firstChild[v];j<firstChild[v+1];j++)
      {
         bcount[v]+=bcount[T[j]];
      }
   }
}

void sort_backEdgesHigh()
{
   int N=m-n+1;
   backEdgesFirst = (int*)malloc(sizeof(int)*(n+1));
   for(int i=0;i<n;i++){backEdgesFirst[i]=0;}
   for(int i=0;i<N;i++){backEdgesFirst[dfs[backEdges[2*i]]+1]++;}
   tempChild[0]=backEdgesFirst[1];
   for(int i=1;i<n;i++){backEdgesFirst[i+1]+=backEdgesFirst[i]; tempChild[i]=backEdgesFirst[i+1];}
   backEdgesTemp = (int*)malloc(sizeof(int)*2*N);
   for(int i=0;i<2*N;i++){backEdgesTemp[i]=backEdges[i];}
   for(int i=0;i<N;i++)
   {
      int x=backEdgesTemp[2*i];
      int y=backEdgesTemp[2*i+1];
      int indx=dfs[x]-1;
      backEdges[2*tempChild[indx]]=x;
      backEdges[2*tempChild[indx]+1]=y;
      tempChild[indx]++;
   }
}

void compute_M()
{
   L = (int*)malloc(sizeof(int)*n);
   R = (int*)malloc(sizeof(int)*n);
   M = (int*)malloc(sizeof(int)*n);
   Mp = (int*)malloc(sizeof(int)*n);

   //calculate all M
   for(int i=n-1;i>0;i--)
   {
      int v=idfs[i];
 
      //initialize L and R to calculate M[v]
      L[v]=-1;
      for(int j=firstChild[v];j<firstChild[v+1];j++)
      {
         if(dfs[low[T[j]]]<dfs[v])
         {
            if(L[v]==-1)
            {
               L[v]=j;
            }
            R[v]=j;
         }
      }

      //compute M
      if(dfs[l[v]]<dfs[v]){M[v]=v;}
      else if(L[v]!=R[v]){M[v]=v;}
      else
      {
         int c=T[L[v]];
         int m=M[c];
         while(1)
         {
            if(dfs[l[m]]<dfs[v]){M[v]=m;break;}
            while(dfs[low[T[L[m]]]]>=dfs[v]){L[m]++;}
            while(dfs[low[T[R[m]]]]>=dfs[v]){R[m]--;}
            if(L[m]!=R[m]){M[v]=m;break;}
            c=T[L[m]];
            m=M[c];
         }
      }             
   }

   //calculate all Mp
   for(int i=n-1;i>1;i--)
   {
      int v=idfs[i];

      //initialize L and R to calculate Mp[v]
      L[v]=-1;
      for(int j=firstChild[v];j<firstChild[v+1];j++)
      {
         if(dfs[low[T[j]]]<dfs[parent[v]])
         {
            if(L[v]==-1)
            {
               L[v]=j;
            }
            R[v]=j;
         }
      }

      //compute Mp
      if(dfs[l[v]]<dfs[parent[v]]){Mp[v]=v;continue;}
      else if(L[v]!=R[v]){Mp[v]=v;}
      else
      {
         int c=T[L[v]];
         int m=Mp[c];
         while(1)
         {
            if(dfs[l[m]]<dfs[parent[v]]){Mp[v]=m;break;}
            while(dfs[low[T[L[m]]]]>=dfs[parent[v]]){L[m]++;}
            while(dfs[low[T[R[m]]]]>=dfs[parent[v]]){R[m]--;}
            if(L[m]!=R[m]){Mp[v]=m;break;}
            c=T[L[m]];
            m=Mp[c];
         }
      }
   }
}

void computeInvM()
{
    //calculate the inverse lists, and have their elements sorted in decreasing order
    invM = new int[n];
    invMfirst = new int[n+1];
    for(int i=0;i<=n;i++){invMfirst[i]=0;}
    for(int i=1;i<n;i++)
    {
       int v=idfs[i];
       invMfirst[M[v]+1]++;
    }
    for(int i=1;i<=n;i++){invMfirst[i]+=invMfirst[i-1];}
    int* invMnext = new int[n+1];
    for(int i=0;i<=n;i++){invMnext[i]=invMfirst[i];}
    for(int i=n-1;i>0;i--)
    {
       int v=idfs[i];
       invM[invMnext[M[v]]++]=v;
    }   

    invMp = new int[n];
    invMpfirst = new int[n+1];
    for(int i=0;i<=n;i++){invMpfirst[i]=0;}
    for(int i=2;i<n;i++)
    {
       int v=idfs[i];
       invMpfirst[Mp[v]+1]++;
    }
    for(int i=1;i<=n;i++){invMpfirst[i]+=invMpfirst[i-1];}
    for(int i=0;i<=n;i++){invMnext[i]=invMpfirst[i];}
    for(int i=n-1;i>1;i--)
    {
       int v=idfs[i];
       invMp[invMnext[Mp[v]]++]=v;
    } 
}

void compute_high()
{
   sort_backEdgesLow(); //sort the back-edges in decreasing order w.r.t. their lower end
   high = (int*)malloc(sizeof(int)*n);
   ufInit();
   for(int i=m-n;i>=0;i--)
   {
      int u=backEdges[2*i];
      int v=backEdges[2*i+1];
      int x = representative[find(u)];
      while (x!=v) 
      {
         high[x] = v;    
         int next = representative[find(parent[x])];
         unite(x,parent[x],next);
         x = next;
      }
   }

   //sort the elements of the inverse lists high^-1 in increasing order
   invHigh = (int*)malloc(sizeof(int)*n);
   firstInvHigh = (int*)malloc(sizeof(int)*(n+1));
   for(int i=0;i<=n;i++){firstInvHigh[i]=0;}
   for(int i=1;i<n;i++){firstInvHigh[high[idfs[i]]+1]++;}
   int temp0=0;
   tempChild[0]=firstInvHigh[1];
   for(int i=1;i<n;i++){firstInvHigh[i+1]+=firstInvHigh[i]; tempChild[i]=firstInvHigh[i+1];}
   for(int i=1;i<n;i++)
   {
      int v=idfs[i];
      if(high[v]==0){invHigh[temp0++]=v;}
      else{invHigh[tempChild[high[v]-1]++]=v;}
   }
}

void sort_backEdgesLow()
{
   int N=m-n+1;
   for(int i=0;i<n;i++){backEdgesFirst[i]=0;}
   for(int i=0;i<N;i++){backEdgesFirst[dfs[backEdges[2*i+1]]+1]++;}
   int temp0=0;
   tempChild[0]=backEdgesFirst[1];
   for(int i=1;i<n;i++){backEdgesFirst[i+1]+=backEdgesFirst[i]; tempChild[i]=backEdgesFirst[i+1];}
   for(int i=0;i<2*N;i++){backEdgesTemp[i]=backEdges[i];}
   for(int i=0;i<N;i++)
   {
      int x=backEdgesTemp[2*i];
      int y=backEdgesTemp[2*i+1];
      if(dfs[y]==0)
      {
         backEdges[2*temp0]=x;
         backEdges[2*temp0+1]=y;  
         temp0++;   
      }
      else
      {
         int indx=dfs[y]-1;
         backEdges[2*tempChild[indx]]=x;
         backEdges[2*tempChild[indx]+1]=y;
         tempChild[indx]++;
      }
   }
}

void compute_highp()
{
   highp = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){ufparent[i]=i; ufrank[i]=0; representative[i]=i;}
   for(int i=m-n;i>=0;i--)
   {
      int u=backEdges[2*i];
      int v=backEdges[2*i+1];
      int x = representative[find(u)];
      while (parent[x]!=v) 
      {
         highp[x] = v;    
         int next = representative[find(parent[x])];
         unite(x,parent[x],next);
         x = next;
      }
   }

   //sort the elements of the inverse lists highp^-1 in increasing order
   invHighp = (int*)malloc(sizeof(int)*n);
   firstInvHighp = (int*)malloc(sizeof(int)*(n+1));
   for(int i=0;i<=n;i++){firstInvHighp[i]=0;}
   for(int i=2;i<n;i++){firstInvHighp[highp[idfs[i]]+1]++;}
   int temp0=0;
   tempChild[0]=firstInvHighp[1];
   for(int i=1;i<n;i++){firstInvHighp[i+1]+=firstInvHighp[i]; tempChild[i]=firstInvHighp[i+1];}
   for(int i=2;i<n;i++)
   {
      int v=idfs[i];
      if(highp[v]==0){invHighp[temp0++]=v;}
      else{invHighp[tempChild[highp[v]-1]++]=v;}
   }
}

void sort_ChildrenHighp()
{
   for(int i=0;i<n;i++){tempChild[i]=firstChild[i];}
   for(int i=n-1;i>=0;i--)
   {
      int v=idfs[i];
      for(int i=firstInvHighp[v];i<firstInvHighp[v+1];i++)
      {
         int c=invHighp[i];
         T[tempChild[parent[c]]++]=c;
      }
   }
}

void ufInit()
{
   ufparent = (int*)malloc(sizeof(int)*n);
   ufrank = (int*)malloc(sizeof(int)*n);
   representative = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){ufparent[i]=i; ufrank[i]=0; representative[i]=i;}
}

int find(int x)
{
   int r=x;
   while(ufparent[r]!=r){r=ufparent[r];}
   while(ufparent[x]!=x){int next=ufparent[x]; ufparent[x]=r; x=next;}
   return r;
}

void unite(int x, int y, int w)
{
   int r1=find(x);
   int r2=find(y);
   if(r1==r2){representative[r1]=w;return;}
   int cmp=ufrank[r1]-ufrank[r2];
   if(cmp<0)
   {
      ufparent[r1]=r2;
      representative[r2]=w;
   }
   else if(cmp>0)
   {
      ufparent[r2]=r1;
      representative[r1]=w;
   }
   else
   {
      ufparent[r1]=r2;
      ufrank[r2]++;
      representative[r2]=w;
   }
}

void create_adj()
{
   G = (int*)malloc(sizeof(int)*4*m);
   firstOut = (int*)malloc(sizeof(int)*(n+1));
   for(int i=0;i<=n;i++){firstOut[i]=0;}
   for(int i=0;i<m;i++){firstOut[edges[2*i]+1]++; firstOut[edges[2*i+1]+1]++;}
   int* nextOut = (int*)malloc(sizeof(int)*(n+1));
   nextOut[0]=0;
   for(int i=1;i<=n;i++){firstOut[i]+=firstOut[i-1]; nextOut[i]=firstOut[i];}
   for(int i=0;i<m;i++)
   {
      int x=edges[2*i]; int y=edges[2*i+1];
      G[nextOut[x]++]=y;
      G[nextOut[y]++]=x;
   }
}

