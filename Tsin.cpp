//usage: ./Tsin <input_graph> <output_file>

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

int n, m;
int* edges;
int* G; int* firstOut;

void create_adj();
void find_cut_edges();
int* cut_edges;
int n_cut_edges;

int* dfs; int* nd; int* lowpt; int* low; int* lowpt2; int* low2; int* tolow; int* parent;
int Nr=1;

void Find_cut_pairs(int);
int* tree_edge_is_cutEdge; int* back_edge_is_cutEdge;
void insert_edge(int,int,char);

int* STACK; int* firstIndx; int* prevIndx;
int SP=0;
void push(int,int,int,int,int);
void pop(int);


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

   struct timeval begin, end;
   gettimeofday(&begin, 0);  

   find_cut_edges();

   //count cut-pairs
   int n_cut_pairs=0;
   for(int i=0;i<n;i++)
   {
      int c=tree_edge_is_cutEdge[i];
      n_cut_pairs+=(c*(c+1))/2;
      c=back_edge_is_cutEdge[i];
      n_cut_pairs+=(c*(c+1))/2;
   }

   gettimeofday(&end, 0);
	long seconds = end.tv_sec - begin.tv_sec;
	long microseconds = end.tv_usec - begin.tv_usec;
	double elapsed = seconds + microseconds*1e-6;
	printf("Total time= %g\n", elapsed);

   fp = fopen(args[2],"w");
   fprintf(fp,"%d\n",n_cut_edges);
   for(int i=0;i<n_cut_edges;i++)
   {
      fprintf(fp,"%d %d\n",cut_edges[2*i]+1,cut_edges[2*i+1]+1);
   }
   fprintf(fp,"%d\n",n_cut_pairs);
   fprintf(fp, "%lf\n", elapsed);
   fclose(fp);
   
   return 0;
}

void push(int v, int x, int y, int p, int q)
{
   prevIndx[SP]=firstIndx[v];
   firstIndx[v]=SP;
   int sp=SP*4;
   STACK[sp]=x;
   STACK[sp+1]=y;
   STACK[sp+2]=p;
   STACK[sp+3]=q;
   SP++;
}

void pop(int v)
{
   firstIndx[v]=prevIndx[firstIndx[v]];
}

void find_cut_edges()
{
   cut_edges = (int*)malloc(sizeof(int)*(4*n-4));
   n_cut_edges = 0;

   tree_edge_is_cutEdge = (int*)malloc(sizeof(int)*n);
   back_edge_is_cutEdge = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++)
   {
      tree_edge_is_cutEdge[i]=0;
      back_edge_is_cutEdge[i]=0;
   }
 
   dfs = (int*)malloc(sizeof(int)*n); 
   for(int i=0;i<n;i++){dfs[i]=-1;}
   parent = (int*)malloc(sizeof(int)*n); 
   nd = (int*)malloc(sizeof(int)*n);   
   lowpt = (int*)malloc(sizeof(int)*n);   
   low = (int*)malloc(sizeof(int)*n);   
   lowpt2 = (int*)malloc(sizeof(int)*n);   
   low2 = (int*)malloc(sizeof(int)*n);   
   tolow = (int*)malloc(sizeof(int)*n);  

   STACK = (int*)malloc(sizeof(int)*m*4);
   firstIndx = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){firstIndx[i]=-1;}
   prevIndx = (int*)malloc(sizeof(int)*m);

   parent[0]=-1; 
   Find_cut_pairs(0);  
}

void Find_cut_pairs(int v)
{
   dfs[v]=Nr++;
   nd[v]=1;
   lowpt[v]=dfs[v]; low[v]=v;
   lowpt2[v]=dfs[v]; low2[v]=v;

   int* temp_vertex = (int*)malloc(sizeof(int)*n);
   int* temp_out = (int*)malloc(sizeof(int)*n);
   char* firstTime = (char*)malloc(sizeof(char)*n);
   for(int i=0;i<n;i++){firstTime[i]=0;}
   int stack_pointer=0;
   
   temp_vertex[0]=v; temp_out[0]=firstOut[v];

while(stack_pointer!=-1)
{
   v=temp_vertex[stack_pointer];
   char descend=0;

   //1
   for(int i=temp_out[stack_pointer];i<firstOut[v+1];i++)
   {
      int w=G[i];
      //1.1
      if(dfs[w]==-1)
      {
         parent[w]=v;
         dfs[w]=Nr++;
         nd[w]=1;
         lowpt[w]=dfs[w]; low[w]=w;
         lowpt2[w]=dfs[w]; low2[w]=w;
         temp_vertex[stack_pointer+1]=w; temp_out[stack_pointer+1]=firstOut[w];
         temp_out[stack_pointer]=i; firstTime[stack_pointer]=1;
         descend=1; break;
      }
      if(firstTime[stack_pointer])
      {
         firstTime[stack_pointer]=0;
         //1.1.1
         if(firstIndx[w]!=-1)
         {
            int sp=firstIndx[w]*4;
            int x=STACK[sp];
            int y=STACK[sp+1];
            int p=STACK[sp+2];
            int q=STACK[sp+3];
            if(w==q)
            {
               pop(w);
               //printf("(%d,%d) (%d,%d)\n",v+1,w+1,x+1,y+1);
               insert_edge(x,y,1);
               insert_edge(v,w,0);
               if(v!=p)
               {
                  push(w,x,y,p,v);
               }
            }
         }
         nd[v]=nd[v]+nd[w];
         //1.1.2
         if(lowpt[w]<lowpt[v])
         {
            lowpt2[v]=lowpt[v]; lowpt[v]=lowpt[w];
            low2[v]=low[v]; low[v]=low[w];
            firstIndx[v]=firstIndx[w];
            tolow[v]=w; 
         }
         //1.1.3
         else if(lowpt[w]<lowpt2[v])
         {
            lowpt2[v]=lowpt[w]; low2[v]=low[w];
            firstIndx[w]=-1;
         }
      }
      //1.2
      else 
      {
         if(dfs[w]<dfs[v] && w!=parent[v])
         {
            //1.2.1
            if(dfs[w]<=lowpt[v])
            {
               lowpt2[v]=lowpt[v]; lowpt[v]=dfs[w];
               low2[v]=low[v]; low[v]=w;
               firstIndx[v]=-1;
               tolow[v]=w;
            }
            else if(dfs[w]<lowpt2[v])
            {
               lowpt2[v]=dfs[w]; low2[v]=w;
            }
         }
      }
   }

   if(descend){stack_pointer++; continue;}   

   //2
   //2.1
   if(firstIndx[v]==-1)
   {   
      //2.1.1
      if(lowpt2[v]>lowpt[v])
      {
         push(v,v,tolow[v],low[v],low2[v]);
      }
   }
   //2.2
   else
   {
      int sp=firstIndx[v]*4;
      int x=STACK[sp];
      int y=STACK[sp+1];
      int p=STACK[sp+2];
      int q=STACK[sp+3];
      //2.2.1
      if(lowpt2[v]>dfs[q])
      {
         push(v,v,tolow[v],q,low2[v]);
      }
      //2.2.2
      else
      {
         while(firstIndx[v]!=-1 && lowpt2[v]<=dfs[p])
         {
            pop(v);
            if(firstIndx[v]!=-1)
            {
               sp=firstIndx[v]*4;
               x=STACK[sp];
               y=STACK[sp+1];
               p=STACK[sp+2];
               q=STACK[sp+3];            
            }         
         }
         if(firstIndx[v]!=-1 && lowpt2[v]<dfs[q])
         {
            push(v,x,y,p,low2[v]);
         }
      }
   } 

   //3
   for(int i=firstOut[v];i<firstOut[v+1]&&firstIndx[v]!=-1;i++)
   {
      int u=G[i];
      if(!(dfs[v]<dfs[u] && v!=parent[u])){continue;}
      int sp=firstIndx[v]*4;
      int x=STACK[sp];
      int y=STACK[sp+1];
      while(firstIndx[v]!=-1 && x==parent[y] && dfs[y]<=dfs[u] && dfs[u]<=dfs[y]+nd[y]-1)
      {
         pop(v);
         if(firstIndx[v]!=-1)
         {
            sp=firstIndx[v]*4;
            x=STACK[sp];
            y=STACK[sp+1];
         }
      }
   }

   stack_pointer--;
}

}

void insert_edge(int a, int b, char generator)
{
   if(generator)
   {
      if(a==parent[b])
      {
         if(!tree_edge_is_cutEdge[b])
         {

            cut_edges[2*n_cut_edges]=a;
            cut_edges[2*n_cut_edges+1]=b;
            n_cut_edges++;           
         }
         tree_edge_is_cutEdge[b]++;
      }
      else
      {
         if(!back_edge_is_cutEdge[a])
         {
            cut_edges[2*n_cut_edges]=a;
            cut_edges[2*n_cut_edges+1]=b;
            n_cut_edges++;  
         }
         back_edge_is_cutEdge[a]++;
      }
   }
   else
   {
      cut_edges[2*n_cut_edges]=a;
      cut_edges[2*n_cut_edges+1]=b;
      n_cut_edges++;
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










