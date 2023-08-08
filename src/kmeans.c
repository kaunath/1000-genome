#include<math.h>

void calc_centroid(double *arr,double *arc,int *cluster_no,double *parc,int points,int centers,int dim)//,int points)
{
	double x1[centers][dim]; //To store data points
	
	double dist[points][centers]; // dist[i][j]: Distance from ith data point to jth cluster
	
	int i,j;
	int a[centers];// a[j]: a to store number of points in the jth cluster
	 // cluster Number the point belongs to
	
	
	for(int i=0;i<centers;i++){
        for(int j=0;j<dim;j++){
            *(parc+dim*i+j)= *(arc+dim*i+j);
        }
    }
	//To compute distances to each point from centres
	for(i =0;i<points; i++){
		for(j=0 ;j<centers;j++){
		 // Calculating distance
            float sum=0;
            for(int k=0;k<dim;k++){
            sum+=pow(*(arr+i*dim+k) - *(arc+j*dim+k),2);
            }
            dist[i][j] =sqrt(sum); ;
	    }
	}

	// to get the cluster number
	for(i=0;i<points;i++) {
	   int j;
	   double min;
	   int index=0;
	   min = dist[i][0];
	   for(j=0;j<centers;j++){
		  if(min > dist[i][j]){
			min = dist[i][j];
			index = j;
		  }
	   }
	
	//return index;
	   cluster_no[i] = index;
	}
	
	for( i=0;i<centers;i++){
	    a[i] = 0;
		for(int j=0;j<dim;j++){
		    x1[i][j]=0;
		}
	}
	for(j=0;j<centers;j++){
	    
		for(i=0;i<points;i++){
		
			if(cluster_no[i]== j){
			    
			 for(int k=0;k<dim;k++){
			     x1[j][k]+=*(arr+dim*i+k);
			    }
			   a[j]++;
	         }
		}
    }
	for(j=0;j<centers;j++){
	   for(int k=0;k<dim;k++){
	       if(a[j]!=0){
	       *(arc+dim*j+k)=x1[j][k]/(float)a[j];
	       }
	   }
	}
}
