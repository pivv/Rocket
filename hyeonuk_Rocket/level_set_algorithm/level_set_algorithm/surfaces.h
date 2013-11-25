#include "operators.h"

struct surface
{
	int dimension;
	int point_num;
	double ** point;		// point[point_num][dimension]
	int edge_num;
	int ** edge;		// edge[edge_num][2]
	int surf_num;
	int ** surf;		// surf[surf_num][3]
	surface(int dim, int pointnum, int edgesurfnum);
	~surface();
};


surface::surface(int dim, int pointnum, int edgesurfnum)
{
	int i;
	dimension = dim;
	point_num = pointnum;
	point = new double *[point_num];
	for(i=0; i<point_num; i++) point[i] = new double[dimension];

	if (dimension == 2)
	{
		edge_num = edgesurfnum;
		edge = new int *[edge_num];
		for(i=0; i<edge_num; i++) edge[i] = new int[2];
	}
	if (dimension == 3)
	{
		surf_num = edgesurfnum;
		surf = new int *[surf_num];
		for(i=0; i<surf_num; i++) surf[i] = new int[3];
	}
}

surface::~surface()
{
	int i;

	if(point!=NULL)
	{
		for(i=0; i<point_num; i++) delete [] point[i];
		delete [] point;
	}

	if (dimension == 2)
	{
		if(edge!=NULL)
		{
			for(i=0; i<edge_num; i++) delete [] edge[i];
			delete [] edge;
		}
	}
	if (dimension == 3)
	{
		if(surf!=NULL)
		{
			for(i=0; i<surf_num; i++) delete [] surf[i];
			delete [] surf;
		}
	}
}

double unsigned_distance_face_point(struct surface* surf, int i0, double *v)
{
	if (surf->dimension==2)
	{
		int i;
		int imin = 0;
		double tempdist;
		double pointdist;

		for(i=0; i<2; i++)
		{
			tempdist = sqrt(pow((v[0]-surf->point[surf->edge[i0][i]][0]),2) + pow((v[1]-surf->point[surf->edge[i0][i]][1]),2));
			if(i==0) pointdist = tempdist;
			else if(tempdist < pointdist)
			{
				pointdist = tempdist; 
				imin = i;
			}
		}

		double *a1;
		double *a2;

		if(imin == 0)
		{
			a1 = surf->point[surf->edge[i0][0]];
			a2 = surf->point[surf->edge[i0][1]];
		}
		else
		{
			a1 = surf->point[surf->edge[i0][1]];
			a2 = surf->point[surf->edge[i0][0]];
		}

		if(vec_inner_product(2,v,a1,a1,a2) > 0) tempdist = pointdist;
		else tempdist = distance_face_point(2,v,a1,a2,NULL);

		return tempdist;
	}

	if (surf->dimension==3)
	{
		int i;
		int imin = 0;
		double tempdist;
		double pointdist;

		for(i=0; i<3; i++)
		{
			tempdist = sqrt(pow((v[0]-surf->point[surf->surf[i0][i]][0]),2) + pow((v[1]-surf->point[surf->surf[i0][i]][1]),2) + pow((v[2]-surf->point[surf->surf[i0][i]][2]),2));
			if(i==0) pointdist = tempdist;
			else if(tempdist < pointdist)
			{
				pointdist = tempdist; 
				imin = i;
			}
		}

		double *a1;
		double *a2;
		double *a3;

		if(imin == 0)
		{
			a1 = surf->point[surf->surf[i0][0]];
			a2 = surf->point[surf->surf[i0][1]];
			a3 = surf->point[surf->surf[i0][2]];
		}
		else if(imin == 1)
		{
			a1 = surf->point[surf->surf[i0][1]];
			a2 = surf->point[surf->surf[i0][2]];
			a3 = surf->point[surf->surf[i0][0]];
		}
		else
		{
			a1 = surf->point[surf->surf[i0][2]];
			a2 = surf->point[surf->surf[i0][0]];
			a3 = surf->point[surf->surf[i0][1]];
		}
			
		double *n = vec_curl(a1, a2, a1, a3);
		double zerovec[3] = {0,};

		double *c1 = vec_curl(v,a1,a1,a2);
		double *c2 = vec_curl(v,a1,a1,a3);

		double rot1 = vec_inner_product(3,zerovec,c1,zerovec,n);
		double rot2 = vec_inner_product(3,zerovec,c2,zerovec,n);

		double inner1 = vec_inner_product(3,v,a1,a1,a2);
		double inner2 = vec_inner_product(3,v,a1,a1,a3);

		if(inner1>0 && inner2>0) tempdist = pointdist;
		else if(rot1>0 && rot2<0) tempdist = distance_face_point(3,v,a1,a2,a3);
		else if(rot1<0) tempdist = distance_face_point(3,v,a1,a2,NULL);
		else tempdist = distance_face_point(3,v,a1,a3,NULL);
		
		delete [] n;
		delete [] c1;
		delete [] c2;

		return tempdist;
	}
	return 0;
}


int line_face_intersecting(struct surface *surf, int i0, double *l, double *v)
{
	if(surf->dimension == 2)
	{
		double *a = surf->point[surf->edge[i0][0]];
		double *b = surf->point[surf->edge[i0][1]];
		double l1[2] = {-l[1], l[0]};
		double zerovec[2] = {0,};

		double a1 = vec_inner_product(2,a,v,zerovec,l1);
		double b1 = vec_inner_product(2,a,b,zerovec,l1);

		double d1 = vec_inner_product(2,a,v,zerovec,l);
		double d2 = vec_inner_product(2,a,b,zerovec,l);

		if(abs(b1)<MINERROR) return -1;

		double s = a1/b1;
		double t = s*d2-d1;

		if(t>0 && s>0 && s<1) return 1;
		else if(s==0 || s==1) return -1;
		else return 0;
	}

	if(surf->dimension == 3)
	{
		double *a = surf->point[surf->surf[i0][0]];
		double *b = surf->point[surf->surf[i0][1]];
		double *c = surf->point[surf->surf[i0][2]];

		double m[3] = {(double) rand()/RAND_MAX, (double) rand()/RAND_MAX, (double) rand()/RAND_MAX};
		
		double zerovec[3] = {0,};
		double *l1 = vec_curl(zerovec,l,zerovec,m);
		vec_normalization(3,l1);

		double *l2 = vec_curl(zerovec,l,zerovec,l1);
		
		double a1 = vec_inner_product(3,a,v,zerovec,l1);
		double a2 = vec_inner_product(3,a,v,zerovec,l2);
		double b1 = vec_inner_product(3,a,b,zerovec,l1);
		double b2 = vec_inner_product(3,a,b,zerovec,l2);
		double c1 = vec_inner_product(3,a,c,zerovec,l1);
		double c2 = vec_inner_product(3,a,c,zerovec,l2);

		delete [] l1;
		delete [] l2;

		double d1 = vec_inner_product(3,a,v,zerovec,l);
		double d2 = vec_inner_product(3,a,b,zerovec,l);
		double d3 = vec_inner_product(3,a,c,zerovec,l);

		double bc = b1*c2-b2*c1;
		if(abs(bc)<MINERROR) return -1;

		double s = (a1*c2-a2*c1)/bc;
		double r = (a1*b2-a2*b1)/(-bc);

		double t = s*d2+r*d3-d1;

		if(t>0 && s>0 && r>0 && s+r<1) return 1;
		else if(s==0 || r==0 || s+r==0) return -1;
		else return 0;
	}

	return 0;
}

double line_surface_intersecting(struct surface *surf, double *v)
{
	if(surf->dimension == 2)
	{
		int i;
		int num;

		double sign = -1;

		while(sign==-1)
		{
			num = 0;
			double l[2] = {(double) rand()/RAND_MAX, (double) rand()/RAND_MAX};
			vec_normalization(2,l);
			for(i=0; i<surf->edge_num; i++)
			{
				sign = line_face_intersecting(surf,i,l,v);
				if(sign==-1) break;
				if(sign==1) num++;
			}
		}

		if(num%2 == 0) return 1; // outer
		else return -1; // inner
	}

	if(surf->dimension == 3)
	{
		int i;
		int num;

		double sign = -1;
		
		while(sign==-1)
		{
			num = 0;
			double l[3] = {(double) rand()/RAND_MAX, (double) rand()/RAND_MAX, (double) rand()/RAND_MAX};
			vec_normalization(3,l);
			for(i=0; i<surf->surf_num; i++)
			{
				sign = line_face_intersecting(surf,i,l,v);
				if(sign==-1) break;
				if(sign==1) num++;
			}
		}
		
		if(num%2 == 0) return 1; // outer
		else return -1; // inner
	}

	return 0;
}

double temp_distance_surf_point(double x, double y, double z)
{
	if(x>y)
	{
		double r = sqrt((x-0.75)*(x-0.75) + (y-0.25)*(y-0.25) + (z-0.5)*(z-0.5));
		return r - 0.1;
	}
	else
	{
		double r = sqrt((x-0.25)*(x-0.25) + (y-0.75)*(y-0.75) + (z-0.5)*(z-0.5));
		return r - 0.1;
	}
}

double distance_surf_point(struct surface* surf, double x, double y, double z)
{
	int i;
	double dist;
	double tempdist;

	if (surf->dimension==2) // counter-clockwise
	{
		double vec[2] = {x,y};
		for(i=0; i<surf->edge_num; i++)
		{
			tempdist = unsigned_distance_face_point(surf, i, vec);
			if(i==0) dist = tempdist;
			else if(tempdist <= dist) dist = tempdist;
		}
		if(dist==0) return dist;
		return line_surface_intersecting(surf,vec) * dist;
	}
	if (surf->dimension==3) // outer-pointing orientation
	{
		//return temp_distance_surf_point(x,y,z);
		double vec[3] = {x,y,z};
		for(i=0; i<surf->surf_num; i++)
		{
			tempdist = unsigned_distance_face_point(surf, i, vec);
			if(i==0) dist = tempdist;
			else if(tempdist <= dist) dist = tempdist;
		}
		if(dist==0) return dist;
		return line_surface_intersecting(surf,vec) * dist;
	}

	return 0;
}

int add_point(struct surface* surf, int point_num, double x, double y, double z)
{
	int i;
	if(surf->point_num < point_num) return 0;
	
	if (surf->dimension==2)
	{
		for(i=0; i<point_num; i++)
		{
			if(surf->point[i][0]==x && surf->point[i][1]==y) break;
		}
		if(i==point_num)
		{
			surf->point[i][0] = x;
			surf->point[i][1] = y;
		}
	}
	if (surf->dimension==3)
	{
		for(i=0; i<point_num; i++)
		{
			if(surf->point[i][0]==x && surf->point[i][1]==y && surf->point[i][2]==z) break;
		}
		if(i==point_num)
		{
			surf->point[i][0] = x;
			surf->point[i][1] = y;
			surf->point[i][2] = z;
		}
	}

	return i;
}

int* find_surface_containing_point(struct surface* surf, int i0)
{
	int i;

	if(surf->dimension == 2)
	{
		int *r = new int[2];
		for(i=0; i<surf->edge_num; i++)
		{
			if(surf->edge[i][0] == i0) r[1] = i;
			else if(surf->edge[i][1] == i0) r[0] = i;
		}
		return r;
	}
	if(surf->dimension == 3)
	{

		int *r = new int[surf->surf_num + 1];
		for(i=0; i<surf->surf_num; i++) r[i] = -1;

		int num = 0;
		for(i=0; i<surf->surf_num; i++)
		{
			if(surf->surf[i][0] == i0 || surf->surf[i][1] == i0 || surf->surf[i][2] == i0) r[num++] = i;
		}

		return r;
	}

	return NULL;
}

int* find_surface_containing_point(int dimension, int face_num, int **faces, int i0, int type)
{
	int i;

	if(dimension == 2)
	{
		int *r = new int[2];
		for(i=0; i<face_num; i++)
		{
			if(faces[0][i] == i0) r[1] = i;
			else if(faces[1][i] == i0) r[0] = i;
		}
		return r;
	}
	if(dimension == 3)
	{
		int *r = new int[face_num + 1];
		for(i=0; i<face_num; i++) r[i] = -1;

		int num = 0;
		for(i=0; i<face_num; i++)
		{
			if(faces[0][i] == i0 || faces[1][i] == i0 || faces[2][i] == i0) r[num++] = i;
			if(type==1 && faces[3][i] == i0) r[num++] = i;
		}

		return r;
	}

	return NULL;
}


/*double signed_distance_face_point(struct surface* surf, int i0, double *v)
{
	if (surf->dimension==2)
	{
		int i;
		int imin = 0;
		double tempdist;
		double pointdist;

		for(i=0; i<2; i++)
		{
			tempdist = sqrt(pow((v[0]-surf->point[surf->edge[i0][i]][0]),2) + pow((v[1]-surf->point[surf->edge[i0][i]][1]),2));
			if(i==0) pointdist = tempdist;
			else if(tempdist < pointdist)
			{
				pointdist = tempdist; 
				imin = i;
			}
		}

		double *a1;
		double *a2;

		if(imin == 0)
		{
			a1 = surf->point[surf->surf[i0][0]];
			a2 = surf->point[surf->surf[i0][1]];
		}
		else
		{
			a1 = surf->point[surf->surf[i0][1]];
			a2 = surf->point[surf->surf[i0][0]];
		}

		if(vec_inner_product(2,v,a1,a1,a2) > 0) tempdist = pointdist;
		else tempdist = distance_face_point(2,v,a1,a2,NULL);

		double r = vec_rot(v,a1,a1,a2);

		if((imin==0 && r < 0) || (imin==1 && r > 0)) tempdist = -tempdist;

		return tempdist;
	}

	if (surf->dimension==3)
	{
		int i;
		int imin = 0;
		double tempdist;
		double pointdist;

		for(i=0; i<3; i++)
		{
			tempdist = sqrt(pow((v[0]-surf->point[surf->surf[i0][i]][0]),2) + pow((v[1]-surf->point[surf->surf[i0][i]][1]),2) + pow((v[2]-surf->point[surf->surf[i0][i]][2]),2));
			if(i==0) pointdist = tempdist;
			else if(tempdist < pointdist)
			{
				pointdist = tempdist; 
				imin = i;
			}
		}

		double *a1;
		double *a2;
		double *a3;

		if(imin == 0)
		{
			a1 = surf->point[surf->surf[i0][0]];
			a2 = surf->point[surf->surf[i0][1]];
			a3 = surf->point[surf->surf[i0][2]];
		}
		else if(imin == 1)
		{
			a1 = surf->point[surf->surf[i0][1]];
			a2 = surf->point[surf->surf[i0][2]];
			a3 = surf->point[surf->surf[i0][0]];
		}
		else
		{
			a1 = surf->point[surf->surf[i0][2]];
			a2 = surf->point[surf->surf[i0][0]];
			a3 = surf->point[surf->surf[i0][1]];
		}
			
		double *n = vec_curl(a1, a2, a1, a3);
		double zerovec[3] = {0,};

		double *c1 = vec_curl(v,a1,a1,a2);
		double *c2 = vec_curl(v,a1,a1,a3);

		double rot1 = vec_inner_product(3,zerovec,c1,zerovec,n);
		double rot2 = vec_inner_product(3,zerovec,c2,zerovec,n);

		double inner1 = vec_inner_product(3,v,a1,a1,a2);
		double inner2 = vec_inner_product(3,v,a1,a1,a3);

		if(inner1>0 && inner2>0) tempdist = pointdist;
		else if(rot1>0 && rot2<0) tempdist = distance_face_point(3,v,a1,a2,a3);
		else if(rot1<0) tempdist = distance_face_point(3,v,a1,a2,NULL);
		else tempdist = distance_face_point(3,v,a1,a3,NULL);
			
		if(vec_inner_product(3,v,a1,zerovec,n)<0) tempdist = -tempdist;
		delete [] n;

		return tempdist;
	}
	return 0;
}

int min_distance_surf_point(struct surface* surf, double x, double y, double z)    // test function
{
	double dist = MAXSIZE * 10;
	double tempdist;
	int imin;

	if (surf->dimension==2)
	{
		for(int i=0; i<surf->point_num; i++)
		{
			tempdist = sqrt(pow((x-surf->point[i][0]),2) + pow((y-surf->point[i][1]),2));
			if(tempdist < dist)
			{
				dist = tempdist;
				imin = i;
			}
		}
		return imin;
	}
	if (surf->dimension==3)
	{
		for(int i=0; i<surf->point_num; i++)
		{
			tempdist = sqrt(pow((x-surf->point[i][0]),2) + pow((y-surf->point[i][1]),2) + pow((z-surf->point[i][2]),2));
			if(tempdist < dist)
			{
				dist = tempdist;
				imin = i;
			}
		}
		return imin;
	}

	return 0;
}

double distance_surf_point_old(struct surface* surf, double x, double y, double z)
{
	int i;
	double dist;
	double tempdist;

	if (surf->dimension==2) // counter-clockwise
	{
		dist = MAXSIZE * 10;
		int imin;
		for(i=0; i<surf->point_num; i++)			// minimal distance point method?
		{
			tempdist = sqrt(pow((x-surf->point[i][0]),2) + pow((y-surf->point[i][1]),2));
			if(tempdist < dist)
			{
				dist = tempdist;
				imin = i;
			}
		}

		double vec[2] = {x,y};

		int *r = find_surface_containing_point(surf, imin);
		if(r==NULL)
		{
			r = NULL; // we should consider the boundary later
		}
		double *a = surf->point[surf->edge[r[0]][0]];
		double *b = surf->point[imin];
		double *c = surf->point[surf->edge[r[1]][1]];
		delete [] r;

		double dist1; double dist2;		// in this case, suff. to consider edge[2][(imin+surf->point_num-1)%surf->point_num] and edge[2][imin]
		if(vec_inner_product(2,vec,b,b,c) > 0) dist1 = dist;
		else dist1 = distance_face_point(2,vec,b,c,NULL);
		if(vec_inner_product(2,vec,b,b,a) > 0) dist2 = dist;
		else dist2 = distance_face_point(2,vec,a,b,NULL);

		dist = min(dist1,dist2);

		double r1 = vec_rot(vec,b,b,c);
		double r2 = vec_rot(vec,b,b,a);

		if(dist==dist1 && r1 > 0) return dist;
		if(dist==dist2 && r2 < 0) return dist;
		return -dist;
	}
	if (surf->dimension==3) // outer-pointing orientation
	{
		dist = MAXSIZE * 10;
		int imin;
		double pointdist = MAXSIZE * 10;			// minimal distance point method?

		for(i=0; i<surf->point_num; i++)
		{
			tempdist = sqrt(pow((x-surf->point[i][0]),2) + pow((y-surf->point[i][1]),2) + pow((z-surf->point[i][2]),2));
			if(tempdist < pointdist)
			{
				pointdist = tempdist; 
				imin = i;
			}
		}

		double vec[3] = {x,y,z};

		int *r = find_surface_containing_point(surf, imin);

		if(r==NULL)
		{
			r = NULL; // we should consider the boundary later
		}

		dist = pointdist;

		i = 0;
		while(r[i]!=-1)
		{
			double *a1;
			double *a2;
			double *a3;

			if(imin == surf->surf[r[i]][0])
			{
				a1 = surf->point[surf->surf[r[i]][0]];
				a2 = surf->point[surf->surf[r[i]][1]];
				a3 = surf->point[surf->surf[r[i]][2]];
			}
			else if(imin == surf->surf[r[i]][1])
			{
				a1 = surf->point[surf->surf[r[i]][1]];
				a2 = surf->point[surf->surf[r[i]][2]];
				a3 = surf->point[surf->surf[r[i]][0]];
			}
			else
			{
				a1 = surf->point[surf->surf[r[i]][2]];
				a2 = surf->point[surf->surf[r[i]][0]];
				a3 = surf->point[surf->surf[r[i]][1]];
			}
			
			double *n = vec_curl(a1, a2, a1, a3);
			double zerovec[3] = {0,};

			double *c1 = vec_curl(vec,a1,a1,a2);
			double *c2 = vec_curl(vec,a1,a1,a3);

			double rot1 = vec_inner_product(3,zerovec,c1,zerovec,n);
			double rot2 = vec_inner_product(3,zerovec,c2,zerovec,n);

			double inner1 = vec_inner_product(3,vec,a1,a1,a2);
			double inner2 = vec_inner_product(3,vec,a1,a1,a3);

			if(inner1>0 && inner2>0) tempdist = pointdist;
			else if(rot1>0 && rot2<0) tempdist = distance_face_point(3,vec,a1,a2,a3);
			else if(rot1<0) tempdist = distance_face_point(3,vec,a1,a2,NULL);
			else tempdist = distance_face_point(3,vec,a1,a3,NULL);
			
			if(vec_inner_product(3,vec,a1,zerovec,n)<0) tempdist = -tempdist;
			delete [] n;

			if(abs(tempdist) <= abs(dist))
			{
				dist = tempdist;
			}

			i++;
		}

		delete [] r;

		return dist;
	}

	return 0;
}

double burning_surf_point(struct surface* surf, double x,double y, double z)
{
	if(surf->burning_rate==NULL) return 0;
	
	double dist = MAXSIZE * 10;
	double tempdist;
	int imin;

	if (surf->dimension==2)
	{
		for(int i=0; i<surf->point_num; i++)
		{
			tempdist = sqrt(pow((x-surf->point[i][0]),2) + pow((y-surf->point[i][1]),2));
			if(tempdist < dist)
			{
				dist = tempdist;
				imin = i;
			}
		}

		double vec[2] = {x,y};

		int *r = find_surface_containing_point(surf, imin);
		double *a = surf->point[surf->edge[r[0]][0]];
		double *b = surf->point[imin];
		double *c = surf->point[surf->edge[r[1]][1]];
		delete [] r;

		double dist1; double dist2;		// in this case, suff. to consider edge[2][(imin+surf->point_num-1)%surf->point_num] and edge[2][imin]
		double inner1 = vec_inner_product(2,b,c,vec,b);
		double inner2 = vec_inner_product(2,a,b,vec,b);
		if(inner1 > 0) dist1 = dist;
		else dist1 = distance_face_point(2,vec,b,c,NULL);
		if(inner2 < 0) dist2 = dist;
		else dist2 = distance_face_point(2,vec,a,b,NULL);

		tempdist = min(dist1,dist2);
		if(tempdist==dist) return surf->burning_rate[imin];
		else if(tempdist==dist1)
		{
			double r = abs(inner1)/vec_inner_product(2,b,c,b,c);
			return (1-r)*surf->burning_rate[imin] + r*surf->burning_rate[(imin+1)%surf->point_num];
		}
		else
		{
			double r = abs(inner2)/vec_inner_product(2,a,b,a,b);
			return (1-r)*surf->burning_rate[imin] + r*surf->burning_rate[(imin+surf->point_num-1)%surf->point_num];
		}
	}
	if (surf->dimension==3)
	{
	}

	return 0;
}*/