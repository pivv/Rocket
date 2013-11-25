#include "levelsets_5reconst.h"

void level_set::tree_remesh()
{
}

double *remesh_flow_velocity(struct quadtree* current, int point_num, double **points, int face_num, int **faces, int i0, int type)
{
	int *r = find_surface_containing_point(2, face_num, faces, i0, type);

	double x = points[0][i0];
	double y = points[1][i0];
		
	double x1 = current->phi_leftbottom->x;
	double x2 = current->phi_rightbottom->x;
	double y1 = current->phi_leftbottom->y;
	double y2 = current->phi_lefttop->y;
	double l = current->length;
			
	double phi_x = quadtree_dx(current->phi_leftbottom)*(x2-x)/l*(y2-y)/l + quadtree_dx(current->phi_lefttop)*(x2-x)/l*(y-y1)/l
		+ quadtree_dx(current->phi_rightbottom)*(x-x1)/l*(y2-y)/l + quadtree_dx(current->phi_righttop)*(x-x1)/l*(y-y1)/l;
	
			
	double phi_y = quadtree_dy(current->phi_leftbottom)*(x2-x)/l*(y2-y)/l + quadtree_dy(current->phi_lefttop)*(x2-x)/l*(y-y1)/l
		+ quadtree_dy(current->phi_rightbottom)*(x-x1)/l*(y2-y)/l + quadtree_dy(current->phi_righttop)*(x-x1)/l*(y-y1)/l;

	double n = sqrt(phi_x*phi_x + phi_y*phi_y);

	double n_x = phi_x/n;
	double n_y = phi_y/n;


	double u = 0.;
	double v = 0.;

	double xx1 = points[0][faces[0][r[0]]];
	double yy1 = points[1][faces[0][r[0]]];
	double r1_2 = (x-xx1)*(x-xx1) + (y-yy1)*(y-yy1);

	u += (x-xx1)/r1_2;
	v += (y-yy1)/r1_2;

	double xx2 = points[0][faces[1][r[1]]];
	double yy2 = points[1][faces[1][r[1]]];
	double r2_2 = (x-xx2)*(x-xx2) + (y-yy2)*(y-yy2);

	u += (x-xx2)/r2_2;
	v += (y-yy2)/r2_2;

	delete [] r;

	double rr = sqrt(u*u+v*v);
	u = u/rr;
	v = v/rr;

	double *velocity = new double[2];

	velocity[0] = u - (u*n_x + v*n_y)*n_x;
	velocity[1] = v - (u*n_x + v*n_y)*n_y;

	return velocity;
}


double *remesh_flow_velocity(struct octree *current, int point_num, double **points, int face_num, int **faces, int i0, int type)
{
	int *r = find_surface_containing_point(3, face_num, faces, i0, type);

	double x = points[0][i0];
	double y = points[1][i0];
	double z = points[2][i0];
	
	double x1 = current->phi_leftbottomback->x;
	double x2 = current->phi_rightbottomback->x;
	double y1 = current->phi_leftbottomback->y;
	double y2 = current->phi_lefttopback->y;
	double z1 = current->phi_leftbottomback->z;
	double z2 = current->phi_leftbottomfront->z;
	double l = current->length;
			
	double phi_x = octree_dx(current->phi_leftbottomback)*(x2-x)/l*(y2-y)/l*(z2-z)/l + octree_dx(current->phi_lefttopback)*(x2-x)/l*(y-y1)/l*(z2-z)/l
		+ octree_dx(current->phi_rightbottomback)*(x-x1)/l*(y2-y)/l*(z2-z)/l + octree_dx(current->phi_righttopback)*(x-x1)/l*(y-y1)/l*(z2-z)/l
		+ octree_dx(current->phi_leftbottomfront)*(x2-x)/l*(y2-y)/l*(z-z1)/l + octree_dx(current->phi_lefttopfront)*(x2-x)/l*(y-y1)/l*(z-z1)/l
		+ octree_dx(current->phi_rightbottomfront)*(x-x1)/l*(y2-y)/l*(z-z1)/l + octree_dx(current->phi_righttopfront)*(x-x1)/l*(y-y1)/l*(z-z1)/l;

	double phi_y = octree_dy(current->phi_leftbottomback)*(x2-x)/l*(y2-y)/l*(z2-z)/l + octree_dy(current->phi_lefttopback)*(x2-x)/l*(y-y1)/l*(z2-z)/l
		+ octree_dy(current->phi_rightbottomback)*(x-x1)/l*(y2-y)/l*(z2-z)/l + octree_dy(current->phi_righttopback)*(x-x1)/l*(y-y1)/l*(z2-z)/l
		+ octree_dy(current->phi_leftbottomfront)*(x2-x)/l*(y2-y)/l*(z-z1)/l + octree_dy(current->phi_lefttopfront)*(x2-x)/l*(y-y1)/l*(z-z1)/l
		+ octree_dy(current->phi_rightbottomfront)*(x-x1)/l*(y2-y)/l*(z-z1)/l + octree_dy(current->phi_righttopfront)*(x-x1)/l*(y-y1)/l*(z-z1)/l;

	double phi_z = octree_dz(current->phi_leftbottomback)*(x2-x)/l*(y2-y)/l*(z2-z)/l + octree_dz(current->phi_lefttopback)*(x2-x)/l*(y-y1)/l*(z2-z)/l
		+ octree_dz(current->phi_rightbottomback)*(x-x1)/l*(y2-y)/l*(z2-z)/l + octree_dz(current->phi_righttopback)*(x-x1)/l*(y-y1)/l*(z2-z)/l
		+ octree_dz(current->phi_leftbottomfront)*(x2-x)/l*(y2-y)/l*(z-z1)/l + octree_dz(current->phi_lefttopfront)*(x2-x)/l*(y-y1)/l*(z-z1)/l
		+ octree_dz(current->phi_rightbottomfront)*(x-x1)/l*(y2-y)/l*(z-z1)/l + octree_dz(current->phi_righttopfront)*(x-x1)/l*(y-y1)/l*(z-z1)/l;

	double n = sqrt(phi_x*phi_x + phi_y*phi_y + phi_z*phi_z);

	double n_x = phi_x/n;
	double n_y = phi_y/n;
	double n_z = phi_z/n;

	double u = 0.;
	double v = 0.;
	double w = 0.;

	int i=0;
	while(r[i]!=-1)
	{
		double xx1 = points[0][faces[0][r[i]]];
		double yy1 = points[1][faces[0][r[i]]];
		double zz1 = points[2][faces[0][r[i]]];
		double r1_2 = (x-xx1)*(x-xx1) + (y-yy1)*(y-yy1) + (z-zz1)*(z-zz1);

		if(r1_2 > MINERROR)
		{
			u += (x-xx1)/r1_2;
			v += (y-yy1)/r1_2;
			w += (z-zz1)/r1_2;
		}

		double xx2 = points[0][faces[1][r[i]]];
		double yy2 = points[1][faces[1][r[i]]];
		double zz2 = points[2][faces[1][r[i]]];
		double r2_2 = (x-xx2)*(x-xx2) + (y-yy2)*(y-yy2) + (z-zz2)*(z-zz2);

		if(r2_2 > MINERROR)
		{
			u += (x-xx2)/r2_2;
			v += (y-yy2)/r2_2;
			w += (z-zz2)/r2_2;
		}

		double xx3 = points[0][faces[2][r[i]]];
		double yy3 = points[1][faces[2][r[i]]];
		double zz3 = points[2][faces[2][r[i]]];
		double r3_2 = (x-xx3)*(x-xx3) + (y-yy3)*(y-yy3) + (z-zz3)*(z-zz3);

		if(r3_2 > MINERROR)
		{
			u += (x-xx3)/r3_2;
			v += (y-yy3)/r3_2;
			w += (z-zz3)/r3_2;
		}

		if(type==1)
		{
			double xx4 = points[0][faces[3][r[i]]];
			double yy4 = points[1][faces[3][r[i]]];
			double zz4 = points[2][faces[3][r[i]]];
			double r4_2 = (x-xx4)*(x-xx4) + (y-yy4)*(y-yy4) + (z-zz4)*(z-zz4);

			if(r4_2 > MINERROR)
			{
				u += (x-xx4)/r4_2;
				v += (y-yy4)/r4_2;
				w += (z-zz4)/r4_2;
			}
		}
	}

	delete [] r;
	
	double rr = sqrt(u*u+v*v);

	u = u/rr;
	v = v/rr;
	w = w/rr;

	double *velocity = new double[3];

	velocity[0] = u - (u*n_x + v*n_y + w*n_z)*n_x;
	velocity[1] = v - (u*n_x + v*n_y + w*n_z)*n_y;
	velocity[2] = w - (u*n_x + v*n_y + w*n_z)*n_z;

	return velocity;
}

void level_set::tree_remesh_flow()
{
	if(dimension==2)
	{
		double delta = MAXSIZE * pow(2,-max_tree_depth);

		int i;
		for(i=0; i<fluid_surface_points_num; i++)
		{
			double x = fluid_surface_points[0][i];
			double y = fluid_surface_points[1][i];
			struct quadtree *current = find_cell_containing_point(x,y);

			double *velocity = remesh_flow_velocity(current, fluid_surface_points_num, fluid_surface_points, fluid_surface_edges_num, fluid_surface_edges, i, 0);
			
			fluid_surface_points[0][i] += delta * velocity[0];
			fluid_surface_points[1][i] += delta * velocity[1];

			delete [] velocity;
		}
	}
	if(dimension==3)
	{
		double delta = MAXSIZE * pow(2,-max_tree_depth);

		int i;
		for(i=0; i<fluid_surface_points_num; i++)
		{
			double x = fluid_surface_points[0][i];
			double y = fluid_surface_points[1][i];
			double z = fluid_surface_points[2][i];
			struct octree *current = find_cell_containing_point(x,y,z);

			double *velocity = remesh_flow_velocity(current, fluid_surface_points_num, fluid_surface_points, fluid_surface_faces_num, fluid_surface_faces, i, 0);
			
			fluid_surface_points[0][i] += delta * velocity[0];
			fluid_surface_points[1][i] += delta * velocity[1];

			delete [] velocity;
		}
	}
}

int add_point_2d(double **point, int *point_num, int type, double x, double y)
{
	// type == 0 : new point, type == 1 : old point

	if(type==0)
	{
		point[*point_num][0] = x;
		point[*point_num][1] = y;

		*point_num++;
	}
	if(type==1)
	{
		int i;
		for(i=*point_num-1; i>=0; i--)
		{
			if(abs(point[i][0]-x)<MINERROR && abs(point[i][1]-y)<MINERROR) break;
		}

		return i;
	}

	return 0;
}

int add_point_3d(double **point, int *point_num, int type, double x, double y, double z)
{
	// type == 0 : new point, type == 1 : old point

	if(type==0)
	{
		point[*point_num][0] = x;
		point[*point_num][1] = y;
		point[*point_num][2] = z;

		*point_num++;
	}
	if(type==1)
	{
		int i;
		for(i=*point_num-1; i>=0; i--)
		{
			if(abs(point[i][0]-x)<MINERROR && abs(point[i][1]-y)<MINERROR &&  abs(point[i][2]-z)<MINERROR) break;
		}

		return i;
	}

	return 0;
}

int *find_pattern_2d(int *zeros)
{
	// zeros : left right bottom top
	// false : - , true : +
	int type1;
	int type2;

	int fixedpt;

	int *r = new int[3];
	r[0] = type1;
	r[1] = type2;
	r[2] = fixedpt;

	return r; // type1, type2, fixedpt
}

int *find_pattern_3d(int *zeros)
{
	// zeros : left right bottom top
	// false : - , true : +
	int type1;
	int type2;

	int fixedpt;
	int fixedface;

	int *r = new int[4];
	r[0] = type1;
	r[1] = type2;
	r[2] = fixedpt;
	r[3] = fixedface;

	return r; // type1, type2, fixedpt
}


void rotate_pattern_2d(int fixedpt, int *zeros)
{
	// zeros : left right bottom top
	// false : - , true : +
}

void rotate_pattern_3d(int fixedpt, int fixedface, int *zeros)
{
	// zeros : left right bottom top
	// false : - , true : +
}


void marching_cube_pattern_2d(int type1, int type2, int **edges, int *num, int *zeros)
{
	// type1 determines pattern type
	// - + : 0 1 type2 determines which sign is at point
	// zeros : left right bottom top
	if(type1==-1);
	if(type1==0)
	{
		if(type2==0)
		{
			edges[0][*num] = zeros[0];
			edges[1][*num] = zeros[2];
			*num ++;
		}
		if(type2==1)
		{
			edges[0][*num] = zeros[2];
			edges[1][*num] = zeros[0];
			*num ++;
		}
	}
	if(type1==1)
	{
		edges[0][*num] = zeros[0];
		edges[1][*num] = zeros[1];
		*num ++;
	}
	if(type1==2)
	{
		edges[0][*num] = zeros[0];
		edges[1][*num] = zeros[2];
		*num ++;
		edges[0][*num] = zeros[1];
		edges[1][*num] = zeros[3];
		*num ++;
	}
}

void level_set::tree_reconstruct_coarse_surface(double csize)
{
	int i; int j;
	int k1; int k2;

	int max_point_num = (int)(5*10*pow(2,max_tree_depth));
	int max_edge_num = (int)(5*10*pow(2,max_tree_depth));
	int max_face_num = (int)(5*10*pow(2,max_tree_depth)*pow(2,max_tree_depth));

	if(dimension==2)
	{
		double **temppoint = new double *[2];
		for(i=0; i<2; i++) temppoint[i] = new double [max_point_num];

		int **tempedge = new int *[2];
		for(i=0; i<2; i++) tempedge[i] = new int [max_edge_num];

		int cnum = floor((double)MAXSIZE/csize);

		int point_num = 0;
		int edge_num = 0;

		for(i=0; i<cnum; i++)
		{
			for(j=0; j<cnum; j++)
			{
				double x1 = MAXSIZE * (double)i/(double)cnum;
				double x2 = MAXSIZE * ((double)i+1.)/(double)cnum;
				double y1 = MAXSIZE * (double)j/(double)cnum;
				double y2 = MAXSIZE * ((double)j+1.)/(double)cnum;

				double phi[4]; // leftbottom, rightbottom, lefttop, righttop

				phi[0] = level_computephi(x1,y1);
				phi[1] = level_computephi(x2,y1);
				phi[2] = level_computephi(x1,y2);
				phi[3] = level_computephi(x2,y2);

				bool pm[4];
				int minus_num = 0;
				int m;
				for (m=0; m<4; m++)
				{
					if(phi[0]<0)
					{
						pm[m] = false;
						minus_num ++;
					}
					else pm[m] = true;
				}

				int zeros[4]; // left, right, top, bottom

				double phi_leftbottom = phi[0];
				double phi_rightbottom = phi[1];
				double phi_lefttop = phi[2];
				double phi_righttop = phi[3];

				if(phi_leftbottom*phi_lefttop<0)
				{
					if(i==0) zeros[0] = add_point_2d(temppoint,&point_num,0,x,y);
					else zeros[0] = add_point_2d(temppoint,&point_num,1,x,y);
				}
				if(phi_leftbottom*phi_rightbottom<0)
				{
					if(j==0) zeros[2] = add_point_2d(temppoint,&point_num,0,x,y);
					else zeros[2] = add_point_2d(temppoint,&point_num,1,x,y);
				}
				if(phi_lefttop*phi_righttop<0)
				{
					zeros[3] = add_point_2d(temppoint,&point_num,0,x,y);
				}
				if(phi_rightbottom*phi_righttop<0)
				{
					zeros[1] = add_point_2d(temppoint,&point_num,0,x,y);
				}

				int *r = find_pattern_2d(zeros); // type1, type2, fixedpt

				rotate_pattern_2d(r[2], zeros);
				marching_cube_pattern_2d(r[0], r[1], tempedge, &edge_num, zeros);

				delete r;
			}
		}
	}
	if(dimension==3)
	{
		double **temppoint = new double *[3];
		for(i=0; i<3; i++) temppoint[i] = new double [max_point_num];

		int **tempface = new int *[3];
		for(i=0; i<3; i++) tempface[i] = new int [max_face_num];

		int cnum = floor((double)MAXSIZE/csize);

		int point_num = 0;
		int face_num = 0;

		for(i=0; i<cnum; i++)
		{
			for(j=0; j<cnum; j++)
			{
				double x1 = MAXSIZE * (double)i/(double)cnum;
				double x2 = MAXSIZE * ((double)i+1.)/(double)cnum;
				double y1 = MAXSIZE * (double)j/(double)cnum;
				double y2 = MAXSIZE * ((double)j+1.)/(double)cnum;

				double phi[4]; // leftbottom, rightbottom, lefttop, righttop

				phi[0] = level_computephi(x1,y1);
				phi[1] = level_computephi(x2,y1);
				phi[2] = level_computephi(x1,y2);
				phi[3] = level_computephi(x2,y2);

				bool pm[4];
				int minus_num = 0;
				int m;
				for (m=0; m<4; m++)
				{
					if(phi[0]<0)
					{
						pm[m] = false;
						minus_num ++;
					}
					else pm[m] = true;
				}

				int zeros[4]; // left, right, top, bottom

				double phi_leftbottom = phi[0];
				double phi_rightbottom = phi[1];
				double phi_lefttop = phi[2];
				double phi_righttop = phi[3];

				if(phi_leftbottom*phi_lefttop<0)
				{
					if(i==0) zeros[0] = add_point_2d(temppoint,&point_num,0,x,y);
					else zeros[0] = add_point_2d(temppoint,&point_num,1,x,y);
				}
				if(phi_leftbottom*phi_rightbottom<0)
				{
					if(j==0) zeros[2] = add_point_2d(temppoint,&point_num,0,x,y);
					else zeros[2] = add_point_2d(temppoint,&point_num,1,x,y);
				}
				if(phi_lefttop*phi_righttop<0)
				{
					zeros[3] = add_point_2d(temppoint,&point_num,0,x,y);
				}
				if(phi_rightbottom*phi_righttop<0)
				{
					zeros[1] = add_point_2d(temppoint,&point_num,0,x,y);
				}

				int *r = find_pattern_2d(zeros); // type1, type2, fixedpt

				rotate_pattern_2d(r[2], zeros);
				marching_cube_pattern_2d(r[0], r[1], tempface, &face_num, zeros);

				delete r;
			}
		}
	}
}


void level_set::iter_tree_reconstruct_surface()
{
	int i;
	int max_point_num = (int)(5*10*pow(2,max_tree_depth));
	int max_edge_num = (int)(5*10*pow(2,max_tree_depth));
	struct surface tempsurf(2, max_point_num, max_edge_num);
	int point_num = 0;
	int edge_num = 0;
	tree_reconstruct_surface_2d(tree_grid_2d, &tempsurf, point_num, edge_num);

	
	for(i=0; i<level_surface->point_num; i++) delete [] level_surface->point[i];
	delete [] level_surface->point;
	for(i=0; i<level_surface->edge_num; i++) delete [] level_surface->edge[i];
	delete [] level_surface->edge;

	level_surface->point_num = point_num;
	level_surface->point = new double *[point_num];
	for(i=0; i<point_num; i++) level_surface->point[i] = new double[dimension];
	level_surface->edge_num = edge_num;
	level_surface->edge = new int *[edge_num];
	for(i=0; i<edge_num; i++) level_surface->edge[i] = new int[2];

	for(i=0; i<point_num; i++)
	{
		level_surface->point[i][0] = tempsurf.point[i][0];
		level_surface->point[i][1] = tempsurf.point[i][1];
	}

	for(i=0; i<edge_num; i++)
	{
		level_surface->edge[i][0] = tempsurf.edge[i][0];
		level_surface->edge[i][1] = tempsurf.edge[i][1];
	}
}

void level_set::tree_reconstruct_surface_2d(struct quadtree *current, struct surface *tempsurf, int &point_num, int &edge_num)
{
	if(current->tree_lefttop!=NULL)
	{
		tree_reconstruct_surface_2d(current->tree_lefttop, tempsurf, point_num, edge_num);
		tree_reconstruct_surface_2d(current->tree_leftbottom, tempsurf, point_num, edge_num);
		tree_reconstruct_surface_2d(current->tree_righttop, tempsurf, point_num, edge_num);
		tree_reconstruct_surface_2d(current->tree_rightbottom, tempsurf, point_num, edge_num);
	}
	else
	{
		int k1; int k2;

		double phi_lefttop = current->phi_lefttop->phi;
		double phi_leftbottom = current->phi_leftbottom->phi;
		double phi_righttop = current->phi_righttop->phi;
		double phi_rightbottom = current->phi_rightbottom->phi;
		
		double x = current->phi_lefttop->x;
		double y = current->phi_lefttop->y;
		double l = current->length;

		double x_left = x;
		double y_left = (abs(phi_lefttop)*(y-l) + abs(phi_leftbottom)*y)/(abs(phi_lefttop) + abs(phi_leftbottom));
		double x_right = x+l;
		double y_right = (abs(phi_righttop)*(y-l) + abs(phi_rightbottom)*y)/(abs(phi_righttop) + abs(phi_rightbottom));
		double x_top = (abs(phi_lefttop)*(x+l) + abs(phi_righttop)*x)/(abs(phi_lefttop) + abs(phi_righttop));
		double y_top = y;
		double x_bottom = (abs(phi_leftbottom)*(x+l) + abs(phi_rightbottom)*x)/(abs(phi_leftbottom) + abs(phi_rightbottom));
		double y_bottom = y-l;

		//double x_left = x;
		//double y_left = y-l/2.;
		//double x_right = x+l;
		//double y_right = y-l/2.;
		//double x_top = x+l/2.;
		//double y_top = y;
		//double x_bottom = x+l/2.;
		//double y_bottom = y-l;

		if(phi_leftbottom>0 && phi_rightbottom>0 && phi_lefttop>0 && phi_righttop>0);			//1
		else if(phi_leftbottom<0 && phi_rightbottom<0 && phi_lefttop<0 && phi_righttop<0);		//2
		else if(l != MAXSIZE * pow(2,-max_tree_depth))
		{
			double d = 5;
		}
		else if(phi_leftbottom>=0 && phi_rightbottom>=0 && phi_lefttop>=0 && phi_righttop<=0)		//3
		{
			// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
			k1 = add_point(tempsurf, point_num, x_right, y_right, 0);
			if(k1==point_num) point_num++;
			k2 = add_point(tempsurf, point_num, x_top, y_top, 0);
			if(k2==point_num) point_num++;
			tempsurf->edge[edge_num][0] = k1;
			tempsurf->edge[edge_num][1] = k2;
			edge_num++;
		}
		else if(phi_leftbottom>=0 && phi_rightbottom>=0 && phi_lefttop<=0 && phi_righttop>=0)		//4
		{
			// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
			k1 = add_point(tempsurf, point_num, x_top, y_top, 0);
			if(k1==point_num) point_num++;
			k2 = add_point(tempsurf, point_num, x_left, y_left, 0);
			if(k2==point_num) point_num++;
			tempsurf->edge[edge_num][0] = k1;
			tempsurf->edge[edge_num][1] = k2;
			edge_num++;
		}
		else if(phi_leftbottom>=0 && phi_rightbottom>=0 && phi_lefttop<=0 && phi_righttop<=0)		//5
		{
			// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
			k1 = add_point(tempsurf, point_num, x_right, y_right, 0);
			if(k1==point_num) point_num++;
			k2 = add_point(tempsurf, point_num, x_left, y_left, 0);
			if(k2==point_num) point_num++;
			tempsurf->edge[edge_num][0] = k1;
			tempsurf->edge[edge_num][1] = k2;
			edge_num++;
		}
		else if(phi_leftbottom>=0 && phi_rightbottom<=0 && phi_lefttop>=0 && phi_righttop>=0)		//6
		{
			// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
			k1 = add_point(tempsurf, point_num, x_bottom, y_bottom, 0);
			if(k1==point_num) point_num++;
			k2 = add_point(tempsurf, point_num, x_right, y_right, 0);
			if(k2==point_num) point_num++;
			tempsurf->edge[edge_num][0] = k1;
			tempsurf->edge[edge_num][1] = k2;
			edge_num++;
		}
		else if(phi_leftbottom>=0 && phi_rightbottom<=0 && phi_lefttop>=0 && phi_righttop<=0)		//7
		{
			// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
			k1 = add_point(tempsurf, point_num, x_bottom, y_bottom, 0);
			if(k1==point_num) point_num++;
			k2 = add_point(tempsurf, point_num, x_top, y_top, 0);
			if(k2==point_num) point_num++;
			tempsurf->edge[edge_num][0] = k1;
			tempsurf->edge[edge_num][1] = k2;
			edge_num++;
		}
		else if(phi_leftbottom>=0 && phi_rightbottom<=0 && phi_lefttop<=0 && phi_righttop>=0)		//8
		{
			// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
			k1 = add_point(tempsurf, point_num, x_bottom, y_bottom, 0);
			if(k1==point_num) point_num++;
			k2 = add_point(tempsurf, point_num, x_left, y_left, 0);
			if(k2==point_num) point_num++;
			tempsurf->edge[edge_num][0] = k1;
			tempsurf->edge[edge_num][1] = k2;
			edge_num++;
					
			k1 = add_point(tempsurf, point_num, x_top, y_top, 0);
			if(k1==point_num) point_num++;
			k2 = add_point(tempsurf, point_num, x_right, y_right, 0);
			if(k2==point_num) point_num++;
			tempsurf->edge[edge_num][0] = k1;
			tempsurf->edge[edge_num][1] = k2;
			edge_num++;
		}
		else if(phi_leftbottom>=0 && phi_rightbottom<=0 && phi_lefttop<=0 && phi_righttop<=0)		//9
		{
			// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
			k1 = add_point(tempsurf, point_num, x_bottom, y_bottom, 0);
			if(k1==point_num) point_num++;
			k2 = add_point(tempsurf, point_num, x_left, y_left, 0);
			if(k2==point_num) point_num++;
			tempsurf->edge[edge_num][0] = k1;
			tempsurf->edge[edge_num][1] = k2;
			edge_num++;
		}
		else if(phi_leftbottom<=0 && phi_rightbottom>=0 && phi_lefttop>=0 && phi_righttop>=0)		//10
		{
			// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
			k1 = add_point(tempsurf, point_num, x_left, y_left, 0);
			if(k1==point_num) point_num++;
			k2 = add_point(tempsurf, point_num, x_bottom, y_bottom, 0);
			if(k2==point_num) point_num++;
			tempsurf->edge[edge_num][0] = k1;
			tempsurf->edge[edge_num][1] = k2;
			edge_num++;
		}
		else if(phi_leftbottom<=0 && phi_rightbottom>=0 && phi_lefttop>=0 && phi_righttop<=0)		//11
		{
			// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
			k1 = add_point(tempsurf, point_num, x_left, y_left, 0);
			if(k1==point_num) point_num++;
			k2 = add_point(tempsurf, point_num, x_bottom, y_bottom, 0);
			if(k2==point_num) point_num++;
			tempsurf->edge[edge_num][0] = k1;
			tempsurf->edge[edge_num][1] = k2;
			edge_num++;
					
			k1 = add_point(tempsurf, point_num, x_right, y_right, 0);
			if(k1==point_num) point_num++;
			k2 = add_point(tempsurf, point_num, x_top, y_top, 0);
			if(k2==point_num) point_num++;
			tempsurf->edge[edge_num][0] = k1;
			tempsurf->edge[edge_num][1] = k2;
			edge_num++;
		}
		else if(phi_leftbottom<=0 && phi_rightbottom>=0 && phi_lefttop<=0 && phi_righttop>=0)		//12
		{
			// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
			k1 = add_point(tempsurf, point_num, x_top, y_top, 0);
			if(k1==point_num) point_num++;
			k2 = add_point(tempsurf, point_num, x_bottom, y_bottom, 0);
			if(k2==point_num) point_num++;
			tempsurf->edge[edge_num][0] = k1;
			tempsurf->edge[edge_num][1] = k2;
			edge_num++;
		}
		else if(phi_leftbottom<=0 && phi_rightbottom>=0 && phi_lefttop<=0 && phi_righttop<=0)		//13
		{
			// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
			k1 = add_point(tempsurf, point_num, x_right, y_right, 0);
			if(k1==point_num) point_num++;
			k2 = add_point(tempsurf, point_num, x_bottom, y_bottom, 0);
			if(k2==point_num) point_num++;
			tempsurf->edge[edge_num][0] = k1;
			tempsurf->edge[edge_num][1] = k2;
			edge_num++;
		}
		else if(phi_leftbottom<=0 && phi_rightbottom<=0 && phi_lefttop>=0 && phi_righttop>=0)		//14
		{
			// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
			k1 = add_point(tempsurf, point_num, x_left, y_left, 0);
			if(k1==point_num) point_num++;
			k2 = add_point(tempsurf, point_num, x_right, y_right, 0);
			if(k2==point_num) point_num++;
			tempsurf->edge[edge_num][0] = k1;
			tempsurf->edge[edge_num][1] = k2;
			edge_num++;
		}
		else if(phi_leftbottom<=0 && phi_rightbottom<=0 && phi_lefttop>=0 && phi_righttop<=0)		//15
		{
			// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
			k1 = add_point(tempsurf, point_num, x_left, y_left, 0);
			if(k1==point_num) point_num++;
			k2 = add_point(tempsurf, point_num, x_top, y_top, 0);
			if(k2==point_num) point_num++;
			tempsurf->edge[edge_num][0] = k1;
			tempsurf->edge[edge_num][1] = k2;
			edge_num++;
		}
		else if(phi_leftbottom<=0 && phi_rightbottom<=0 && phi_lefttop<=0 && phi_righttop>=0)		//16
		{
			// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
			k1 = add_point(tempsurf, point_num, x_top, y_top, 0);
			if(k1==point_num) point_num++;
			k2 = add_point(tempsurf, point_num, x_right, y_right, 0);
			if(k2==point_num) point_num++;
			tempsurf->edge[edge_num][0] = k1;
			tempsurf->edge[edge_num][1] = k2;
			edge_num++;
		}
	}
}
