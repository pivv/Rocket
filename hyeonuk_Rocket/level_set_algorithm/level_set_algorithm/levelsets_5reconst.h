#include "levelsets_4reinitial.h"

void level_set::iter_tree_reconstruct_surface()
{
	int i;
	int max_point_num = (int)(5*10*pow(2,max_tree_depth));
	int max_edge_num = (int)(5*10*pow(2,max_tree_depth));
	struct surface tempsurf(2, max_point_num, max_edge_num, false);
	int point_num = 0;
	int edge_num = 0;
	tree_reconstruct_surface_2d(tree_grid_2d, &tempsurf, point_num, edge_num);

	
	for(i=0; i<level_surface->point_num; i++) delete [] level_surface->point[i];
	delete [] level_surface->point;
	for(i=0; i<level_surface->edge_num; i++) delete [] level_surface->edge[i];
	delete [] level_surface->edge;
	if(level_surface->burning_rate!=NULL) delete [] level_surface->burning_rate;

	level_surface->point_num = point_num;
	level_surface->point = new double *[point_num];
	for(i=0; i<point_num; i++) level_surface->point[i] = new double[dimension];
	level_surface->edge_num = edge_num;
	level_surface->edge = new int *[edge_num];
	for(i=0; i<edge_num; i++) level_surface->edge[i] = new int[2];
	level_surface->burning_rate = NULL;

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
