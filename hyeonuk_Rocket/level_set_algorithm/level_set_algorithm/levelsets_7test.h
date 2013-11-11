#include "levelsets_6saving.h"

double level_set::test_volume(struct surface *after_surface)
{
	double volume1 = 0.;
	double volume2 = 0.;

	/*int i;
	int j;
	for(i=0; i<pow(2,max_tree_depth); i++)
	{
		for(j=0; j<pow(2,max_tree_depth); j++)
		{
			double x = ((double)i+1./2.) * pow(2,-max_tree_depth) * MAXSIZE;
			double y = ((double)j+1./2.) * pow(2,-max_tree_depth) * MAXSIZE;
			double v[2] = {x,y};

			if(line_surface_intersecting(after_surface, v) <= 0) volume2 += pow(pow(2,-max_tree_depth) * MAXSIZE, 2);
		}
	}*/

    //volume2 = PI * (0.15 * MAXSIZE) * (0.15 * MAXSIZE);

	//test_volume_recursive(tree_grid_2d, &volume1, &volume2);

	volume2 = PI * (0.45 * MAXSIZE) * (0.45 * MAXSIZE);

	volume1 = 0.;

	double *point = level_surface->point[0];
	int num = 0;
	int num2;
	for(int i=0; i<level_surface->edge_num; i++)
	{
		int *r = find_surface_containing_point(level_surface, num);
		if(level_surface->edge[r[0]][0]==num) num2 = level_surface->edge[r[0]][1];
		else num2 = level_surface->edge[r[1]][1];
		volume1 += level_surface->point[num][0] * level_surface->point[num2][1] - level_surface->point[num][1] * level_surface->point[num2][0];
		num = num2;
		delete r;
	}

	volume1 = abs(volume1)/2.;

	return 100. * (volume2 - volume1)/volume2;
}

void level_set::test_volume_recursive(struct quadtree *current, double *d, double *e)
{
	if(current->tree_lefttop!=NULL)
	{
		test_volume_recursive(current->tree_lefttop, d, e);
		test_volume_recursive(current->tree_leftbottom, d, e);
		test_volume_recursive(current->tree_righttop, d, e);
		test_volume_recursive(current->tree_rightbottom, d, e);
	}
	else
	{
		double x = current->phi_lefttop->x;
		double y = current->phi_lefttop->y;
		double l = current->length;
		
		double phi_middle = computephi(current,x+l/2.,y-l/2.);
		if(phi_middle <= 0) *d = *d + pow(l,2);

		double r1 = distance_circle_point(x, y, 0.15, 0.5, 0.5);
		double r2 = distance_circle_point(x+l/2., y, 0.15, 0.5, 0.5);
		double r3 = distance_circle_point(x, y-l/2., 0.15, 0.5, 0.5);
		double r4 = distance_circle_point(x+l/2., y-l/2., 0.15, 0.5, 0.5);
        if(r1 <= 0 || r2 <=0 || r3 <=0 || r4 <=0) *e = *e + pow(l,2);
	}
}

double level_set::test_l1(struct surface *after_surface)
{
    double l1 = 0.;
	double v = 0.;
    test_l1_recursive(tree_grid_2d, after_surface, &l1, &v);
    return l1/v;
}

void level_set::test_l1_recursive(struct quadtree *current, struct surface *surf, double *d, double *e)
{
	if(current->tree_lefttop!=NULL)
	{
		test_l1_recursive(current->tree_lefttop, surf, d, e);
		test_l1_recursive(current->tree_leftbottom, surf, d, e);
		test_l1_recursive(current->tree_righttop, surf, d, e);
		test_l1_recursive(current->tree_rightbottom, surf, d, e);
	}
	else
	{
		double x = current->phi_lefttop->x;
		double y = current->phi_lefttop->y;
		double l = current->length;
		double phi_middle = computephi(current,x+l/2.,y-l/2.);
		//double dist = distance_circle_point(x+l/2.,y-l/2., 0.45, 0.5, 0.5);
		double dist = distance_surf_point(surf,x+l/2.,y-l/2.,0);
        if(abs(phi_middle) < sqrt(2) * l)
		{
			*d = *d + abs(phi_middle-dist) * pow(l,2);
			*e = *e + pow(l,2);
		}
	}
}

double level_set::test_linf(struct surface *after_surface)
{
    double linf = 0.;
	double v = 0.;
    test_linf_recursive(tree_grid_2d, after_surface, &linf, &v);
    return linf/v;
}

void level_set::test_linf_recursive(struct quadtree *current, struct surface *surf, double *d, double *e)
{
	if(current->tree_lefttop!=NULL)
	{
		test_linf_recursive(current->tree_lefttop, surf, d, e);
		test_linf_recursive(current->tree_leftbottom, surf, d, e);
		test_linf_recursive(current->tree_righttop, surf, d, e);
		test_linf_recursive(current->tree_rightbottom, surf, d, e);
	}
	else
	{
		double x = current->phi_lefttop->x;
		double y = current->phi_lefttop->y;
		double l = current->length;
		double phi_middle = computephi(current,x+l/2.,y-l/2.);
		//double dist = distance_circle_point(x+l/2.,y-l/2., 0.45, 0.5, 0.5);
		double dist = distance_surf_point(surf,x+l/2.,y-l/2., 0);
        if(abs(phi_middle) < sqrt(2) * l && *d < abs(phi_middle-dist))
		{
			*d = abs(phi_middle-dist);
			*e = *e + pow(l,2);
		}
	}
}