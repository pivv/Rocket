#include "levelsets_2extra.h"

double level_set::level_computephi(double x, double y)
{
	struct quadtree *current1 = find_left_top_face(tree_grid_2d,x,y);
	struct quadtree *current2 = find_right_bottom_face(tree_grid_2d,x,y);
	if(current1->length < current2->length)	return computephi(current1,x,y);
	else return computephi(current2,x,y);
}

double level_set::level_computephi(double x, double y, double z)
{
	struct octree *current1 = find_left_bottom_back_face(tree_grid_3d,x,y,z);
	struct octree *current2 = find_right_top_back_face(tree_grid_3d,x,y,z);
	struct octree *current3 = find_left_bottom_front_face(tree_grid_3d,x,y,z);
	struct octree *current4 = find_right_top_front_face(tree_grid_3d,x,y,z);

	int mini = 1;
	double minleng = current1->length;
	if(current2->length < minleng)
	{
		mini = 2;
		minleng = current2->length;
	}
	if(current3->length < minleng)
	{
		mini = 3;
		minleng = current3->length;
	}
	if(current4->length < minleng)
	{
		mini = 4;
		minleng = current4->length;
	}

	if(mini==1) return computephi(current1,x,y,z);
	else if(mini==2) return computephi(current2,x,y,z);
	else if(mini==3) return computephi(current3,x,y,z);
	else return computephi(current4,x,y,z);
}


//here



double external_velocity_x_2d(double x, double y, int num)
{
	return -(y-(double)MAXSIZE/2.);
	if(num < 256) return -pow(sin(PI*x/(double)MAXSIZE),2) * sin(2*PI*y/(double)MAXSIZE);
	else return pow(sin(PI*x/(double)MAXSIZE),2) * sin(2*PI*y/(double)MAXSIZE);
}

double external_velocity_y_2d(double x, double y, int num)
{
	return (x-(double)MAXSIZE/2.);
	if(num < 256) return pow(sin(PI*y/(double)MAXSIZE),2) * sin(2*PI*x/(double)MAXSIZE);
	else return -pow(sin(PI*y/(double)MAXSIZE),2) * sin(2*PI*x/(double)MAXSIZE);
}


double external_velocity_x_3d(double x, double y, double z, int num)
{
	return 0.1;
}

double external_velocity_y_3d(double x, double y, double z, int num)
{
	return 0.1;
}

double external_velocity_z_3d(double x, double y, double z, int num)
{
	return 0.1;
}


void level_set::tree_propagate_surface_velocity(int num)
{
	if(dimension==2)
	{
		pointphi_2d *temp = tree_point_phi_2d_start;
		while(temp!=NULL)
		{
			double xhat = temp->x - time_step/2. * external_velocity_x_2d(temp->x, temp->y, num);
			double yhat = temp->y - time_step/2. * external_velocity_y_2d(temp->x, temp->y, num);
			double xd = temp->x - time_step * external_velocity_x_2d(xhat, yhat, num);
			double yd = temp->y - time_step * external_velocity_y_2d(xhat, yhat, num);
			temp->tempphi = level_computephi(xd,yd);
			temp = temp->next;
		}

		temp = tree_point_phi_2d_start;
		while(temp!=NULL)
		{
			temp->phi = temp->tempphi;
			temp = temp->next;
		}

		temp = tree_point_phi_2d_start;
		while(temp!=NULL)
		{
			quadtree_coeffs_saving(tree_grid_2d, temp);
			temp = temp->next;
		}

		iter_reinitial_scheme();
		generate_tree(tree_grid_2d, false);
	}

	if(dimension==3)
	{
		pointphi_3d *temp = tree_point_phi_3d_start;
		while(temp!=NULL)
		{
			double xhat = temp->x - time_step/2. * external_velocity_x_3d(temp->x, temp->y, temp->z, num);
			double yhat = temp->y - time_step/2. * external_velocity_y_3d(temp->x, temp->y, temp->z, num);
			double zhat = temp->z - time_step/2. * external_velocity_z_3d(temp->x, temp->y, temp->z, num);
			double xd = temp->x - time_step * external_velocity_x_3d(xhat, yhat, zhat, num);
			double yd = temp->y - time_step * external_velocity_y_3d(xhat, yhat, zhat, num);
			double zd = temp->z - time_step * external_velocity_z_3d(xhat, yhat, zhat, num);
			temp->tempphi = level_computephi(xd,yd,zd);
			temp = temp->next;
		}

		temp = tree_point_phi_3d_start;
		while(temp!=NULL)
		{
			temp->phi = temp->tempphi;
			temp = temp->next;
		}

		iter_reinitial_scheme();
		generate_tree(tree_grid_3d, false);
	}

}


void level_set::tree_propagate_surface()
{
	pointphi_2d *temp = tree_point_phi_2d_start;
	while(temp!=NULL)
	{
		//temp->tempphi = temp->phi + time_step * tree_hg(temp, false,0,0,0,0);
		temp->tempphi = temp->phi - time_step * tree_hg(temp, true,0,0,0,0);
		temp = temp->next;
	}

	temp = tree_point_phi_2d_start;
	while(temp!=NULL)
	{
		temp->tempphi2 = temp->phi;
		temp->phi = temp->tempphi;
		temp = temp->next;
	}

	temp = tree_point_phi_2d_start;
	while(temp!=NULL)
	{
		quadtree_coeffs_saving(tree_grid_2d, temp);
		temp = temp->next;
	}

	temp = tree_point_phi_2d_start;
	while(temp!=NULL)
	{
		//temp->tempphi = temp->phi + time_step * tree_hg(temp, false,0,0,0,0);
		temp->tempphi = temp->phi - time_step * tree_hg(temp, true,0,0,0,0);
		temp = temp->next;
	}

	temp = tree_point_phi_2d_start;
	while(temp!=NULL)
	{
		temp->phi = (temp->tempphi + temp->tempphi2)/2.;
		temp = temp->next;
	}

	temp = tree_point_phi_2d_start;
	while(temp!=NULL)
	{
		quadtree_coeffs_saving(tree_grid_2d, temp);
		temp = temp->next;
	}
	
	//iter_reinitial_scheme_2d();
	generate_tree(tree_grid_2d, false);

}
