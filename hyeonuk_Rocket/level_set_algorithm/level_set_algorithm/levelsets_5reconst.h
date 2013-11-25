#include "levelsets_4reinitial.h"

void level_set::tree_propagate_velocity_fluid_points(int num)
{
	if(dimension==2)
	{
		int i;
		for(i=0; i<fluid_surface_points_num; i++)
		{
			double x = fluid_surface_points[0][i];
			double y = fluid_surface_points[1][i];

			double u = external_velocity_x_2d(x,y,num);
			double v = external_velocity_y_2d(x,y,num);

			double K1[2] = {u,v};

			K1[0] *= time_step;
			K1[1] *= time_step;

			double x_temp = x;
			double y_temp = y;

			x = x + K1[0];
			y = y + K1[1];

			u = external_velocity_x_2d(x,y,num);
			v = external_velocity_y_2d(x,y,num);

			K1[0] = u;
			K1[1] = v;

			K1[0] *= time_step;
			K1[1] *= time_step;
			
			x = 3./4.*x_temp + 1./4.*x + 1./4.*K1[0];
			y = 3./4.*y_temp + 1./4.*y + 1./4.*K1[1];

			u = external_velocity_x_2d(x,y,num);
			v = external_velocity_y_2d(x,y,num);
			
			K1[0] = u;
			K1[1] = v;

			K1[0] *= time_step;
			K1[1] *= time_step;
			
			fluid_surface_points[0][i] = 1./3.*x_temp + 2./3.*x + 2./3.*K1[0];
			fluid_surface_points[1][i] = 1./3.*y_temp + 2./3.*y + 2./3.*K1[1];
		}
	}
	if(dimension==3)
	{
		int i;
		for(i=0; i<fluid_surface_points_num; i++)
		{
			double x = fluid_surface_points[0][i];
			double y = fluid_surface_points[1][i];
			double z = fluid_surface_points[2][i];

			double u = external_velocity_x_3d(x,y,z,num);
			double v = external_velocity_y_3d(x,y,z,num);
			double w = external_velocity_z_3d(x,y,z,num);

			double K1[3] = {u,v,w};

			K1[0] *= time_step;
			K1[1] *= time_step;
			K1[2] *= time_step;

			double x_temp = x;
			double y_temp = y;
			double z_temp = z;

			x = x + K1[0];
			y = y + K1[1];
			z = z + K1[2];

			u = external_velocity_x_3d(x,y,z,num);
			v = external_velocity_y_3d(x,y,z,num);
			w = external_velocity_z_3d(x,y,z,num);

			K1[0] = u;
			K1[1] = v;
			K1[2] = w;

			K1[0] *= time_step;
			K1[1] *= time_step;
			K1[2] *= time_step;
			
			x = 3./4.*x_temp + 1./4.*x + 1./4.*K1[0];
			y = 3./4.*y_temp + 1./4.*y + 1./4.*K1[1];
			z = 3./4.*z_temp + 1./4.*z + 1./4.*K1[2];

			u = external_velocity_x_3d(x,y,z,num);
			v = external_velocity_y_3d(x,y,z,num);
			w = external_velocity_z_3d(x,y,z,num);
			
			K1[0] = u;
			K1[1] = v;
			K1[2] = w;

			K1[0] *= time_step;
			K1[1] *= time_step;
			K1[2] *= time_step;
			
			fluid_surface_points[0][i] = 1./3.*x_temp + 2./3.*x + 2./3.*K1[0];
			fluid_surface_points[1][i] = 1./3.*y_temp + 2./3.*y + 2./3.*K1[1];
			fluid_surface_points[2][i] = 1./3.*z_temp + 2./3.*z + 2./3.*K1[2];
		}
	}
}

double *convection_normal_point(struct quadtree* current, double x, double y, double r)
{
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

	double u = r*n_x;
	double v = r*n_y;

	double *velocity = new double[2];
	velocity[0] = u;
	velocity[1] = v;

	return velocity;
}

double *convection_normal_point(struct octree* current, double x, double y, double z, double r)
{
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

	double u = r*n_x;
	double v = r*n_y;
	double w = r*n_z;

	double *velocity = new double[3];
	velocity[0] = u;
	velocity[1] = v;
	velocity[2] = w;

	return velocity;
}

void level_set::tree_propagate_fluid_points(double *burning_rate)
{
	if(dimension==2)
	{
		int i;
		for(i=0; i<fluid_surface_points_num; i++)
		{
			double x = fluid_surface_points[0][i];
			double y = fluid_surface_points[1][i];
			struct quadtree *current = find_cell_containing_point(x,y);

			double r;

			if(burning_rate!=NULL) r = burning_rate[i];
			else r = 1.;

			if(rcase_surface!=NULL && level_computephi_type(x,y,0) > 0) r = 0.;

			double *K1 = convection_normal_point(current, x, y, r);
			K1[0] *= time_step;
			K1[1] *= time_step;

			double x_temp = x;
			double y_temp = y;

			x = x + K1[0];
			y = y + K1[1];

			delete K1;

			K1 = convection_normal_point(current, x, y, r);
			K1[0] *= time_step;
			K1[1] *= time_step;
			
			x = 3./4.*x_temp + 1./4.*x + 1./4.*K1[0];
			y = 3./4.*y_temp + 1./4.*y + 1./4.*K1[1];

			delete K1;

			K1 = convection_normal_point(current, x, y, r);
			K1[0] *= time_step;
			K1[1] *= time_step;
			
			fluid_surface_points[0][i] = 1./3.*x_temp + 2./3.*x + 2./3.*K1[0];
			fluid_surface_points[1][i] = 1./3.*y_temp + 2./3.*y + 2./3.*K1[1];
		}
	}
	if(dimension==3)
	{
		int i;
		for(i=0; i<fluid_surface_points_num; i++)
		{
			double x = fluid_surface_points[0][i];
			double y = fluid_surface_points[1][i];
			double z = fluid_surface_points[2][i];
			struct octree *current = find_cell_containing_point(x,y,z);

			double r;

			if(burning_rate!=NULL) r = burning_rate[i];
			else r = 1.;

			if(rcase_surface!=NULL && level_computephi_type(x,y,z,0) > 0) r = 0.;

			double *K1 = convection_normal_point(current, x, y, z, r);
			K1[0] *= time_step;
			K1[1] *= time_step;
			K1[2] *= time_step;

			double x_temp = x;
			double y_temp = y;
			double z_temp = z;

			x = x + K1[0];
			y = y + K1[1];
			z = z + K1[2];

			delete K1;

			K1 = convection_normal_point(current, x, y, z, r);
			K1[0] *= time_step;
			K1[1] *= time_step;
			K1[2] *= time_step;
			
			x = 3./4.*x_temp + 1./4.*x + 1./4.*K1[0];
			y = 3./4.*y_temp + 1./4.*y + 1./4.*K1[1];
			z = 3./4.*z_temp + 1./4.*z + 1./4.*K1[2];

			delete K1;

			K1 = convection_normal_point(current, x, y, z, r);
			K1[0] *= time_step;
			K1[1] *= time_step;
			K1[2] *= time_step;
			
			fluid_surface_points[0][i] = 1./3.*x_temp + 2./3.*x + 2./3.*K1[0];
			fluid_surface_points[1][i] = 1./3.*y_temp + 2./3.*y + 2./3.*K1[1];
			fluid_surface_points[2][i] = 1./3.*z_temp + 2./3.*z + 2./3.*K1[2];
		}
	}
}