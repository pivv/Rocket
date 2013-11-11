#include "levelsets_5reconst.h"

void level_set::savingtree()
{
	if(dimension==2)
	{
		int i;
		if(mat_savingtree == NULL)
		{
			mat_savingtree = new double *[MAXSAVING];
			for(i=0; i<MAXSAVING; i++) mat_savingtree[i] = new double [3];
		}
		numnode = 0;
		savingtree_recursive(tree_grid_2d);
	}
	if(dimension==3)
	{
		int i; int j; int k;
		int num = (int)pow(2,max_tree_depth);
		if(array_savingtree == NULL)
		{
			array_savingtree = new double [(num+1)*(num+1)*(num+1)];
		}

		int n = 0;
		for(i=0; i<=num; i++)
		{
			for(j=0; j<=num; j++)
			{
				for(k=0; k<=num; k++)
				{
					array_savingtree[n++] = level_computephi((double)MAXSIZE * (double)i/(double)num, (double)MAXSIZE * (double)j/(double)num, (double)MAXSIZE * (double)k/(double)num);
				}
			}
		}
		
		/*if(mat_savingtree == NULL)
		{
			mat_savingtree = new double *[MAXSAVING];
			for(i=0; i<MAXSAVING; i++) mat_savingtree[i] = new double [4];
		}

		numnode = 0;
		savingtree_recursive(tree_grid_3d);*/
	}
}

void level_set::savingtree_recursive(struct quadtree *current)
{
	if(current->tree_lefttop!=NULL)
	{
		savingtree_recursive(current->tree_lefttop);
		savingtree_recursive(current->tree_leftbottom);
		savingtree_recursive(current->tree_righttop);
		savingtree_recursive(current->tree_rightbottom);
	}
	else
	{
		mat_savingtree[numnode][0] = current->phi_lefttop->x;
		mat_savingtree[numnode][1] = current->phi_lefttop->y;
		mat_savingtree[numnode++][2] = current->phi_lefttop->phi;
		mat_savingtree[numnode][0] = current->phi_leftbottom->x;
		mat_savingtree[numnode][1] = current->phi_leftbottom->y;
		mat_savingtree[numnode++][2] = current->phi_leftbottom->phi;
		mat_savingtree[numnode][0] = current->phi_rightbottom->x;
		mat_savingtree[numnode][1] = current->phi_rightbottom->y;
		mat_savingtree[numnode++][2] = current->phi_rightbottom->phi;
		mat_savingtree[numnode][0] = current->phi_righttop->x;
		mat_savingtree[numnode][1] = current->phi_righttop->y;
		mat_savingtree[numnode++][2] = current->phi_righttop->phi;
		mat_savingtree[numnode][0] = current->phi_lefttop->x;
		mat_savingtree[numnode][1] = current->phi_lefttop->y;
		mat_savingtree[numnode++][2] = current->phi_lefttop->phi;
	}
}

void level_set::savingtree_recursive(struct octree *current)
{
	if(current->tree_leftbottomback!=NULL)
	{
		savingtree_recursive(current->tree_leftbottomback);
		savingtree_recursive(current->tree_lefttopback);
		savingtree_recursive(current->tree_rightbottomback);
		savingtree_recursive(current->tree_righttopback);
		savingtree_recursive(current->tree_leftbottomfront);
		savingtree_recursive(current->tree_lefttopfront);
		savingtree_recursive(current->tree_rightbottomfront);
		savingtree_recursive(current->tree_righttopfront);
	}
	else
	{
		mat_savingtree[numnode][0] = current->phi_leftbottomback->x;
		mat_savingtree[numnode][1] = current->phi_leftbottomback->y;
		mat_savingtree[numnode][2] = current->phi_leftbottomback->z;
		mat_savingtree[numnode++][2] = current->phi_leftbottomback->phi;
		mat_savingtree[numnode][0] = current->phi_lefttopback->x;
		mat_savingtree[numnode][1] = current->phi_lefttopback->y;
		mat_savingtree[numnode][2] = current->phi_lefttopback->z;
		mat_savingtree[numnode++][2] = current->phi_lefttopback->phi;
		mat_savingtree[numnode][0] = current->phi_rightbottomback->x;
		mat_savingtree[numnode][1] = current->phi_rightbottomback->y;
		mat_savingtree[numnode][2] = current->phi_rightbottomback->z;
		mat_savingtree[numnode++][2] = current->phi_rightbottomback->phi;
		mat_savingtree[numnode][0] = current->phi_righttopback->x;
		mat_savingtree[numnode][1] = current->phi_righttopback->y;
		mat_savingtree[numnode][2] = current->phi_righttopback->z;
		mat_savingtree[numnode++][2] = current->phi_righttopback->phi;
		mat_savingtree[numnode][0] = current->phi_leftbottomfront->x;
		mat_savingtree[numnode][1] = current->phi_leftbottomfront->y;
		mat_savingtree[numnode][2] = current->phi_leftbottomfront->z;
		mat_savingtree[numnode++][2] = current->phi_leftbottomfront->phi;
		mat_savingtree[numnode][0] = current->phi_lefttopfront->x;
		mat_savingtree[numnode][1] = current->phi_lefttopfront->y;
		mat_savingtree[numnode][2] = current->phi_lefttopfront->z;
		mat_savingtree[numnode++][2] = current->phi_lefttopfront->phi;
		mat_savingtree[numnode][0] = current->phi_rightbottomfront->x;
		mat_savingtree[numnode][1] = current->phi_rightbottomfront->y;
		mat_savingtree[numnode][2] = current->phi_rightbottomfront->z;
		mat_savingtree[numnode++][2] = current->phi_rightbottomfront->phi;
		mat_savingtree[numnode][0] = current->phi_righttopfront->x;
		mat_savingtree[numnode][1] = current->phi_righttopfront->y;
		mat_savingtree[numnode][2] = current->phi_righttopfront->z;
		mat_savingtree[numnode++][2] = current->phi_righttopfront->phi;
	}
}

