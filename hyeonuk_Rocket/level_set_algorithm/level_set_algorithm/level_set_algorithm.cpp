#include "levelsets_7test.h"

#include <cstring>
#include <fstream>
#include <sstream>

using namespace std;
using std::cout;
using std::cin;


string int2str(int n)
{
	stringstream ss;
	ss << n;
	return ss.str();
}


int main()
{
	int dimension;

	cout << "dimension is : ";
	cin >> dimension;
	cout << endl;

	int i; int j; //int k;
	int surface_point_num;

	srand((unsigned int) time(NULL));
	rand()/RAND_MAX;
	rand()/RAND_MAX;
	
	
	if(dimension==2)
	{
		surface_point_num = 500;

		struct surface surface1(dimension, surface_point_num, surface_point_num, false);

		surface1.dimension = dimension;
		surface1.point_num = surface_point_num;
		surface1.edge_num = surface_point_num;

		for(i=0; i<surface_point_num;i++)
		{
			surface1.point[i][0] = 0.5 + 0.15*cos(2*PI*i/(double)surface_point_num);
			surface1.point[i][1] = 0.75 + 0.15*sin(2*PI*i/(double)surface_point_num);
		}

		for(i=0; i<surface_point_num; i++)
		{
			surface1.edge[i][0] = i;
			surface1.edge[i][1] = (i+1)%surface_point_num;
		}

		//for(i=0; i<surface_point_num; i++) surface1.burning_rate[i] = 3 + 2*cos(2*PI*i/(double)surface_point_num);




		//for(i=0; i<surface_point_num; i++) surface1.burning_rate[i] = -1;

		//for(i=0; i<surface_point_num; i++) cout<<surface1.point[0][i]<<"   "<<surface1.point[1][i]<<endl;

		//cout<<distance_surf_point(&surface1,60,60)<<endl;



		struct surface surface2(dimension, surface_point_num, surface_point_num, false);
	

		surface2.dimension = dimension;
		surface2.point_num = surface_point_num;
		surface2.edge_num = surface_point_num;

		for(i=0; i<surface_point_num/4;i++)
		{
			surface2.point[i][0] = 0.25-MINERROR + (0.75-0.25)*(double)i*4./(double)surface_point_num;
			surface2.point[i][1] = 0.25-MINERROR;

			surface2.point[surface_point_num/4 + i][0] = 0.75-MINERROR;
			surface2.point[surface_point_num/4 + i][1] = 0.25-MINERROR + (0.75-0.25)*(double)i*4./(double)surface_point_num;

			surface2.point[surface_point_num/2 + i][0] = 0.75-MINERROR + (0.25-0.75)*(double)i*4./(double)surface_point_num;
			surface2.point[surface_point_num/2 + i][1] = 0.75-MINERROR;

			surface2.point[3*surface_point_num/4 + i][0] = 0.25-MINERROR;
			surface2.point[3*surface_point_num/4 + i][1] = 0.75-MINERROR + (0.25-0.75)*(double)i*4./(double)surface_point_num;
		}

		for(i=0; i<surface_point_num; i++)
		{
			surface2.edge[i][0] = i;
			surface2.edge[i][1] = (i+1)%surface_point_num;
		}


	
		struct surface surface3(dimension, surface_point_num, surface_point_num, false);

		surface3.dimension = dimension;
		surface3.point_num = surface_point_num;
		surface3.edge_num = surface_point_num;

		for(i=0; i<surface_point_num;i++)
		{
			if(i<surface_point_num/2)
			{
				surface3.point[i][0] = 0.4 + 0.1*cos(2*PI*i/(double)(surface_point_num/2));
				surface3.point[i][1] = 0.4 + 0.1*sin(2*PI*i/(double)(surface_point_num/2));
			}
			else
			{
				surface3.point[i][0] = 0.6 + 0.1*cos(2*PI*(i-surface_point_num/2)/(double)(surface_point_num/2));
				surface3.point[i][1] = 0.6 + 0.1*sin(2*PI*(i-surface_point_num/2)/(double)(surface_point_num/2));
			}
		}

		for(i=0; i<surface_point_num; i++)
		{
			if(i<surface_point_num/2)
			{
				surface3.edge[i][0] = i;
				surface3.edge[i][1] = (i+1)%(surface_point_num/2);
			}
			else
			{
				surface3.edge[i][0] = i;
				surface3.edge[i][1] = surface_point_num/2 + (i+1)%(surface_point_num/2);
			}
		}

	
		//for(i=0; i<surface_point_num; i++)
		//{
		//	if(i<surface_point_num/2) surface3.burning_rate[i] = 3;
		//	else surface3.burning_rate[i] = 2;
		//}

		struct surface surface4(dimension, surface_point_num, surface_point_num, false);

		surface4.dimension = dimension;
		surface4.point_num = surface_point_num;
		surface4.edge_num = surface_point_num;

		for(i=0; i<surface_point_num;i++)
		{
			surface4.point[i][0] = 0.5 + (0.15 + 0.1 * sin(7.*2.*PI*i/(double)surface_point_num))*cos(2*PI*i/(double)surface_point_num);
			surface4.point[i][1] = 0.5 + (0.15 + 0.1 * sin(7.*2.*PI*i/(double)surface_point_num))*sin(2*PI*i/(double)surface_point_num);
		}

		for(i=0; i<surface_point_num; i++)
		{
			surface4.edge[i][0] = i;
			surface4.edge[i][1] = (i+1)%surface_point_num;
		}

		
		struct surface surface5(dimension, surface_point_num, surface_point_num, false);

		surface5.dimension = dimension;
		surface5.point_num = surface_point_num;
		surface5.edge_num = surface_point_num;

		for(i=0; i<surface_point_num;i++)
		{
			surface5.point[i][0] = 0.5 + 0.15*cos(2*PI*i/(double)surface_point_num);
			surface5.point[i][1] = 0.5 + 0.15*sin(2*PI*i/(double)surface_point_num);
		}

		for(i=0; i<surface_point_num; i++)
		{
			surface5.edge[i][0] = i;
			surface5.edge[i][1] = (i+1)%surface_point_num;
		}






		int max_dim;

		cout << "maximum dimension is : ";
		cin >> max_dim;
		cout << endl;

		double time_step = MAXSIZE * pow(2,-max_dim) / 2.;
		//double time_step = MAXSIZE * pow(2,-max_dim) * 4.;
		//double time_step = MAXSIZE * pow(2,-max_dim) * 0.2;
		//double time_step = MAXSIZE * pow(2,-max_dim) * 2.*PI;
		//double time_step = 4./500.;

		class level_set level1(dimension, max_dim, time_step, &surface2, NULL);

		cout << "level set is constructed" << endl;
	
	
		


		ofstream myfile;

		size_t currentSize = getPeakRSS();
		size_t peakSize = getCurrentRSS();

		cout << currentSize << "   " << peakSize << endl;

		string myString1 = "tree_level";
		string myString2 = ".dat";
		string myString3 = "_surface";
		string tempString;


		
		myfile.open("tree_level_initial.dat");

		level1.savingtree();

		for(i=0; i<level1.numnode; i++)
		{
			myfile << level1.mat_savingtree[i][0] << " " << level1.mat_savingtree[i][1] << " " << level1.mat_savingtree[i][2];
			myfile << endl;
		}

		myfile.close();
		
		myfile.open("tree_level_surface_initial.dat");

		for(j=0; j<level1.level_surface->edge_num; j++)
		{
			myfile << level1.level_surface->point[level1.level_surface->edge[j][0]][0] << " " << level1.level_surface->point[level1.level_surface->edge[j][0]][1];
			myfile << endl;			
			myfile << level1.level_surface->point[level1.level_surface->edge[j][1]][0] << " " << level1.level_surface->point[level1.level_surface->edge[j][1]][1];
			myfile << endl;
		}

		myfile.close();


		
		//level1.signphi();

		/*level1.iter_reinitial_scheme();
		cout << "level set is reinitialized" << endl;

		myfile.open("tree_level_initial.dat");

		level1.savingtree();

		for(i=0; i<level1.numnode; i++)
		{
			myfile << level1.mat_savingtree[i][0] << " " << level1.mat_savingtree[i][1] << " " << level1.mat_savingtree[i][2];
			myfile << endl;
		}

		myfile.close();*/


		/*level1.iter_tree_reconstruct_surface();
		cout << "surface is reconstructed" << endl;

		myfile.open("tree_level_surface_initial.dat");

		for(j=0; j<level1.level_surface->edge_num; j++)
		{
			myfile << level1.level_surface->point[level1.level_surface->edge[j][0]][0] << " " << level1.level_surface->point[level1.level_surface->edge[j][0]][1];
			myfile << endl;			
			myfile << level1.level_surface->point[level1.level_surface->edge[j][1]][0] << " " << level1.level_surface->point[level1.level_surface->edge[j][1]][1];
			myfile << endl;
		}

		myfile.close();*/


		float time;
		BOOL err;

		CHECK_TIME_START;


		myfile.open("tree_reinitial.dat");

		myfile.close();

		for(i=0; i<=100; i++) //pow(2,max_dim)
		{
			if(i%1==0)
			{
				level1.savingtree();

				myfile.open(myString1+int2str(i/1)+myString2);

				for(j=0; j<level1.numnode; j++)
				{
					myfile << level1.mat_savingtree[j][0] << " " << level1.mat_savingtree[j][1] << " " << level1.mat_savingtree[j][2];
					myfile << endl;
				}

				myfile.close();


				level1.iter_tree_reconstruct_surface();
				cout << "surface is reconstructed" << endl;

				currentSize = getPeakRSS();
				peakSize = getCurrentRSS();

				cout << currentSize << "   " << peakSize << endl;

				myfile.open(myString1+myString3+int2str(i/1)+myString2);

				for(j=0; j<level1.level_surface->edge_num; j++)
				{
					myfile << level1.level_surface->point[level1.level_surface->edge[j][0]][0] << " " << level1.level_surface->point[level1.level_surface->edge[j][0]][1];
					myfile << endl;			
					myfile << level1.level_surface->point[level1.level_surface->edge[j][1]][0] << " " << level1.level_surface->point[level1.level_surface->edge[j][1]][1];
					myfile << endl;
				}

				myfile.close();
			}

			level1.tree_propagate_surface_velocity(i);
			//level1.tree_propagate_surface();
			cout << "level set is propagated " << i+1 << " times" << endl;

			//currentSize = getPeakRSS();
			//peakSize = getCurrentRSS();

			//cout << currentSize << "   " << peakSize << endl;

		}
		

		CHECK_TIME_END(time, err);


		/*myfile.open("tree_level_result.dat");

		level1.savingtree();

		for(i=0; i<level1.numnode; i++)
		{
			myfile << level1.mat_savingtree[i][0] << " " << level1.mat_savingtree[i][1] << " " << level1.mat_savingtree[i][2];
			myfile << endl;
		}

		myfile.close();
		
		level1.iter_tree_reconstruct_surface();
		myfile.open("tree_level_surface_result.dat");

		for(j=0; j<level1.level_surface->edge_num; j++)
		{
			myfile << level1.level_surface->point[level1.level_surface->edge[j][0]][0] << " " << level1.level_surface->point[level1.level_surface->edge[j][0]][1];
			myfile << endl;			
			myfile << level1.level_surface->point[level1.level_surface->edge[j][1]][0] << " " << level1.level_surface->point[level1.level_surface->edge[j][1]][1];
			myfile << endl;
		}

		myfile.close();*/


		double volume_error = level1.test_volume(&surface4);
		double l1_error = level1.test_l1(&surface4);
		double linf_error = level1.test_linf(&surface4);

		cout << endl;
		cout << "resolution is " << max_dim << endl;
		cout << "volume error is " << volume_error << endl;
		cout << "l1 error is " << l1_error << endl;
		cout << "linf error is " << linf_error << endl;
		cout << "time is " << time << endl;
		

		double end_key;
		cout<<"계속하시려면 아무 숫자나 누르세요."<<endl;
		cin >> end_key;


		/*myfile.open("level1.dat");
	
		for(j=level1.mesh_y_num; j>0; j--)
		{
			for(i=1; i<level1.mesh_x_num + 1; i++)
			{
				myfile << level1.level_phi_2d[i][j] << " ";
			}
			myfile << endl;
		}

		myfile.close();



		level1.extrapolate_burning_rate();

		for(i=0; i<200; i++)
		{
			level1.propagate_surface();
			cout << "level set is propagated " << i+1 << " times" << endl;
		}

		level1.reconstruct_surface();
		cout << "surface is reconstructed" << endl;
		level1.construct_phi();
		cout << "level set is constructed" << endl;


		myfile.open("level2.dat");
	
		for(j=level1.mesh_y_num; j>0; j--)
		{
			for(i=1; i<level1.mesh_x_num + 1; i++)
			{
				myfile << level1.level_phi_2d[i][j] << " ";
			}
			myfile << endl;
		}

		myfile.close();*/
	}

	if(dimension==3)
	{
		surface_point_num = 4;
		struct surface surface1(dimension, surface_point_num, surface_point_num, false);

		surface1.dimension = dimension;
		surface1.point_num = surface_point_num;
		surface1.surf_num = surface_point_num;

		/*for(i=0; i<surface_point_num;i++)
		{
			surface1.point[i][0] = 50 + 10*cos(2*PI*i/(double)surface_point_num);
			surface1.point[i][1] = 50 + 30*sin(2*PI*i/(double)surface_point_num);
		}

		for(i=0; i<surface_point_num; i++)
		{
			surface1.edge[i][0] = i;
			surface1.edge[i][1] = (i+1)%surface_point_num;
		}*/

		surface1.point[0][0] = 0.3;
		surface1.point[0][1] = 0.3;
		surface1.point[0][2] = 0.3;
		
		surface1.point[1][0] = 0.7;
		surface1.point[1][1] = 0.3;
		surface1.point[1][2] = 0.3;
		
		surface1.point[2][0] = 0.3;
		surface1.point[2][1] = 0.7;
		surface1.point[2][2] = 0.3;
		
		surface1.point[3][0] = 0.3;
		surface1.point[3][1] = 0.3;
		surface1.point[3][2] = 0.7;

		surface1.surf[0][0] = 0;
		surface1.surf[0][1] = 2;
		surface1.surf[0][2] = 1;
		
		surface1.surf[1][0] = 0;
		surface1.surf[1][1] = 1;
		surface1.surf[1][2] = 3;
		
		surface1.surf[2][0] = 0;
		surface1.surf[2][1] = 3;
		surface1.surf[2][2] = 2;
		
		surface1.surf[3][0] = 1;
		surface1.surf[3][1] = 2;
		surface1.surf[3][2] = 3;




		int max_dim;

		cout << "maximum dimension is : ";
		cin >> max_dim;
		cout << endl;

		double time_step = pow(2,-max_dim) / 2.;

		class level_set level1(dimension, max_dim, time_step, &surface1, NULL);

		cout << "level set is constructed" << endl;
		


		string myString1 = "tree_3dlevel";
		string myString2 = ".dat";
		string myString3 = "_surface";


		level1.savingtree();


		ofstream myfile;

		myfile.open("level_3d_initial.dat");

		int num = (int)pow(2,max_dim);
		int maxnum = (num+1)*(num+1)*(num+1);

		for(i=0; i<maxnum; i++)
		{
			myfile << level1.array_savingtree[i] << " ";
		}

		myfile.close();
		

		size_t currentSize = getPeakRSS();
		size_t peakSize = getCurrentRSS();
		cout << currentSize << "   " << peakSize << endl;


		for(i=0; i<=100; i++)
		{
			if(i%5==0)
			{
				level1.savingtree();
				myfile.open(myString1+int2str(i/5)+myString2);

				/*for(j=0; j<level1.numnode; j++)
				{
					myfile << level1.mat_savingtree[j][0] << " " << level1.mat_savingtree[j][1] << " " << level1.mat_savingtree[j][2] << " " << level1.mat_savingtree[j][3];
					myfile << endl;
				}*/

				num = (int)pow(2,max_dim);
				maxnum = (num+1)*(num+1)*(num+1);

				for(j=0; j<maxnum; j++)
				{
					myfile << level1.array_savingtree[j] << " ";
				}

				myfile.close();
			}
			level1.tree_propagate_surface_velocity(i);
			cout << "level set is propagated " << i+1 << " times" << endl;

			currentSize = getPeakRSS();
			peakSize = getCurrentRSS();

			cout << currentSize << "   " << peakSize << endl;
		}


		level1.savingtree();


		myfile.open("level_3d_final.dat");

		num = (int)pow(2,max_dim);
		maxnum = (num+1)*(num+1)*(num+1);

		for(i=0; i<maxnum; i++)
		{
			myfile << level1.array_savingtree[i] << " ";
		}


		/*for(i=1; i<level1.mesh_x_num + 1; i++)
		{
			for(j=1; j<level1.mesh_y_num + 1; j++)
			{
				for(k=1; k<level1.mesh_z_num + 1; k++)
				{
					myfile << level1.level_phi_3d[i][j][k] << " ";
				}
			}
		}*/

		myfile.close();

	}

	return 0;
}