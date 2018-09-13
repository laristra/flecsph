#include <math.h>
#include <iostream>
#include <vector>
#include <iterator>
#include <cassert>
#include <hdf5.h>
#include <algorithm>

/**
 * Binary Neutron Star Initial Data generator 
 */
const int64_t nneighbs = 6;
const double solarMass = 1.9885 * pow(10,30) * pow(10,3);
const double gravConst = 1;
const double ksi_1 = M_PI;
const double radiusKM = pow(10,-5);
const double mass = 1.; 
const double radiusCM = radiusKM*pow(10,5);
const double K_ratio = 1./(radiusCM*radiusCM*radiusCM*gravConst);
const double Kconstant = (radiusCM*radiusCM*2*M_PI*gravConst)/(ksi_1*ksi_1) * K_ratio;
const double Aconstant = sqrt(2.*M_PI*gravConst/Kconstant);
const double VolSphere = 4.*M_PI*pow(radiusCM,3.)/3.;
const double centralDensity = (mass*Aconstant*Aconstant*Aconstant)/
	(4.*M_PI*(sin(Aconstant*radiusCM)-Aconstant*radiusCM*cos(Aconstant*radiusCM)));
const double tRelax = pow(gravConst*centralDensity,-1./2.);
const double volPart = (4./3.*M_PI*pow(radiusCM,3.));
const double dist_stars = 2.9;

double H_ratio = 1./radiusCM;
double Hconstant = 0.;
int64_t nparticles = 0;
double QZZ; 
double angularMomentum; 
double Omega; 


typedef struct particle{
	double x_, y_, z_; 
	double mass_; 
	double density_;

	particle(double x, double y, double z){
		x_ = x; y_ = y; z_ = z; 
		mass_ = 0.;
		density_ = 0.;
	};
}particle; 

double 
radiusSphere(double x, double y , double z){
	return sqrt(x*x+y*y+z*z);
}

void 
generate_position(std::vector<particle>& v, int64_t ntarget,double part_dist){
	// First generate particles inside sphere of radius 1 
	for(int64_t i = 0 ; i < ntarget; ++i)
	{
		double x = -1 + i * part_dist;
		for(int64_t j = 0; j < ntarget; ++j)
		{
			double y = -1 + j * part_dist;
			for(int64_t k = 0 ; k < ntarget; ++k)
			{
				double z = -1 + k * part_dist;
				double distance = radiusSphere(x,y,z);
				if(distance <= 1.0)
				{
					v.push_back(particle(x,y,z));
				}
			}
		}
	}
}

void 
compute_smoothing_length(int64_t nparticles)
{
	Hconstant = 1./2.*radiusCM*pow(nneighbs*1./nparticles,1./3.) * H_ratio;
}

void 
printConst()
{
	std::cout<<std::endl;
	std::cout<<"Kconstant = "<<Kconstant<<std::endl;
	std::cout<<"radiusCM  = "<< radiusCM<<std::endl;
	std::cout<<"Aconstant = "<<Aconstant<<std::endl;
	std::cout<<"VolSphere = "<<VolSphere<<std::endl;
	std::cout<<"Hconstant = "<<Hconstant<<std::endl;
	std::cout<<"CentDens. = "<<centralDensity<<std::endl;
	std::cout<<"trelax.   = "<<tRelax<<std::endl;
	std::cout<<"volPart   = "<<volPart<<std::endl;
}

// Lane emden equation stuff 
double sinc(double ksi)
{
	return sin(ksi)/ksi;
}

double density(double radius)
{
	return sinc(radius*Aconstant)*centralDensity;
}

// Mass of particle from its density 
double massDensity(double density)
{
	double n = nparticles/volPart;
	return density/n;
}

void 
compute_particles_mass(std::vector<particle>& v)
{
	double totalMass = 0.;
	// Compute the density and then mass for the particles 
	for(int i = 0; i < nparticles; ++i)
	{
		if(v[i].x_ == 0. &&  v[i].y_ == 0. && v[i].z_ == 0.)
		{
			v[i].density_ = centralDensity;
		}else
		{
			v[i].density_ = density(
				radiusSphere(v[i].x_,v[i].y_,v[i].z_));
		}
		// Then associate mass and other stuff and normalize
		v[i].mass_ = massDensity(v[i].density_) / mass;
		totalMass += v[i].mass_;
	}
	std::cout<<std::endl;
	std::cout<<"Total mass = "<<totalMass<<std::endl;
}


void 
write_dataset(
	hid_t file, 
	const char * name,
	std::vector<double> data)
{
	hsize_t size = data.size();
	hid_t space_id = H5Screate_simple(1,&size,NULL);
	hid_t dataset = H5Dcreate(
		file, name, H5T_IEEE_F64LE,space_id,
		H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

	herr_t status = H5Dwrite(dataset,H5T_IEEE_F64LE,
		H5S_ALL,H5S_ALL,H5P_DEFAULT,&data[0]);
	status = H5Sclose (space_id);
    status = H5Dclose (dataset);
}

void 
write_attribute(
	hid_t file, 
	const char * name,
	int value)
{
	hsize_t size = 1;
	hid_t space_id = H5Screate_simple(1,&size,NULL);
	hid_t attribute = H5Acreate(
		file, name, H5T_NATIVE_INT,space_id,
		H5P_DEFAULT,H5P_DEFAULT);
	herr_t status = H5Awrite(attribute,H5T_NATIVE_INT,&value);
	status = H5Sclose (space_id);
    status = H5Aclose (attribute);
}

void 
write_attribute(
	hid_t file, 
	const char * name,
	double value)
{
	hsize_t size = 1;
	hid_t space_id = H5Screate_simple(1,&size,NULL);
	hid_t attribute = H5Acreate(
		file, name, H5T_IEEE_F64LE,space_id,
		H5P_DEFAULT,H5P_DEFAULT);
	herr_t status = H5Awrite(attribute,H5T_IEEE_F64LE,&value);
	status = H5Sclose (space_id);
    status = H5Aclose (attribute);
}

void 
compute_angular_momentum(
	std::vector<particle>& v)
{
	QZZ = 0.;
	for(auto& p: v){
		QZZ += p.mass_*(p.x_*p.x_+p.y_*p.y_);
	}
	std::cout<<"QZZ   = "<<QZZ<<std::endl;
	Omega = sqrt(gravConst*2/(dist_stars*dist_stars*dist_stars));
	std::cout<<"Omega = "<<Omega<<std::endl;
	angularMomentum = QZZ*Omega;
	std::cout<<"angularMoment = "<<angularMomentum<<std::endl;
}

int main(int argc, char* argv[])
{
	int64_t ntarget = 10;
	if(argc == 2)
	{
		ntarget = atoi(argv[1]);
	}
	
	double part_dist = 2./(ntarget-1);
	
	std::cout<<"ntarget = "<<ntarget<<" dist = "<<part_dist<<std::endl;
	std::vector<particle> particles; 
	generate_position(particles,ntarget,part_dist);
	nparticles = particles.size();
	std::cout<<"Position generated = "<<particles.size()<<std::endl;
	// Compute the smoothing length 
	compute_smoothing_length(nparticles);
	printConst();
	compute_particles_mass(particles);

	// Suppress elements with mass less than 1.e-10
	particles.erase(remove_if(
			particles.begin(),particles.end(),
			[](const auto& e) -> bool
			{ 
				return e.mass_ < 1.0e-10; 
			})
		,particles.end());

	nparticles = particles.size();
	std::cout<<std::endl<<"Particles after cleaning = "<<nparticles<<std::endl;

	std::cout<<"Generating two stars = "<<nparticles*2<<std::endl;

	compute_angular_momentum(particles);


	char name[128]; 
	sprintf(name,"hdf5_bns_3D_%d.h5part",nparticles*2);

	auto dataFile = H5Fcreate(
		name,
		H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

	auto step_id = H5Gcreate(dataFile, "/Step#0", 
	 	H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);



	write_attribute(dataFile,"dimension",3);
	write_attribute(dataFile,"angularMomentum",angularMomentum);


	// Output in HDF5 format 
	std::vector<double> data1(nparticles*2); 
	std::vector<double> data2(nparticles*2); 
	std::vector<double> data3(nparticles*2);

	// Positions 1st star 
	for(int64_t i = 0 ; i < 2*nparticles; ++i){
		data1[i] = particles[i%nparticles].x_ + dist_stars/2.*pow(-1,i/nparticles); 
		data2[i] = particles[i%nparticles].y_; 
		data3[i] = particles[i%nparticles].z_; 
	}

	write_dataset(dataFile, "/Step#0/x",data1);
	write_dataset(dataFile, "/Step#0/y",data2);
	write_dataset(dataFile, "/Step#0/z",data3);

	std::cout<<"X    ["<<*std::min_element(data1.begin(),data1.end())<<","
	<< *std::max_element(data1.begin(),data1.end()) <<"]"<<std::endl;
	
	std::cout<<"Y    ["<<*std::min_element(data2.begin(),data2.end())<<","
	<< *std::max_element(data2.begin(),data2.end()) <<"]"<<std::endl;
	
	std::cout<<"Z    ["<<*std::min_element(data3.begin(),data3.end())<<","
	<< *std::max_element(data3.begin(),data3.end()) <<"]"<<std::endl;

	// Mass 1st star 
	for(int64_t i = 0 ; i < nparticles*2; ++i){
		data1[i] = particles[i%nparticles].mass_; 
		data2[i] = Hconstant;
		data3[i] = particles[i%nparticles].density_;
	}

	write_dataset(dataFile, "/Step#0/m",data1);
	write_dataset(dataFile, "/Step#0/h",data2);
	write_dataset(dataFile, "/Step#0/rho",data3);

	std::cout<<"mass ["<<*std::min_element(data1.begin(),data1.end())<<","
	<< *std::max_element(data1.begin(),data1.end()) <<"]"<<std::endl;

	std::cout<<"H    ["<<*std::min_element(data2.begin(),data2.end())<<","
	<< *std::max_element(data2.begin(),data2.end()) <<"]"<<std::endl;

	std::cout<<"dens ["<<*std::min_element(data3.begin(),data3.end())<<","
	<< *std::max_element(data3.begin(),data3.end()) <<"]"<<std::endl;

	// Empty data sets 
	std::fill(data1.begin(),data1.end(),0);
	write_dataset(dataFile, "/Step#0/ax",data1);
	write_dataset(dataFile, "/Step#0/ay",data1);
	write_dataset(dataFile, "/Step#0/az",data1);
	write_dataset(dataFile, "/Step#0/vx",data1);
	write_dataset(dataFile, "/Step#0/vy",data1);
	write_dataset(dataFile, "/Step#0/vz",data1);
	write_dataset(dataFile, "/Step#0/P",data1);

	H5Gclose(step_id);
	H5Fclose(dataFile);


	return EXIT_SUCCESS; 
}
