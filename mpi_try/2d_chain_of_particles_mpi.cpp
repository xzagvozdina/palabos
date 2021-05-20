/* This file is part of the Palabos library.
 *
 * The Palabos softare is developed since 2011 by FlowKit-Numeca Group Sarl
 * (Switzerland) and the University of Geneva (Switzerland), which jointly
 * own the IP rights for most of the code base. Since October 2019, the
 * Palabos project is maintained by the University of Geneva and accepts
 * source code contributions from the community.
 * 
 * Contact:
 * Jonas Latt
 * Computer Science Department
 * University of Geneva
 * 7 Route de Drize
 * 1227 Carouge, Switzerland
 * jonas.latt@unige.ch
 *
 * The most recent release of Palabos can be downloaded at 
 * <https://palabos.unige.ch/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "palabos2D.h"
#include "palabos2D.hh"
#include "offLattice/immersedWalls2D.h"
#include "offLattice/immersedWalls2D.hh"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>


using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
#define DESCRIPTOR D2Q9Descriptor

T pi = std::acos((T) -1);

plint xDirection  = 0;
plint yDirection  = 1;

plint n0 = 10; // Reference resolution. dx in particle diameter
plint ly = 8; // Resolution (cylinder diameter).
plint lx = ly * 250;
plint a = n0; // inter-particle distance to add

plint nx = lx * n0, ny = ly * n0;

T length = (T)nx; // Cylinder length (another good choice: 4*nx).


T radius = ny / 2.;

T k = 1; // repulsive forces params
T delta = 1; // repulsive forces params

T uAverage = 0.06;
T Reynolds = 1;

plint maxIter  = 4000001;//100000; // Maximum number of iterations for the simulation.
plint outIter_vtk  = 20e2;//maxIter+1; // Number of iterations for printing the average kinetic energy on the screen.
plint outIter_particles = 1; // Number of iterations for saving data in the disk.

/// A functional, used to create an initial condition for with zero velocity,
///   and linearly decreasing pressure.
template<typename T>
class Poiseuille {
public:
    Poiseuille(T _uAve, T _rho, T _rad, T _cy)
            : uAve(_uAve), rho(_rho), rad(_rad), cy(_cy)
    {
    }
    void operator()(plint iX, plint iY,  T& rho_, Array<T,2>& u) const {
        //rho = poiseuilleDensity(iX,parameters);
        T r = std::fabs(iY - cy)/rad;
        rho_ = rho;
        if ( r <= 1.){
            u[0] = 2 * uAve * (1 - r * r);
        }
        else
            u[0] = T();
        u[1] = T();
    }
private:
    T uAve; T rho;	T rad;	T cy;
};

template<typename T>
class PoiseuilleProfile {
public:
    PoiseuilleProfile(T _uAve, T _rad, T _cy)
            : uAve(_uAve), rad(_rad), cy(_cy)
    {
    }
    void operator()(plint iX, plint iY, Array<T,2>& u)const {
        T r = std::fabs(iY - cy) / rad;
        if ( r <= 1.){
            u[0] = 2 * uAve * (1 - r * r);
        }
        else
            u[0] = T();
        u[1] = T();
    }
private:
    T uAve; T rad;	T cy;
};

/// A functional, used to initialize a pressure boundary to constant density
template<typename T>
class ConstantDensity {
public:
    ConstantDensity(T density_)
            : density(density_)
    { }
    T operator()(plint iX, plint iY) const {
        return density;
    }
private:
    T density;
};

class SurfaceVelocity {
public:
    SurfaceVelocity(T dt_, T m_, T velLB_, T cx_, T cy_, T rad_, T angvel_=0.)
    {
        dt = dt_;
        m = m_;
        //velocity = {inletvelLB_, 0.};
        center_velocity = {velLB_, 0.};
        center = {cx_, cy_};
        omega = angvel_;
        rad = rad_;
        rad2 = rad_ * rad_;
    }

    Array<T,2> operator()(Array<T,2> const& pos)
    {
        Array<T, 2> r = pos - center;
		return center_velocity + Array<T, 2>({-r[1] * omega, r[0] * omega});
    }

    void update(const Array<T, 2>& g, const Array<T, 3>& torque){
		omega += 2 * torque[2] / m / rad2;
        center_velocity += g / m;
    }

    void updateRepulsive(const Array<T, 2>& f){
        center_velocity += f / m;
    }

    void moveCenter(){
		
        center += center_velocity;
    }
	
	Array<T, 2> getCenterVelocity(){
        return center_velocity;
    }

    Array<T, 2> getCenter(){
        return center;
    }
	
	T getOmega(){
		return omega;
	}

private:
    //Array<T, 2> velocity;
    Array<T, 2> center;
    Array<T, 2> center_velocity;
    T omega;
    T dt;
    T m;
    T rad;
    T rad2;
};


void cylinderSetup( MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
                    IncomprFlowParam<T> const& parameters,
                    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>& boundaryCondition )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();

    pcout << "nx " << nx << " ny " << ny << endl;

    // Create Velocity boundary conditions everywhere
    boundaryCondition.setVelocityConditionOnBlockBoundaries (
            lattice, Box2D(0, 0, 1, ny-2) );

    defineDynamics(lattice,  Box2D(0, nx-1, ny-1, ny-1),
            new BounceBack<T, DESCRIPTOR>(1.));

    defineDynamics(lattice,  Box2D(0, nx-1, 0, 0),
                   new BounceBack<T, DESCRIPTOR>(1.));
      
    //boundaryCondition.addVelocityBoundary0N(Box2D(0, 0, 1, ny-2), lattice);               
    //setBoundaryVelocity(lattice, Box2D(0, 0, 1, ny-2), Array<T,2>(uAverage, (T) 0));               

    //integrateProcessingFunctional(new CopyUnknownPopulationsFunctional2D<T,DESCRIPTOR, 0, +1>,
    //                              Box2D(nx-1,nx-1, 1,ny-2), lattice);
                                  
    boundaryCondition.setVelocityConditionOnBlockBoundaries (
            lattice, Box2D(nx-1,nx-1, 1,ny-2) );


    initializeAtEquilibrium (
            lattice, lattice.getBoundingBox(), 
            //1., Array<T,2>(uAverage, (T) 0));
            Poiseuille<T>(uAverage, 1., radius, (T) ny / 2. ) );

    //lattice.initialize();
}

void writeGif(MultiBlockLattice2D<T,DESCRIPTOR>& lattice, plint iter)
{
    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif(createFileName("u", iter, 6),
                               *computeVelocityNorm(lattice) );
}

void writeVTK(MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
              IncomprFlowParam<T> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    //T dt = parameters.getDeltaT();
    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, 10), dx);
    vtkOut.writeData<2,float>(*computeVelocity(lattice), "velocity");
}

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);
    char outputDir[] = "./tmp_mpi/";
    global::directories().setOutputDir(outputDir);

    IncomprFlowParam<T> parameters(
            (T) uAverage,//u_drop
            //(T) uAverage/5,  // uMax
            (T) Reynolds,    // Re
            n0,        // N
            lx,        // lx
            ly         // ly
    );

    plint ibIter = 4;

    writeLogFile(parameters, "Chain of drops. cetral is higher than others");

    bool incompressibleModel = true;
    //MultiBlockLattice2D<T, DESCRIPTOR> lattice (
    //        parameters.getNx(), parameters.getNy(),
    //       new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()) );
    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
            parameters.getNx(), parameters.getNy(),
            new IncBGKdynamics<T,DESCRIPTOR>(parameters.getOmega()) );

    MultiScalarField2D<T> *rhoBar = generateMultiScalarField<T>((MultiBlock2D&) lattice, (plint)4).release();
    rhoBar->toggleInternalStatistics(false);

    MultiTensorField2D<T,2> *j = generateMultiTensorField<T,2>((MultiBlock2D&) lattice, (plint)4).release();
    j->toggleInternalStatistics(false);

    std::vector<MultiBlock2D*> rhoBarJarg;
    rhoBarJarg.push_back(&lattice);
    rhoBarJarg.push_back(rhoBar);
    rhoBarJarg.push_back(j);

    lattice.periodicity().toggleAll(false);
    rhoBar->periodicity().toggleAll(false);
    j->periodicity().toggleAll(false);

    integrateProcessingFunctional(
            new ExternalRhoJcollideAndStream2D<T,DESCRIPTOR>(),
            lattice.getBoundingBox(), rhoBarJarg, 0);
    integrateProcessingFunctional(
            new BoxRhoBarJfunctional2D<T,DESCRIPTOR>(),
            lattice.getBoundingBox(), rhoBarJarg, 2);

    // The next container block is necessary for the immersed-wall algorithm.
    MultiContainerBlock2D container(*rhoBar);

    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
            boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();

    cylinderSetup(lattice, parameters, *boundaryCondition);
    
    
	pcout << "Re " << Reynolds << endl;
    pcout << "nu " << parameters.getLatticeNu() << endl;
    pcout << "omega " << parameters.getOmega() << endl;
    pcout << "tau " << parameters.getTau() << endl;
    pcout << "dx " << parameters.getDeltaX() << endl;
    pcout << "dx^2 " << parameters.getDeltaX()*parameters.getDeltaX() << endl;
    pcout << "dt " << parameters.getDeltaT() << endl;
    pcout << "ny " << ny << endl;
    pcout << "Ny " << parameters.getNy() << endl;
	
    Box2D inlet(0,0, 1, parameters.getNy()-2);
    Box2D outlet(parameters.getNx()-1,parameters.getNx()-1, 1, parameters.getNy()-2);

	plb::global::MpiManager& mpiman = plb::global::mpi();
    //T sigma = period / (T)10;
    //T m, signal, uAve;
    T cy = parameters.getNy() / (T)2;
    T dx = parameters.getDeltaX();

    pcout << "Creating the immersed sphere surface." << std::endl;


   std::vector<Array<T,4> > centers_rad;

    T xc, yc, r;
    double count_particles = 0.;
    // cx, cy, rad    
    double rad = n0/2;
    
    T cx0 = 250.;//(T)2 * rad;
    T cy0 = (T)ny / 2.;
    
    int num_particles = 8;
    
    std::mt19937 gen(12345);
    std::uniform_real_distribution<> dist_float(-0.5, 0.5);
      
    
    for(int i = 0; i < num_particles; i++, count_particles++){
    	centers_rad.push_back({cx0 + i * (a + 2 * rad) + dist_float(gen), cy0 + dist_float(gen), 
    							rad, count_particles});
    }
    
    //centers_rad.push_back({cx0, cy0, rad, count_particles});
    //count_particles = 1.;

    pcout << "number of particles "<< num_particles <<std::endl;
    //pcout << "cx " << n0 / 2 + 0.5 * flag_asym << " cy = " << n0 / 2 + 0.5 * flag_asym << std::endl;
    //return 0;

    vector<vector<Array<T,2>>> all_vertices;
    vector<vector<T>> all_areas;
    vector<SurfaceVelocity> all_velocities;
    T n_vert = 100;
    T rho = 1000./970.;
    T mass_s;

    for (int ind = 0; ind < num_particles; ind++){
        vector<Array<T,2> > vertices;
        vector<T> areas;
        xc = centers_rad[ind][0];
        yc = centers_rad[ind][1];
        r = centers_rad[ind][2];

        n_vert = n_vert;//(2 * pi * particle_r) * 3;//(r == 4) ? 50 : 150;
        pcout << "ind " << ind <<" n_vert " << n_vert << "cx = " << xc 
        	<< " cy = " << yc << std::endl;

        mass_s = pi * r * r * rho;
        for (pluint iVertex = 0; iVertex < (pluint) n_vert; iVertex++) {
            Array<T, 2> vert =  {xc + r * cos(iVertex * 2 * pi / n_vert),
                                 yc + r * sin(iVertex * 2 * pi / n_vert)};
            vertices.push_back(vert);
            areas.push_back(2 * 2 * pi * r / n_vert);
        }

        all_vertices.push_back(vertices);
        all_areas.push_back(areas);

        //T rad = fabs(cy-yc) / cy;
        SurfaceVelocity vel = SurfaceVelocity( parameters.getDeltaT(), mass_s, 2*uAverage, xc, yc, r);
        //SurfaceVelocity vel = SurfaceVelocity( parameters.getDeltaT(), mass_s, 0, xc, yc, r);
        all_velocities.push_back(vel);
    }
    
    
	//mpiman.init();		//mpiman.init();
	
	int mpirank = mpiman.getRank();	
	int mpisize = mpiman.getSize();	
	int particles_per_process = num_particles / mpisize + int(num_particles%mpisize > 0);
	
	pcout << "mpiman size " << mpisize << endl;	
	cout << "mpiman rank " << mpirank << endl;	
	pcout << "particles_per_process " << particles_per_process << endl;	
    
    char buffer[256];

    T previousIterationTime = T();
    // Main loop over time iterations.
    global::timer("mainLoop").restart();
    for (plint iT=0; iT < 200; ++iT) {
        global::timer("mainLoop").restart();

        setBoundaryVelocity(lattice, inlet, PoiseuilleProfile<T>(uAverage, radius, cy));
        setBoundaryVelocity(lattice, outlet, PoiseuilleProfile<T>(uAverage, radius, cy));

        // Lattice Boltzmann iteration step.
        lattice.executeInternalProcessors();
		
		//if (iT > 14e5){
		if (true){
        plint iVertex;
        
	    //for (unsigned int ind = 0; ind < all_vertices.size(); ind++) {
	    for (unsigned int ind = mpirank*particles_per_process; ind < (mpirank+1)*particles_per_process; ind++) {
	    	if (ind >= num_particles) continue;
	    	
	        all_velocities[ind].moveCenter();
	        for (iVertex = 0; iVertex < (plint) all_vertices[ind].size(); iVertex++) {
	            all_vertices[ind][iVertex] += all_velocities[ind](all_vertices[ind][iVertex]); 
	        }
	        
	        // particles floated from channel 
	        if (all_vertices[ind][0][0] > parameters.getNx() + 2 * centers_rad[ind][2])  {
	        	//cout << all_vertices[ind][0][0] << endl;
	        	int mpirank2 = mpiman.getRank();	
				int mpisize2 = mpiman.getSize();	
				mpiman.abort(13);				
				cout << "1 mpiman size " << mpisize2 << endl;	
				cout << "1 mpiman rank " << mpirank2 << endl;	
				mpiman.abort(13);
	        	return 0;
	        	
	        }
	        // particle crosses wall
	        Array<T, 2> center= all_velocities[ind].getCenter();
	        if (center[1] > ny  || center[1] < 0){
	        	pcout << iT << " particle " << ind << "crossed walls, y = " << center[1] << endl;
	        }   
	    }


	    Array<T, 2> ff; Array<T, 3> torque;

	    Array<T, 2> f12, dist;
	    Array<T, 2> center1, center2;
	    
	    int env = 2;
	    int nx = 2 * rad + 2 * env;
	    int ny = 2 * rad + 2 * env;
	    
	    int scount = 4 + nx * ny * 3; // 1 rhoBar 2 j[0] 3 j[1]
        int rcount = scount * mpisize;
    	double* rbuf = new double[rcount];
    	//cout << "nx=" << nx << " ny=" << ny << " scount=" <<scount << " rcount=" <<rcount <<endl;
	    
		for (unsigned int ind = mpirank*particles_per_process; ind < (mpirank+1)*particles_per_process; ind++) {
			//if (ind >= num_particles) continue;
	    		
	    	center1 = all_velocities[ind].getCenter();
	    	
	    	Box2D envelope((int)center1[0]-rad-env, (int)center1[0]-rad+env, (int)center1[1]-rad-env, (int)center1[1]+rad+env);
			
			static std::vector< Array<T,2> > dummyNormals;
			std::vector<MultiBlock2D*> args;
			args.push_back(&container);
			applyProcessingFunctional (
            	new InstantiateImmersedWallData2D<T>(all_vertices[ind], all_areas[ind], dummyNormals), envelope, args );
			
	        //instantiateImmersedWallData(all_vertices[ind], all_areas[ind], container);
	        for (int i = 0; i < ibIter; i++) {
			    args.resize(0);
			    
			    //std::vector<MultiBlock2D*> args;
				args.push_back(rhoBar);
				args.push_back(j);
				args.push_back(&container);
				applyProcessingFunctional (
					new InamuroIteration2D<T,SurfaceVelocity>(all_velocities[ind], (T) 1.0 / parameters.getOmega(), incompressibleModel), envelope, args);
				
	        }
	        
	        //mpimap.allToAll_double(sbuf, scount, rbuf, rcount);
	        /*int sx = (int)center1[0]-rad-env, ex = sx + nx;
	    	int sy = (int)center1[1]-rad-env, ey = sy + ny; 
			
	    	double* arr = rbuf + scount * mpirank;
	        arr[0] = sx; arr[1] = ex; arr[2] = sy; arr[3] = ey;
	    	// rhoBar broadcasting to other processes
	    	for (int ix = 0; ix < nx; ix++){
	    		for(int iy = 0; iy < ny; iy++){
	    			arr[4 + ix * ny + iy] = rhoBar->get(sx+ix, sy+iy);
	    			arr[4 + ix * ny + iy + nx*ny] = j->get(sx+ix, sy+iy)[0];
	    			arr[4 + ix * ny + iy + 2*nx*ny] = j->get(sx+ix, sy+iy)[1];
	    			if(ix == 17 && iy == 17){
	    			cout << mpirank << "  " << ix << "  "<<iy <<  "  "<< arr[4 + ix * ny + iy]<<
	    			"  " << arr[4 + ix * ny + iy + nx*ny] << 
	    			"  " << arr[4 + ix * ny + iy + 2*nx*ny] << "  " <<endl;}
	    		}
	    	}
	    	mpiman.barrier();
	    	for (int r = 0; r < mpisize; r++){
	        	mpiman.bCast(rbuf + scount*r, scount, r);
	        }
	        
	        
	        mpiman.barrier();
	        
	        for (int r = 0; r < mpisize; r++){
	        	arr = rbuf + scount * r;
	        	
	        	cout << mpirank << " from " << r <<"  " << arr[0] << "  " << arr[1] << "  "
	        	<< arr[2] << "  " << arr[3] << endl;
	        	
	        	for (int ix = 0; ix < arr[1]-arr[0]; ix++){
	    			for(int iy = 0; iy < arr[3]-arr[2]; iy++){
	    				rhoBar->get(arr[0]+ix, arr[2]+iy) = arr[4+ix*(int)(arr[3]-arr[2])+iy];
	    				j->get(arr[0]+ix, arr[2]+iy)[0] = arr[4+ix*(int)(arr[3]-arr[2])+iy + nx*ny];
	    				j->get(arr[0]+ix, arr[2]+iy)[1] = arr[4+ix*(int)(arr[3]-arr[2])+iy + 2*nx*ny];
	    				 
	    				//cout << arr[4+ix*(int)(arr[3]-arr[2])+iy] << "  ";
	    				if(false){
	    				//if (ix == 17 && iy == 17){
	    				cout << mpirank << " after " <<r<<"  " << ix << "  "<<iy <<  "  "<< rhoBar->get(arr[0]+ix, arr[2]+iy)<<
	    			"  " << j->get(arr[0]+ix, arr[2]+iy)[0] << 
	    			"  " << j->get(arr[0]+ix, arr[2]+iy)[1] << "  " <<endl;}
	    			}
	    		}
	        }*/
	        mpiman.barrier();
	        //pcout << j->get(40,100)[0] << " " << j->get(40,100)[1] << endl;
	        args.resize(0);
	        args.push_back(&container);
			ReduceImmersedTorque2D<T> functional_torque(center1, 0);
			applyProcessingFunctional(functional_torque, envelope, args);
			torque = functional_torque.getSumTorque();
			
			ReduceImmersedForce2D<T> functional_force(0);
			applyProcessingFunctional(functional_force, envelope, args);
			ff = functional_force.getSumG();
	        //ff = reduceImmersedForce<T>(container, 0);

	        //torque = reduceImmersedTorque(container, center1, 0);
	        
			all_velocities[ind].update(-ff, -torque);
			
			
            
			if (iT % 1 == 0) {
	            sprintf(buffer, "%ssphere_coords_8_%05d.txt", outputDir, (int)centers_rad[ind][3]);
	            ofstream file_coord(buffer, std::ofstream::out | std::ofstream::app);

				Array<T,2> vel_ = all_velocities[ind].getCenterVelocity();			   
	            file_coord << std::setprecision(12) << iT << " " 
							<< center1[0] << " "
							<< center1[1] << " "
							<< all_vertices[ind][0][0] << " " 
	                        << all_vertices[ind][0][1] << " "
							<< vel_[0]  << " "
							<< vel_[1]  << " "
							<< ff[0] << " "
							<< ff[1] << " "
							<< torque[2] << " "
							<< all_velocities[ind].getOmega() << std::endl;
	            file_coord.close();
	        }
	    }
        mpiman.barrier();
        delete rbuf;
        }
        
    }
    previousIterationTime = global::timer("mainLoop").stop();
    pcout << "previousIterationTime " << previousIterationTime << std::endl;


    delete boundaryCondition;
}
