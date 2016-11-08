// Make sure the header is only included once
#ifndef SPHERE_H_
#define SPHERE_H_

#include <vector>

//#include "eigen-nvcc/Eigen/Core"

#include "datatypes.h"

// DEM class
class DEM {

    // Values and functions only accessible from the class internally
    private:

        // Input filename (full path)
        std::string inputbin;

        // Simulation ID
        std::string sid;

        // Output level
        int verbose;

        // Number of dimensions
        unsigned int nd;

        // Number of particles
        unsigned int np;

        // HOST STRUCTURES
        // Structure containing individual particle kinematics
        Kinematics k;

        // Structure containing energy values
        Energies e;

        // Structure of global parameters
        Params params;

        // Structure containing spatial parameters
        Grid grid;

        // Structure of temporal parameters
        Time time;

        // Structure of wall parameters
        Walls walls;

        // Image structure (red, green, blue, alpa)
        rgba* img;
        unsigned int width;
        unsigned int height;

        // Device management
        int ndevices;     // number of CUDA GPUs
        int device;       // primary GPU
        int* domain_size; // elements per GPU


        // DEVICE ARRAYS

        // Particle kinematics arrays
        Float4        *dev_x;
        Float4        *dev_xyzsum;
        Float4        *dev_vel;
        Float4        *dev_vel0;
        Float4        *dev_acc;
        Float4        *dev_force;
        Float4        *dev_angpos;
        Float4        *dev_angvel;
        Float4        *dev_angvel0;
        Float4        *dev_angacc;
        Float4        *dev_torque;
        unsigned int  *dev_contacts;
        Float4        *dev_distmod;
        Float4        *dev_delta_t;
        Float         *dev_es_dot;
        Float         *dev_es;
        Float         *dev_ev_dot;
        Float         *dev_ev;
        Float         *dev_p;

        // Sorted kinematics arrays
        Float4        *dev_x_sorted;
        Float4        *dev_vel_sorted;
        Float4        *dev_angvel_sorted;

        // Sorting grid arrays
        unsigned int  *dev_gridParticleCellID;
        unsigned int  *dev_gridParticleIndex;
        unsigned int  *dev_cellStart;
        unsigned int  *dev_cellEnd;

        // Wall arrays
        int           *dev_walls_wmode;
        Float4        *dev_walls_nx;        // normal, pos.
        Float4        *dev_walls_mvfd;      // mass, velo., force, dev. stress
        Float         *dev_walls_tau_x;     // wall shear stress
        Float         *dev_walls_force_partial; // Pre-sum per wall
        Float         *dev_walls_force_pp;  // Force per particle per wall
        Float         *dev_walls_acc;       // Wall acceleration
        Float         *dev_walls_tau_eff_x_pp;      // Shear force per particle
        Float         *dev_walls_tau_eff_x_partial; // Pre-sum of shear force

        // Bond arrays
        uint2         *dev_bonds;           // Particle bond pairs
        Float4        *dev_bonds_delta;     // Particle bond displacement
        Float4        *dev_bonds_omega;     // Particle bond rotation

        // Raytracer arrays
        unsigned char *dev_img;
        float4        *dev_ray_origo;       // Ray data always single precision
        float4        *dev_ray_direction;


        // GPU initialization, must be called before startTime()
        void initializeGPU();

        // Copy all constant data to constant device memory
        void transferToConstantDeviceMemory();
        void rt_transferToConstantDeviceMemory();

        // Check for CUDA errors
        void checkForCudaErrors(const char* checkpoint_description,
                const int run_diagnostics = 1);
        void checkForCudaErrorsIter(const char* checkpoint_description,
                const unsigned int iteration,
                const int run_diagnostics = 1);

        // Check values stored in constant device memory
        void checkConstantMemory();

        // Initialize camera values and transfer to constant device memory
        void cameraInit(const float3 eye,
                const float3 lookat,
                const float imgw,
                const float focalLength);

        // Allocate global device memory to hold data
        void allocateGlobalDeviceMemory();
        void rt_allocateGlobalDeviceMemory();

        // Allocate global memory on helper devices
        void allocateHelperDeviceMemory();
        void freeHelperDeviceMemory();

        // Free dynamically allocated global device memory
        void freeGlobalDeviceMemory();
        void rt_freeGlobalDeviceMemory();

        // Copy non-constant data to global GPU memory
        void transferToGlobalDeviceMemory(int status = 1);

        // Copy non-constant data from global GPU memory to host RAM
        void transferFromGlobalDeviceMemory();
        void rt_transferFromGlobalDeviceMemory();

        // Find and return the max. radius
        Float r_max();

        // Write porosities found in porosity() to text file
        void writePorosities(
                const char *target,
                const int z_slices,
                const Float *z_pos,
                const Float *porosity);

        // Lattice-Boltzmann data arrays (D3Q19)
        Float  *f;          // Fluid distribution (f0..f18)
        Float  *f_new;      // t+deltaT fluid distribution (f0..f18)
        Float  *dev_f;      // Device equivalent
        Float  *dev_f_new;  // Device equivalent
        Float4 *v_rho;      // Fluid velocity v (xyz), and pressure rho (w)
        Float4 *dev_v_rho;  // Device equivalent

        //// Porous flow
        int fluid;      // 0: no, 1: yes
        int cfd_solver; // 0: Navier Stokes, 1: Darcy

        // Navier Stokes values, host
        NavierStokes ns;

        // Navier Stokes values, device
        Float*  dev_ns_p;            // Cell hydraulic pressure
        Float3* dev_ns_v;            // Cell fluid velocity
        Float*  dev_ns_v_x;          // Cell fluid velocity in staggered grid
        Float*  dev_ns_v_y;          // Cell fluid velocity in staggered grid
        Float*  dev_ns_v_z;          // Cell fluid velocity in staggered grid
        Float3* dev_ns_v_p;          // Averaged predicted cell fluid velocity
        Float*  dev_ns_v_p_x;        // Predicted cell fluid velocity in st. gr.
        Float*  dev_ns_v_p_y;        // Predicted cell fluid velocity in st. gr.
        Float*  dev_ns_v_p_z;        // Predicted cell fluid velocity in st. gr.
        Float3* dev_ns_vp_avg;       // Average particle velocity in cell
        Float*  dev_ns_d_avg;        // Average particle diameter in cell
        Float3* dev_ns_F_pf;         // Interaction force on fluid
        //Float*  dev_ns_F_pf_x;       // Interaction force on fluid
        //Float*  dev_ns_F_pf_y;       // Interaction force on fluid
        //Float*  dev_ns_F_pf_z;       // Interaction force on fluid
        Float*  dev_ns_phi;          // Cell porosity
        Float*  dev_ns_dphi;         // Cell porosity change
        //Float3* dev_ns_div_phi_v_v;  // Divegence used in velocity prediction
        Float*  dev_ns_epsilon;      // Pressure difference
        Float*  dev_ns_epsilon_new;  // Pressure diff. after Jacobi iteration
        Float*  dev_ns_epsilon_old;  // Pressure diff. before Jacobi iteration
        Float*  dev_ns_norm;         // Normalized residual of epsilon values
        Float*  dev_ns_f;            // Values of forcing function
        Float*  dev_ns_f1;           // Constant terms in forcing function
        Float3* dev_ns_f2;           // Constant slopes in forcing function
        Float*  dev_ns_v_prod;       // Outer product of fluid velocities
        //Float*  dev_ns_tau;          // Fluid stress tensor
        Float3* dev_ns_div_phi_vi_v; // div(phi*vi*v)
        //Float3* dev_ns_div_phi_tau;  // div(phi*tau)
        //Float3* dev_ns_div_tau;      // div(tau)
        Float*  dev_ns_div_tau_x;    // div(tau) on x-face
        Float*  dev_ns_div_tau_y;    // div(tau) on y-face
        Float*  dev_ns_div_tau_z;    // div(tau) on z-face
        Float3* dev_ns_f_pf;         // Interaction force on particles
        Float4* dev_ns_f_d;          // Drag force on particles
        Float4* dev_ns_f_p;          // Pressure force on particles
        Float4* dev_ns_f_v;          // Viscous force on particles
        Float4* dev_ns_f_sum;        // Total force on particles

        // Helper device arrays, input
        unsigned int** hdev_gridParticleIndex;
        unsigned int** hdev_gridCellStart;
        unsigned int** hdev_gridCellEnd;
        Float4** hdev_x;
        Float4** hdev_x_sorted;
        Float4** hdev_vel;
        Float4** hdev_vel_sorted;
        Float4** hdev_angvel;
        Float4** hdev_angvel_sorted;
        Float4** hdev_walls_nx;
        Float4** hdev_walls_mvfd;
        Float4** hdev_distmod;

        // Helper device arrays, output
        Float4** hdev_force;
        Float4** hdev_torque;
        Float4** hdev_delta_t;
        Float** hdev_es_dot;
        Float** hdev_es;
        Float** hdev_ev_dot;
        Float** hdev_ev;
        Float** hdev_p;
        Float** hdev_walls_force_pp;
        unsigned int** hdev_contacts;


        //// Navier Stokes functions

        // Memory allocation
        void initNSmem();
        void freeNSmem();

        // Returns the number of fluid cells
        unsigned int NScells();         // Pressure and other centered nodes
        unsigned int NScellsVelocity(); // Inter-cell nodes (velocity)

        // Returns the mean particle radius
        Float meanRadius();

        // Get linear (1D) index from 3D coordinate
        unsigned int idx(const int x, const int y, const int z); // pres. nodes
        unsigned int vidx(const int x, const int y, const int z); // vel. nodes

        // Initialize Navier Stokes values and arrays
        void initNS();

        // Clean up Navier Stokes arrays
        void endNS();
        void endNSdev();

        // Check for stability in the FTCS solution
        void checkNSstability();

        // Returns the average value of the normalized residual norm in host mem
        double avgNormResNS();

        // Returns the maximum value of the normalized residual norm in host mem
        double maxNormResNS();

        // Allocate and free memory for NS arrays on device
        void initNSmemDev();
        void freeNSmemDev();

        // Transfer array values between GPU and CPU
        void transferNStoGlobalDeviceMemory(int statusmsg);
        void transferNSfromGlobalDeviceMemory(int statusmsg);
        void transferNSnormFromGlobalDeviceMemory();
        void transferNSepsilonFromGlobalDeviceMemory();
        void transferNSepsilonNewFromGlobalDeviceMemory();

        // Darcy values, host
        Darcy darcy;

        // Darcy values, device
        Float*  dev_darcy_p_old;     // Previous cell hydraulic pressure
        Float*  dev_darcy_dp_expl;   // Explicit integrated pressure change
        Float*  dev_darcy_p;         // Cell hydraulic pressure
        Float*  dev_darcy_p_new;     // Updated cell hydraulic pressure
        Float3* dev_darcy_v;         // Cell fluid velocity
        Float*  dev_darcy_phi;       // Cell porosity
        Float*  dev_darcy_dphi;      // Cell porosity change
        Float*  dev_darcy_div_v_p;   // Cell particle velocity divergence
        //Float*  dev_darcy_v_p_x;     // Cell particle velocity
        //Float*  dev_darcy_v_p_y;     // Cell particle velocity
        //Float*  dev_darcy_v_p_z;     // Cell particle velocity
        Float*  dev_darcy_norm;      // Normalized residual of epsilon values
        Float4* dev_darcy_f_p;       // Pressure gradient force on particles
        Float*  dev_darcy_k;         // Cell hydraulic permeability
        Float3* dev_darcy_grad_k;    // Spatial gradient of permeability
        Float3* dev_darcy_grad_p;    // Spatial gradient of fluid pressure
        Float3* dev_darcy_vp_avg;    // Average particle velocity in cell

        // Darcy functions
        void initDarcyMem();
        Float largestDarcyPermeability();
        Float smallestDarcyPorosity();
        Float3 largestDarcyVelocities();
        void initDarcyMemDev();
        unsigned int darcyCells();
        unsigned int darcyCellsVelocity();
        void transferDarcyToGlobalDeviceMemory(int statusmsg);
        void transferDarcyFromGlobalDeviceMemory(int statusmsg);
        void transferDarcyNormFromGlobalDeviceMemory();
        void transferDarcyPressuresFromGlobalDeviceMemory();
        void freeDarcyMem();
        void freeDarcyMemDev();
        unsigned int d_idx(const int x, const int y, const int z);
        unsigned int d_vidx(const int x, const int y, const int z);
        void checkDarcyStability();
        void printDarcyArray(FILE* stream, Float* arr);
        void printDarcyArray(FILE* stream, Float* arr, std::string desc);
        void printDarcyArray(FILE* stream, Float3* arr);
        void printDarcyArray(FILE* stream, Float3* arr, std::string desc);
        double avgNormResDarcy();
        double maxNormResDarcy();
        void initDarcy();
        void writeDarcyArray(Float* arr, const char* filename);
        void writeDarcyArray(Float3* arr, const char* filename);
        void endDarcy();
        void endDarcyDev();


    public:
        // Values and functions accessible from the outside

        // Constructor, some parameters with default values
        DEM(std::string inputbin,
                const int verbosity = 1,
                const int checkVals = 1,
                const int dry = 0,
                const int initCuda = 1,
                const int transferConstMem = 1,
                const int fluidFlow = 0,
                const int exclusive = 0);

        // Destructor
        ~DEM(void);

        // Read binary input file
        void readbin(const char *target);
        void readbinl(const char *target);
        // Write binary output file
        void writebin(const char *target);
        void writebinl(const char *target);

        // Check numeric values of selected parameters
        void diagnostics();
        void checkValues();

        // Report key parameter values to stdout
        void reportValues();

        // Iterate through time, using temporal limits
        // described in "time" struct.
        void startTime();

        // Render particles using raytracing
        void render(
                const int method = 1,
                const float maxval = 1.0e3f,
                const float lower_cutoff = 0.0f,
                const float focalLength = 1.0f,
                const unsigned int img_width = 800,
                const unsigned int img_height = 800);

        // Write image data to PPM file
        void writePPM(const char *target);

        // Calculate porosity with depth and save as text file
        void porosity(const int z_slices = 10);

        // find and return the min. position of any particle in each dimension
        Float3 minPos();

        // find and return the max. position of any particle in each dimension
        Float3 maxPos();

        // Find particle-particle intersections, saves the indexes
        // and the overlap sizes
        void findOverlaps(
                std::vector< std::vector<unsigned int> > &ij,
                std::vector< Float > &delta_n_ij);

        // Calculate force chains and save as Gnuplot script
        void forcechains(
                const std::string format = "interactive",
                const int threedim = 1,
                const double lower_cutoff = 0.0,
                const double upper_cutoff = 1.0e9);

        // Print all particle-particle contacts to stdout
        void printContacts();


        ///// Porous flow functions

        // Print fluid arrays to file stream
        void printNSarray(FILE* stream, Float* arr);
        void printNSarray(FILE* stream, Float* arr, std::string desc);
        void printNSarray(FILE* stream, Float3* arr);
        void printNSarray(FILE* stream, Float3* arr, std::string desc);

        // Write fluid arrays to file
        void writeNSarray(Float* array, const char* filename);
        void writeNSarray(Float3* array, const char* filename);
};

#endif
// vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
