/**
 * @file      nbody.h
 *
 * @author    Peter Lukac \n
 *            Faculty of Information Technology \n
 *            Brno University of Technology \n
 *            xlukac11@stud.fit.vutbr.cz
 *
 * @brief     PCG Assignment 2
 *            N-Body simulation in ACC
 *
 * @version   2021
 *
 * @date      11 November  2020, 11:22 (created) \n
 * @date      11 November  2020, 11:37 (revised) \n
 *
 */

#ifndef __NBODY_H__
#define __NBODY_H__

#include <cstdlib>
#include <cstdio>
#include  <cmath>
#include <cstring>
#include "h5Helper.h"

/// Gravity constant
constexpr float G = 6.67384e-11f;

/// Collision distance threshold
constexpr float COLLISION_DISTANCE = 0.01f;

/**
 * @struct float4
 * Structure that mimics CUDA float4
 */
struct float4
{
  float x;
  float y;
  float z;
  float w;
};

/// Define sqrtf from CUDA libm library
#pragma acc routine(sqrtf) seq

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                       Declare following structs / classes                                          //
//                                  If necessary, add your own classes / routines                                     //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Structure with particle data
 */
struct Particles
{
    float* pos_x;
    float* pos_y;
    float* pos_z;
    float* vel_x;
    float* vel_y;
    float* vel_z;
    float* weight;

    int N;
    int swap_count = 0;

    Particles (int N){
      this->N = N;

      pos_x = (float*)malloc(sizeof(float) * N);
      pos_y = (float*)malloc(sizeof(float) * N);
      pos_z = (float*)malloc(sizeof(float) * N);

      vel_x = (float*)malloc(sizeof(float) * N);
      vel_y = (float*)malloc(sizeof(float) * N);
      vel_z = (float*)malloc(sizeof(float) * N);

      weight = (float*)malloc(sizeof(float) * N);

      #pragma acc enter data copyin(this)
      #pragma acc enter data create(pos_x[0:N], pos_y[0:N], pos_z[0:N], vel_x[0:N], vel_y[0:N], vel_z[0:N], weight[0:N])
    }
    

    /**
     * Copies data to the other struct on CPU
     * @param [in]  p     - Destination Particles
     */
    void copyTo(Particles& p){
      memcpy(p.pos_x, pos_x, sizeof(float) * N);
      memcpy(p.pos_y, pos_y, sizeof(float) * N);
      memcpy(p.pos_z, pos_z, sizeof(float) * N);

      memcpy(p.vel_x, vel_x, sizeof(float) * N);
      memcpy(p.vel_y, vel_y, sizeof(float) * N);
      memcpy(p.vel_z, vel_z, sizeof(float) * N);

      memcpy(p.weight, weight, sizeof(float) * N);
    }

    void updateSelf(){
      #pragma acc update self(pos_x[0:N], pos_y[0:N], pos_z[0:N], vel_x[0:N], vel_y[0:N], vel_z[0:N], weight[0:N])
    }

    void updateDevice(){
      #pragma acc update device(pos_x[0:N], pos_y[0:N], pos_z[0:N], vel_x[0:N], vel_y[0:N], vel_z[0:N], weight[0:N])
    }

    ~Particles (){
      #pragma acc exit data delete(pos_x[0:N], pos_y[0:N], pos_z[0:N], vel_x[0:N], vel_y[0:N], vel_z[0:N], weight[0:N])
      #pragma acc exit data delete(this)

      free(pos_x);
      free(pos_y);
      free(pos_z);

      free(vel_x);
      free(vel_y);
      free(vel_z);

      free(weight);
    }

};// end of Particles
//----------------------------------------------------------------------------------------------------------------------

/**
 * @struct Velocities
 * Velocities of the particles
 */
struct Velocities
{
    float x;
    float y;
    float z;
    
};// end of Velocities
//----------------------------------------------------------------------------------------------------------------------

/**
 * Compute gravitation velocity
 * @param [in]  p_in     - Input Particles
 * @param [out] p_out    - Output Particles
 * @param [in ] N        - Number of particles
 * @param [in]  dt       - Time step size
 */
void calculate_velocity(Particles& p_in, Particles& p_out, const int N, const float dt);


/**
 * Compute gravitation velocity
 * @param [in]  p        - Particles
 * @param [out] tmp_vel  - Temporal velocity
 * @param [in ] N        - Number of particles
 * @param [in]  dt       - Time step size
 */
void calculate_gravitation_velocity(const Particles& p,
                                    Velocities&      tmp_vel,
                                    const int        N,
                                    const float      dt);

/**
 * Calculate collision velocity
 * @param [in]  p        - Particles
 * @param [out] tmp_vel  - Temporal velocity
 * @param [in ] N        - Number of particles
 * @param [in]  dt       - Time step size
 */
void calculate_collision_velocity(const Particles& p,
                                  Velocities&      tmp_vel,
                                  const int        N,
                                  const float      dt);

/**
 * Update particle position
 * @param [in]  p        - Particles
 * @param [out] tmp_vel  - Temporal velocity
 * @param [in ] N        - Number of particles
 * @param [in]  dt       - Time step size
 */
void update_particle(Particles& p,
                     const Velocities&      tmp_vel,
                     const int        N,
                     const float      dt);



/**
 * Compute center of gravity - implement in steps 3 and 4.
 * @param [in] p - Particles
 * @param [in] N - Number of particles
 * @return Center of Mass [x, y, z] and total weight[w]
 */
float4 centerOfMassGPU(const Particles& p,
                       const int        N);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Compute center of mass on CPU
 * @param memDesc
 * @return centre of gravity
 */
float4 centerOfMassCPU(MemDesc& memDesc);

#endif /* __NBODY_H__ */
