/**
 * @file      nbody.cpp
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

#include <math.h>
#include <cfloat>
#include "nbody.h"
#include <chrono>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                       Declare following structs / classes                                          //
//                                  If necessary, add your own classes / routines                                     //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void calculate_velocity(Particles& restrict p_in, Particles& restrict p_out, const int N, const float dt){
  #pragma acc parallel loop gang present(p_in, p_out)
  for (int i = 0; i < N; i++) {

    float tmp_vel_x = 0;
    float tmp_vel_y = 0;
    float tmp_vel_z = 0;

    float pos_x = p_in.pos_x[i];
    float pos_y = p_in.pos_y[i];
    float pos_z = p_in.pos_z[i];

    float vel_x = p_in.vel_x[i];
    float vel_y = p_in.vel_y[i];
    float vel_z = p_in.vel_z[i];

    float weight = p_in.weight[i];

    #pragma acc loop seq
    for (int j = 0; j < N; j++) {
      // Same For both
      float dx = pos_x - p_in.pos_x[j];
      float dy = pos_y - p_in.pos_y[j];
      float dz = pos_z - p_in.pos_z[j];

      // non coliding velocities
      float r = sqrt(dx*dx + dy*dy + dz*dz);

      float r3 = r * r * r + FLT_MIN;

      float G_dt_r3 = -G * dt / r3;
      float Fg_dt_m2_r = G_dt_r3 * p_in.weight[j];

      float vx = Fg_dt_m2_r * dx;
      float vy = Fg_dt_m2_r * dy;
      float vz = Fg_dt_m2_r * dz;

      // non coliding velocity
      tmp_vel_x += (r > COLLISION_DISTANCE) ? vx : 0.0f;
      tmp_vel_y += (r > COLLISION_DISTANCE) ? vy : 0.0f;
      tmp_vel_z += (r > COLLISION_DISTANCE) ? vz : 0.0f;


      // coliding velocities
      vx = ((weight * vel_x - p_in.weight[j] * vel_x + 2 * p_in.weight[j] * p_in.vel_x[j]) /
        (weight + p_in.weight[j])) - vel_x;
      vy = ((weight * vel_y - p_in.weight[j] * vel_y + 2 * p_in.weight[j] * p_in.vel_y[j]) /
        (weight + p_in.weight[j])) - vel_y;
      vz = ((weight * vel_z - p_in.weight[j] * vel_z + 2 * p_in.weight[j] * p_in.vel_z[j]) /
        (weight + p_in.weight[j])) - vel_z;


      if (r > 0.0f && r < COLLISION_DISTANCE) {
        tmp_vel_x += vx;
        tmp_vel_y += vy;
        tmp_vel_z += vz;
      }
    }
    // update particle
    p_out.vel_x[i] = p_in.vel_x[i] + tmp_vel_x;
    p_out.vel_y[i] = p_in.vel_y[i] + tmp_vel_y;
    p_out.vel_z[i] = p_in.vel_z[i] + tmp_vel_z;

    p_out.pos_x[i] = p_in.pos_x[i] + p_out.vel_x[i] * dt;
    p_out.pos_y[i] = p_in.pos_y[i] + p_out.vel_y[i] * dt;
    p_out.pos_z[i] = p_in.pos_z[i] + p_out.vel_z[i] * dt;

  }
}


/// Compute gravitation velocity
void calculate_gravitation_velocity(const Particles& restrict p,
                                    Velocities&      restrict tmp_vel,
                                    const int        N,
                                    const float      dt)
{

}// end of calculate_gravitation_velocity
//----------------------------------------------------------------------------------------------------------------------

void calculate_collision_velocity(const Particles& restrict p,
                                  Velocities&      restrict tmp_vel,
                                  const int        N,
                                  const float      dt)
{
 
}// end of calculate_collision_velocity
//----------------------------------------------------------------------------------------------------------------------

/// Update particle position
void update_particle(Particles& p,
                     const Velocities&      tmp_vel,
                     const int        N,
                     const float      dt)
{

}// end of update_particle
//----------------------------------------------------------------------------------------------------------------------


/// Compute center of gravity
float4 centerOfMassGPU(const Particles& restrict p,
                       const int        N)
{ 
  // N_2 is N rounded up to the next power of 2
  int N_2 = N;
  N_2--;
  for (int i = 1; i < sizeof(int); i <<= 1)
    N_2 |= N_2 >> i;
  N_2++;

  // copy of the particle data, because they will be overwritten in reduction
  float* x = (float*)malloc(sizeof(float) * N);
  float* y = (float*)malloc(sizeof(float) * N);
  float* z = (float*)malloc(sizeof(float) * N);
  float* w = (float*)malloc(sizeof(float) * N);

  #pragma acc enter data create(x[0:N], y[0:N], z[0:N], w[0:N])

  #pragma acc parallel loop gang present(p)
  for (int i = 0; i < N; i++){
    x[i] = p.pos_x[i];
    y[i] = p.pos_y[i];
    z[i] = p.pos_z[i];
    w[i] = p.weight[i];
  }

  // pre-scan inspired implementation of COM reduction
  // main loop is called log(N) times
  // main loop has no pragma, therefore creates implicit barrier for the threads
  for(int i = 1; i < N_2; i <<= 1) {

    // parallel section calculates one step in tree
    #pragma acc parallel loop gang present(x, y, z, w)
    for (int j = 0; j < N;  j+=2){
      if(i + j < N) {
        const float dx = x[j + i] - x[j];
        const float dy = y[j + i] - y[j];
        const float dz = z[j + i] - z[j];

        float dw = ((w[j + i] + w[j]) > 0.0f)
                        ? ( w[j + i] / (w[j + i] + w[j])) : 0.0f;

        x[j] += dx * dw;
        y[j] += dy * dw;
        z[j] += dz * dw;
        w[j] += w[j + i];
      }
    }
  }

  #pragma acc exit data copyout(x[0:1], y[0:1], z[0:1], w[0:1])

  return {x[0], y[0], z[0], w[0]};

}// end of centerOfMassGPU
//----------------------------------------------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Compute center of mass on CPU
float4 centerOfMassCPU(MemDesc& memDesc)
{
  float4 com = {0 ,0, 0, 0};

  for(int i = 0; i < memDesc.getDataSize(); i++)
  {
    // Calculate the vector on the line connecting points and most recent position of center-of-mass
    const float dx = memDesc.getPosX(i) - com.x;
    const float dy = memDesc.getPosY(i) - com.y;
    const float dz = memDesc.getPosZ(i) - com.z;

    // Calculate weight ratio only if at least one particle isn't massless
    const float dw = ((memDesc.getWeight(i) + com.w) > 0.0f)
                          ? ( memDesc.getWeight(i) / (memDesc.getWeight(i) + com.w)) : 0.0f;

    // Update position and weight of the center-of-mass according to the weight ration and vector
    com.x += dx * dw;
    com.y += dy * dw;
    com.z += dz * dw;
    com.w += memDesc.getWeight(i);
  }
  return com;
}// end of centerOfMassCPU
//----------------------------------------------------------------------------------------------------------------------
