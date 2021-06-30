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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                       Declare following structs / classes                                          //
//                                  If necessary, add your own classes / routines                                     //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Compute gravitation velocity
void calculate_gravitation_velocity(const Particles& restrict p,
                                    Velocities&      restrict tmp_vel,
                                    const int        N,
                                    const float      dt)
{
  #pragma acc parallel loop gang present(p, tmp_vel)
  for (int i = 0; i < N; i++) {

    float pos_x = (&p)[i].pos_x;
    float pos_y = (&p)[i].pos_y;
    float pos_z = (&p)[i].pos_z;
    float weight = (&p)[i].weight;

    float tmp_vel_x = 0.0;
    float tmp_vel_y = 0.0;
    float tmp_vel_z = 0.0;

    #pragma acc loop seq
    for (int j = 0; j < N; j++) {

      float dx = pos_x - (&p)[j].pos_x;
      float dy = pos_y - (&p)[j].pos_y;
      float dz = pos_z - (&p)[j].pos_z;

      float r = sqrt(dx*dx + dy*dy + dz*dz);


      float F = -G * weight * (&p)[j].weight / (r * r + FLT_MIN);

      float nx = F * dx/ (r + FLT_MIN);
      float ny = F * dy/ (r + FLT_MIN);
      float nz = F * dz/ (r + FLT_MIN);

      float vx = nx * dt / weight;
      float vy = ny * dt / weight;
      float vz = nz * dt / weight;

      tmp_vel_x += (r > COLLISION_DISTANCE) ? vx : 0.0f;
      tmp_vel_y += (r > COLLISION_DISTANCE) ? vy : 0.0f;
      tmp_vel_z += (r > COLLISION_DISTANCE) ? vz : 0.0f;
    }

    (&tmp_vel)[i].x = tmp_vel_x;
    (&tmp_vel)[i].y = tmp_vel_y;
    (&tmp_vel)[i].z = tmp_vel_z;
  }
  

}// end of calculate_gravitation_velocity
//----------------------------------------------------------------------------------------------------------------------

void calculate_collision_velocity(const Particles& restrict p,
                                  Velocities&      restrict tmp_vel,
                                  const int        N,
                                  const float      dt)
{
  
  #pragma acc parallel loop gang present(p, tmp_vel)
  for (int i = 0; i < N; i++) {

    float pos_x = (&p)[i].pos_x;
    float pos_y = (&p)[i].pos_y;
    float pos_z = (&p)[i].pos_z;

    float vel_x = (&p)[i].vel_x;
    float vel_y = (&p)[i].vel_y;
    float vel_z = (&p)[i].vel_z;

    float weight = (&p)[i].weight;

    float tmp_vel_x = 0.0;
    float tmp_vel_y = 0.0;
    float tmp_vel_z = 0.0;

    #pragma acc loop seq
    for (int j = 0; j < N; j++) {

      float dx = pos_x - (&p)[j].pos_x;
      float dy = pos_y - (&p)[j].pos_y;
      float dz = pos_z - (&p)[j].pos_z;

      float r = sqrtf(dx*dx + dy*dy + dz*dz);

      float vx = ((weight* vel_x - (&p)[j].weight *vel_x + 2* (&p)[j].weight* (&p)[j].vel_x) /
              (weight + (&p)[j].weight)) - vel_x ;
      float vy = ((weight* vel_y - (&p)[j].weight *vel_y + 2* (&p)[j].weight* (&p)[j].vel_y) /
              (weight + (&p)[j].weight)) - vel_y ;
      float vz = ((weight* vel_z - (&p)[j].weight *vel_z + 2* (&p)[j].weight* (&p)[j].vel_z) /
              (weight + (&p)[j].weight)) - vel_z ;

      if (r > 0.0f && r < COLLISION_DISTANCE) {
          tmp_vel_x += vx;
          tmp_vel_y += vy;
          tmp_vel_z += vz;
      }
    }
    (&tmp_vel)[i].x += tmp_vel_x;
    (&tmp_vel)[i].y += tmp_vel_y;
    (&tmp_vel)[i].z += tmp_vel_z;
  }

}// end of calculate_collision_velocity
//----------------------------------------------------------------------------------------------------------------------

/// Update particle position
void update_particle(Particles& restrict p,
                     const Velocities& restrict tmp_vel,
                     const int        N,
                     const float      dt)
{

  #pragma acc parallel loop gang present(p, tmp_vel)
  for (int i = 0; i < N; i++) {
    (&p)[i].vel_x += (&tmp_vel)[i].x;
    (&p)[i].vel_y += (&tmp_vel)[i].y;
    (&p)[i].vel_z += (&tmp_vel)[i].z;

    (&p)[i].pos_x += (&p)[i].vel_x * dt;
    (&p)[i].pos_y += (&p)[i].vel_y * dt;
    (&p)[i].pos_z += (&p)[i].vel_z * dt;
  }

}// end of update_particle
//----------------------------------------------------------------------------------------------------------------------



/// Compute center of gravity
float4 centerOfMassGPU(const Particles& p,
                       const int        N)
{

  return {0.0f, 0.0f, 0.0f, 0.0f};
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
