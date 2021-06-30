/**
 * @file      main.cpp
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

#include <chrono>
#include <cstdio>
#include <cmath>

#include "nbody.h"
#include "h5Helper.h"

/**
 * Main routine of the project
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char **argv)
{
  // Parse command line parameters
  if (argc != 7)
  {
    printf("Usage: nbody <N> <dt> <steps> <write intesity> <input> <output>\n");
    exit(EXIT_FAILURE);
  }

  const int   N         = std::stoi(argv[1]);
  const float dt        = std::stof(argv[2]);
  const int   steps     = std::stoi(argv[3]);
  const int   writeFreq = (std::stoi(argv[4]) > 0) ? std::stoi(argv[4]) : 0;

  printf("N: %d\n", N);
  printf("dt: %f\n", dt);
  printf("steps: %d\n", steps);

  const size_t recordsNum = (writeFreq > 0) ? (steps + writeFreq - 1) / writeFreq : 0;


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //                                         Code to be implemented                                                   //
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // 1.  Memory allocation on CPU
  Particles* particles_in = new Particles(N);
  Particles* particles_out = new Particles(N);

  // 2. Create memory descriptor
  /*
   * Caution! Create only after CPU side allocation
   * parameters:
   *                                    Stride of two               Offset of the first
   *             Data pointer           consecutive elements        element in floats,
   *                                    in floats, not bytes        not bytes
  */
  MemDesc md(
              particles_in->pos_x,                1,                          0,            // Position in X
              particles_in->pos_y,                1,                          0,            // Position in Y
              particles_in->pos_z,                1,                          0,            // Position in Z
              particles_in->vel_x,                1,                          0,            // Velocity in X
              particles_in->vel_y,                1,                          0,            // Velocity in Y
              particles_in->vel_z,                1,                          0,            // Velocity in Z
              particles_in->weight,               1,                          0,            // Weight
              N,                                                                // Number of particles
              recordsNum);                                                      // Number of records in output file



  H5Helper h5Helper(argv[5], argv[6], md);

  // Read data
  try
  {
    h5Helper.init();
    h5Helper.readParticleData();
  } catch (const std::exception& e)
  {
    std::cerr<<e.what()<<std::endl;
    return EXIT_FAILURE;
  }

  // 3. Copy data to GPU
  // This copies particles to the target struct (on CPU)
  particles_in->copyTo(*particles_out);
  
  particles_in->updateDevice();
  particles_out->updateDevice();


  // Start the time
  auto startTime = std::chrono::high_resolution_clock::now();


  // 4. Run the loop - calculate new Particle positions.
  for (int s = 0; s < steps; s++)
  {

    calculate_velocity(*particles_in, *particles_out, N, dt);

    // swap pointers
    Particles* tmp = particles_in;
    particles_in = particles_out;
    particles_out = tmp;

    /// In step 4 - fill in the code to store Particle snapshots.
    if (writeFreq > 0 && (s % writeFreq == 0))
    {
      
      /*float4 comOnGPU = {0.0f, 0.0f, 0.0f, 0.f};
      
      h5Helper.writeParticleData(s / writeFreq);
      h5Helper.writeCom(comOnGPU.x, comOnGPU.y, comOnGPU.z, comOnGPU.w, s / writeFreq);*/
    }
  }// for s ...

  // 5. In steps 3 and 4 -  Compute center of gravity
  float4 comOnGPU = {0.0f, 0.0f, 0.0f, 0.f};



  // Stop watchclock
  const auto   endTime = std::chrono::high_resolution_clock::now();
  const double time    = (endTime - startTime) / std::chrono::milliseconds(1);
  printf("Time: %f s\n", time / 1000);


  // 5. Copy data from GPU back to CPU.
  particles_in->updateSelf();
  particles_out->updateSelf();

  // if the count of interations is odd, then particles_in is particles_out and we have copy data from the last iteration
  if ( steps % 2 == 1 )
    // in fact it is particles_out -> particles_in
    particles_in->copyTo(*particles_out);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /// Calculate center of gravity
  float4 comOnCPU = centerOfMassCPU(md);


    std::cout<<"Center of mass on CPU:"<<std::endl
      << comOnCPU.x <<", "
      << comOnCPU.y <<", "
      << comOnCPU.z <<", "
      << comOnCPU.w
      << std::endl;

    std::cout<<"Center of mass on GPU:"<<std::endl
      << comOnGPU.x <<", "
      << comOnGPU.y <<", "
      << comOnGPU.z <<", "
      << comOnGPU.w
      <<std::endl;

  // Store final positions of the particles into a file
  h5Helper.writeComFinal(comOnGPU.x, comOnGPU.y, comOnGPU.z, comOnGPU.w);
  h5Helper.writeParticleDataFinal();

  return EXIT_SUCCESS;
}// end of main
//----------------------------------------------------------------------------------------------------------------------

