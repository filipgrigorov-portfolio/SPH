#include "sph.hpp"

#include <cmath>
#include <iostream>

namespace {
constexpr float GRAVITY_ACC = 9.81f;
float clamp(float currentValue, float minValue, float maxValue) {
    return std::max(minValue, std::min(currentValue, maxValue));
}

// Note: According to "SPH Based Shallow Water Simulation" by Solenthaler
/// @brief Smoothing kernels
auto poly6Lambda = [](float lSquared, float h) -> float {
    const float hSquared = h * h;
    return (lSquared <= hSquared) ? 4.f / (M_PI * std::pow(h, 8)) * std::pow(hSquared - lSquared, 3.f) : 0.f;
};
auto spikyLambda = [](float l, float h) -> float {
    return (l <= h) ? -10.f / (M_PI * std::pow(h, 5)) * std::pow(h - l, 3) : 0.f;
};
auto viscosityLambda = [](float l, float h) -> float {
    return (l <= h) ? 40.f / (M_PI * std::pow(h, 5)) * (h - l) : 0.f;
};
} // namespace

namespace fluid {
/// @brief This class represents the "average" fluid attributes for each
/// particle based on its neighbours in a given radius R
SPH::SPH(float camWidth, float camHeight, uint32_t cubeSize, float particleSize, 
  std::unordered_map<std::string, float>& fluidProps, std::vector<Particle> &fluidParticles)
: mParticleSize(particleSize) {
  // Note: Initialize the global fluid properties
  mRestDensity = fluidProps["rest_density"];
  mGasConst = fluidProps["gas_const"];
  mViscosity = fluidProps["viscosity"];

  fluidParticles.reserve(cubeSize * cubeSize);
  const auto startRowPos = 0.2 * camHeight;
  const auto startColPos = 0.1 * camWidth;
  const float& increment = particleSize;
  const auto endRowPos = startRowPos + (cubeSize * increment);
  const auto endColPos = startColPos + (cubeSize * increment);
  for (float row = startRowPos; row < endRowPos; row += increment) {
    for (float col = startColPos; col < endColPos; col += increment) {
      const float rowUpdated = row / startRowPos;
      const float colUpdated = col / startColPos;
      const uint32_t idx = uint32_t(rowUpdated * cubeSize + colUpdated);
      float jitter = static_cast<float>(random()) / static_cast<float>(RAND_MAX);
      fluidParticles.push_back({Eigen::Vector2f(col + jitter * 0.5, row),
                                Eigen::Vector2f(0.0, 0.0),
                                Eigen::Vector2f(0.0, 0.0), 0.f, 0.f});
      std::cout << fluidParticles.back().position(0) << " ";//debug
    }
  }

  std::cout << "Generated " << fluidParticles.size() << " particles" << std::endl;
}

void SPH::accumulateForces(std::vector<Particle> &fluidParticles) {
    const float& h = mParticleSize;
    for (auto& particle : fluidParticles) {
        Eigen::Vector2f pressureForce{0.f, 0.f};
        Eigen::Vector2f viscousForce{0.f, 0.f};
        Eigen::Vector2f gravityForce{0.f, -GRAVITY_ACC};
        Eigen::Vector2f surfaceTensionForce{0.f, 0.f};

        for (const auto& neighbour : fluidParticles) {
            if (neighbour.position == particle.position) {
                continue;
            }

            const auto distVector = neighbour.position - particle.position;
            const float l = distVector.norm();
            const Eigen::Vector2f unitDirection = distVector.normalized();

            // pressure forces
            const float averagePressure = (particle.pressure + neighbour.pressure) / (2.f * neighbour.density);
            pressureForce -= neighbour.mass * averagePressure * spikyLambda(l, h) * unitDirection;
            
            // viscous forces
            viscousForce += neighbour.mass * ((neighbour.velocityField - particle.velocityField) / neighbour.density) * viscosityLambda(l, h);

            // surface tension
        }

        viscousForce *= mViscosity;
        particle.forceField = pressureForce + viscousForce + gravityForce * particle.mass / particle.density + surfaceTensionForce;
    }
}

void SPH::accumulateParticlesDensity(std::vector<Particle> &fluidParticles) {
    const float& h = mParticleSize;
    for (auto &particle : fluidParticles) {
        particle.density = 0.f; // init at every step        
        for (const auto &neighbour : fluidParticles) {
            const float lSquared = (neighbour.position - particle.position).squaredNorm();
            particle.density += neighbour.mass * poly6Lambda(lSquared, h);
        }

        // Update pressure (Tait equation, suggested by Desbrun)
        particle.pressure = mGasConst * (particle.density - mRestDensity);
    }
}

// Note: We need to address boundary conditions, here
// Reverse the velocity fields if they hit the walls of the container
/*
  * a_t+1 = F_i / density_i = dv_i/dt
  * v_t+1 = v_t + dt * F_tot / density
  * p_t+1 = p_t + dt * v_t+1
*/
void SPH::timeIntegrate(std::vector<Particle> &fluidParticles, float dt, float dampingScale, float boundaryWidth, float boundaryHeight) {
  for (auto &particle : fluidParticles) {
    particle.velocityField += dt * particle.forceField / particle.density;
    particle.position += dt * particle.velocityField;

    const float minValue = mParticleSize;
    const float maxXValue = boundaryWidth - minValue;
    const float maxYValue = boundaryHeight - minValue;
    // Left or right
    if (particle.position(0) <= minValue || particle.position(0) >= maxXValue) {
      // Note: Damping is applied as energy is lost to the boundary
      particle.velocityField(0) *= -dampingScale;
      particle.position(0) = clamp(particle.position(0), minValue, maxXValue);
    }

    // Top or bottom
    if (particle.position(1) <= minValue || particle.position(1) >= maxYValue) {
      // Note: Damping is applied as energy is lost to the boundary
      particle.velocityField(1) *= -dampingScale;
      particle.position(1) = clamp(particle.position(1), minValue, maxXValue);
    }
  }
}
} // namespace fluid
