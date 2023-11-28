#include "sph.hpp"

#include <cmath>
#include <iostream>

namespace {
static constexpr float PI = 3.14159265359f;
struct Poly6 {
  float operator()(const Eigen::Vector2f &ri, const Eigen::Vector2f &rj,
                   float h) {
    const auto distance = (ri - rj).squaredNorm();
    static const float kPoly6Const = 315.f / (64.f * PI * std::pow(h, 6));
    return kPoly6Const * std::pow(distance, 3.f);
  }
};

template <typename KernelType>
float W(const Eigen::Vector2f &ri, const Eigen::Vector2f &rj, float h) {
  return KernelType().operator()(ri, rj, h);
}

std::vector<fluid::Particle>
findNeighbouringParticles(const fluid::Particle &particle,
                          const std::vector<fluid::Particle> &fluidParticles,
                          float radius) {
  std::vector<fluid::Particle> neighbours;
  neighbours.reserve(20);
  // Naive implementation
  for (const auto &potentialNeighbour : fluidParticles) {
    if ((particle.position - potentialNeighbour.position).norm() <= radius) {
      neighbours.emplace_back(particle);
    }
  }
  return neighbours;
}
} // namespace

namespace fluid {
/// @brief This class represents the "average" fluid attributes for each
/// particle based on its neighbours in a given radius R
SPH::SPH(double camWidth, double camHeight, uint32_t cubeSize,
         float particleDiameter, std::vector<Particle> &fluidParticles) {
  fluidParticles.reserve(cubeSize * cubeSize);
  const auto startRowPos = 0.2 * camHeight;
  const auto startColPos = 0.1 * camWidth;
  const float increment = particleDiameter;
  const auto endRowPos = startRowPos + (cubeSize * increment);
  const auto endColPos = startColPos + (cubeSize * increment);
  for (float row = startRowPos; row < endRowPos; row += increment) {
    for (float col = startColPos; col < endColPos; col += increment) {
      const float rowUpdated = row / startRowPos;
      const float colUpdated = col / startColPos;
      const uint32_t idx = uint32_t(rowUpdated * cubeSize + colUpdated);
      // std::cout << col << "," << row << std::endl;
      fluidParticles.push_back({Eigen::Vector2f(col, row),
                                Eigen::Vector2f(0.0, 0.0),
                                Eigen::Vector2f(0.0, 0.0), 0.f, 0.f});
    }
  }
  std::cout << "Generated " << fluidParticles.size() << " particles"
            << std::endl;
}

void SPH::updateParticlesDensity(std::vector<Particle> &fluidParticles,
                                 float density, float mass, float h,
                                 float gasConstant, float radius) {
  for (auto &particle : fluidParticles) {
    particle.density = 0.f; // init at every step
    std::vector<Particle> neighbours =
        findNeighbouringParticles(particle, fluidParticles, radius);
    std::cout << "Found " << neighbours.size() << " neighbours" << std::endl;
    for (auto neighbour : neighbours) {
      particle.density +=
          mass * W<Poly6>(particle.position, neighbour.position, h);
    }
    // Update pressure (Tait equation)
    particle.pressure = gasConstant * (particle.density - density);
  }
}

// Note: We need to address boundary conditions, here
// Reverse the velocity fields if they hit the walls of the container
// TODO:
void SPH::timeIntegrate(std::vector<Particle> &fluidParticles, float dt,
                        float mass, float containerWidth,
                        float containerHeight) {
  for (auto &particle : fluidParticles) {
    Eigen::Vector2f newVelocity =
        particle.velocityField + dt * particle.forceField / mass;
    Eigen::Vector2f newPosition =
        particle.position + dt * particle.velocityField;

    // Left or right
    if (newPosition(0) <= 0 || newPosition(0) >= containerWidth - 1) {
      newVelocity(0) = -newVelocity(0);
    }

    // Top or bottom
    if (newPosition(1) <= 0 || newPosition(1) >= containerHeight - 1) {
      newVelocity(1) = -newVelocity(1);
    }
    particle.velocityField = newVelocity;
    particle.position = newPosition;
  }
}
} // namespace fluid
