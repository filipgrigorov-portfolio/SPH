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
    mSurfaceTension = fluidProps["surface_tension"];

    fluidParticles.reserve(cubeSize * cubeSize);
    const auto startRowPos = 0.1 * camHeight;
    const auto startColPos = 0.4 * camWidth;
    const float& increment = particleSize;
    const auto endRowPos = startRowPos + (cubeSize * increment);
    const auto endColPos = startColPos + (cubeSize * increment);
    for (float row = startRowPos; row < endRowPos; row += increment) {
      for (float col = startColPos; col < endColPos; col += increment) {
        const float rowUpdated = row / startRowPos;
        const float colUpdated = col / startColPos;
        const uint32_t idx = uint32_t(rowUpdated * cubeSize + colUpdated);
        const Eigen::Vector2f randomPerturbation = {((float)rand() / (RAND_MAX)) + 1, ((float)rand() / (RAND_MAX)) + 1};
        fluidParticles.push_back({Eigen::Vector2f(col, row) + randomPerturbation,
                                  Eigen::Vector2f(0.0, 0.0),
                                  Eigen::Vector2f(0.0, 0.0), 0.f, 0.f});
      }
    }

    std::cout << "Generated " << fluidParticles.size() << " particles" << std::endl;
}

void SPH::setSemiCircularBoundary(const std::string& cliArg) {
    if (!cliArg.empty()) {
        if (cliArg == "on") {
            mIsSemiCircularBoundaryEnabled = true;
        }
    }
}

void SPH::setSemiCircularBoundary(const std::string& cmdArg1, const std::string& cmdArg2) {
    if (!cmdArg1.empty() && !cmdArg2.empty()) {
        if (cmdArg1 == "on") {
            mIsSemiCircularBoundaryEnabled = true;
        }
        
        if (cmdArg2 == "with_penetration") {
            mIsSemiCircularBoundaryRunWithPenetration = true;
        }
    }
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
            surfaceTensionForce += -mSurfaceTension * viscosityLambda(l, h) * unitDirection;
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

        // Note: Not that stable, it requires further work but I wanted to experiment
        if (mIsSemiCircularBoundaryEnabled) {
            // BCs: center location, o = [0.5 * boundaryWidth, 0.01 * boundaryHeight] -> circular BC
            static constexpr float X_OFFSET = 0.5f;
            static constexpr float Y_OFFSET = 0.01f;
            static const Eigen::Vector2f bcOrigin{X_OFFSET * boundaryWidth, Y_OFFSET * boundaryHeight};
            static float bcRadius = 5 * mParticleSize;
            Eigen::Vector2f normal = bcOrigin - particle.position;
            Eigen::Vector2f unitNormal = -normal.normalized() * dampingScale;
            if (normal.squaredNorm() <= bcRadius * bcRadius) {
                particle.velocityField(0) *= unitNormal(0);
                particle.velocityField(1) *= unitNormal(1);

                if (mIsSemiCircularBoundaryRunWithPenetration) {
                    particle.position += dt * particle.velocityField;
                } else {
                    const float x = std::fabs(particle.position(0) - bcOrigin(0));
                    const float y = std::fabs(particle.position(1) - bcOrigin(1));
                    const float theta = std::atan2(y, x);

                    // Note: Before the 90 degree mark
                    if (theta < M_PI / 2) {
                        particle.position(0) += bcRadius * std::cos(theta);
                        particle.position(1) = bcRadius + mParticleSize * std::sin(theta);
                    }

                    // Note: After the 90 degree mark
                    if (theta > M_PI / 2) {
                        particle.position(0) = bcRadius * std::cos(theta);
                        particle.position(1) = bcRadius + mParticleSize * std::sin(theta);
                    }

                    // Note: At the 90 degree mark
                    if (theta == M_PI / 2) {
                        particle.position(1) = bcRadius;
                    }

                    // Note: At theta = 0 or pi
                    if (theta == 0 || theta == M_PI) {
                        particle.position(0) = bcRadius;
                    }
                }
            }
        }
    }
}
} // namespace fluid
