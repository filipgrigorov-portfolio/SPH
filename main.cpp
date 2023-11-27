#include <iostream>

#include "sph.hpp"

#include <GL/glut.h>

namespace {
static constexpr int WINDOW_WIDTH = 800;
static constexpr int WINDOW_HEIGHT = 600;
static constexpr double FRUSTRUM_WIDTH = 1.5 * WINDOW_WIDTH;
static constexpr double FRUSTRUM_HEIGHT = 1.5 * WINDOW_HEIGHT;

static constexpr float particleDiameter = 30.f;

static constexpr float RED = 0.05490196078f;
static constexpr float GREEN = 0.5294117647f;
static constexpr float BLUE = 1.0f;
static constexpr float ALPHA = 1;

static constexpr uint32_t CUBE_SIZE = 20;
static std::vector<fluid::Particle> fluidParticles;
static fluid::SPH sph(FRUSTRUM_HEIGHT, FRUSTRUM_WIDTH, CUBE_SIZE,
                      particleDiameter, fluidParticles);
} // namespace

void initOpenGL() {
  // Note: select background color
  glClearColor(0.9f, 0.9f, 0.9f, 1);

  // Note: Remove aliasing
  glEnable(GL_POINT_SMOOTH);

  glPointSize(particleDiameter / 2.f);

  // Note: initialize viewing values
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  // Note: multiply the current matrix with an orthographic matrix
  glOrtho(0.0, FRUSTRUM_WIDTH, 0.0, FRUSTRUM_HEIGHT, 0.0, 1.0);
}

// Note: Incovenient but it uses global static variables
void renderFunction() {
  // Clear all pixels
  glClear(GL_COLOR_BUFFER_BIT);

  // Note: Color of the fluid particles, R,G,B,Alpha
  glColor4f(RED, GREEN, BLUE, ALPHA);

  glBegin(GL_POINTS);
  {
    // Note: Go through the particles
    std::cout << "Rendering " << fluidParticles.size() << " particles"
              << std::endl;
    for (const auto &fluidParticle : fluidParticles) {
      std::cout << fluidParticle.position(0) << "," << fluidParticle.position(1)
                << std::endl;
      glVertex2d(fluidParticle.position(0), fluidParticle.position(1));
    }
  }
  glEnd();

  glutSwapBuffers();
}

int main(int argc, char **argv) {
  glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
  glutInit(&argc, argv);
  glutCreateWindow("SPH 3D");

  // Note: Intialize the openGL parameters
  initOpenGL();

  // Note: Display and refresh of window
  glutDisplayFunc(renderFunction);
  // glutIdleFunc(udateStep);
  // glutKeyboardFunc(keyboardEvent);

  glutMainLoop();

  return EXIT_SUCCESS;
}
