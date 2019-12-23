#define OSX

#ifdef _WIN32
#   include "windows.h"
#else
#   include <sys/time.h>
#endif

#ifdef OSX
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#endif

#include <chrono>
#include <thread>

#include <cstdlib>
#include <time.h>
#include <math.h>
#include <vector>
#include <iostream>
#include "algebra3.h"
#define M_PI 3.14159265

#define min(x, y) (x < y ? x : y)
#define max(x, y) (x > y ? x : y)
#define abs(x) (x < 0 ? -x : x)

typedef std::chrono::steady_clock Clock;
typedef std::chrono::duration<float> Duration;
typedef vec2 Point;

class Timer
{
public:
  Clock::time_point start_clock, tmp_clock, last_clock;
  Duration duration;

  Timer()
  {
    reset_timer();
  }

  void reset_timer ()
  {
    start_clock = last_clock = Clock::now();
  }

  double get_gap ()
  {
    tmp_clock = Clock::now();
    duration = tmp_clock - last_clock;
    last_clock = tmp_clock;

    return duration.count();
  }

  double get_interval ()
  {
    duration = Clock::now() - start_clock;
    return duration.count();
  }
};

Timer timer;

class ColorProperty
{
public:
  double Kr, Kg, Kb;
  double Sr, Sg, Sb;
  double rho;     // density
  double omega;   // stain power
  double gamma;   // granularity

  ColorProperty (double Kr, double Kg, double Kb, double Sr, double Sg, double Sb, double rho, double omega, double gamma):
    Kr(Kr), Kg(Kg), Kb(Kb), Sr(Sr), Sg(Sg), Sb(Sb), rho(rho), omega(omega), gamma(gamma) {}
};

class Color
{
public:
  double g;       // shallow water pigment
  double d;       // deposited pigment

  Color (): g(0), d(0) {}
};

const double C_MIN = 0.5;
const double C_MAX = 0.9;
const int COLOR = 12;
class Paper
{
public:
  double h;       // height
  double c;       // capacity
  double p;       // pressure
  double u;       // velocity x
  double v;       // velocity y
  int M;          // wetness mask
  double s;       // water saturation
  Color color[COLOR]; // pigment
  double R;
  double G;
  double B;

  Paper ()
  {
    h = (rand() % 9999) / 10000.0 + 0.0001;
    c = C_MIN + h * (C_MAX - C_MIN);
    M = 0;
  }
};

/**********************************************
 *                Parameters
 **********************************************/

const int BEZIER = 25;
const int PAPER_SIZE = 250; //500;
const int POINT_SIZE = 8; //4;
const ColorProperty COLOR_PROPERTY[] = {
  ColorProperty(0.22, 1.47, 0.57, 0.05, 0.003, 0.03, 0.02, 5.5, 0.81),
  ColorProperty(0.46, 1.07, 1.50, 1.28, 0.38, 0.21, 0.05, 7.0, 0.40),
  ColorProperty(0.10, 0.36, 3.45, 0.97, 0.65, 0.007, 0.05, 3.4, 0.81),
  ColorProperty(1.62, 0.61, 1.64, 0.01, 0.012, 0.003, 0.09, 1.0, 0.41),
  ColorProperty(1.52, 0.32, 0.25, 0.06, 0.26, 0.40, 0.01, 1.0, 0.31),
  ColorProperty(0.74, 1.54, 2.10, 0.09, 0.09, 0.004, 0.09, 9.3, 0.90),
  ColorProperty(0.14, 1.08, 1.68, 0.77, 0.015, 0.018, 0.02, 1.0, 0.63),
  ColorProperty(0.13, 0.81, 3.45, 0.005, 0.009, 0.007, 0.01, 1.0, 0.14),
  ColorProperty(0.06, 0.21, 1.78, 0.50, 0.88, 0.009, 0.06, 1.0, 0.08),
  ColorProperty(1.55, 0.47, 0.63, 0.01, 0.05, 0.035, 0.02, 1.0, 0.12),
  ColorProperty(0.86, 0.86, 0.06, 0.005, 0.005, 0.09, 0.01, 3.1, 0.91),
  ColorProperty(0.08, 0.11, 0.07, 1.25, 0.42, 1.43, 0.06, 1.0, 0.08)
};

int window_width = 600;
int window_height = 600;
int display_size, diff;

int brush_size = 5;       // brush's radius
double brush_ppt = 5;     // brush's pressure per time
double brush_g = 1;       // brush's pigment concentation
int brush_color = 0;      // brush's color id
Paper paper[PAPER_SIZE][PAPER_SIZE];

/**********************************************
 *                   Utility
 **********************************************/

// double w1, w2;
// int ii1, ii2, jj1, jj2;
// double get_u (double i, double j)
// {
//   ii1 = i;
//   w1 = 1 - (i - ii1);
//   ii2 = ii1 + 1;

//   jj1 = j;
//   w2 = 1 - (j - jj1);
//   jj2 = jj1 + 1;

//   ii1 = min(max(ii1, 0.0), PAPER_SIZE - 1.0);
//   ii2 = min(max(ii2, 0.0), PAPER_SIZE - 1.0);
//   jj1 = min(max(jj1, 0.0), PAPER_SIZE - 1.0);
//   jj2 = min(max(jj2, 0.0), PAPER_SIZE - 1.0);

//   return w1 * (w2 * paper[ii1][jj1].u + (1 - w2) * paper[ii1][jj2].u) +
//    (1 - w1) * (w2 * paper[ii2][jj1].u + (1 - w2) * paper[ii2][jj2].u);
// }

// double get_v (double i, double j)
// {
//   ii1 = i;
//   w1 = 1 - (i - ii1);
//   ii2 = ii1 + 1;

//   jj1 = j;
//   w2 = 1 - (j - jj1);
//   jj2 = jj1 + 1;

//   ii1 = min(max(ii1, 0), PAPER_SIZE - 1);
//   ii2 = min(max(ii2, 0), PAPER_SIZE - 1);
//   jj1 = min(max(jj1, 0), PAPER_SIZE - 1);
//   jj2 = min(max(jj2, 0), PAPER_SIZE - 1);

//   return w1 * (w2 * paper[ii1][jj1].v + (1 - w2) * paper[ii1][jj2].v) +
//    (1 - w1) * (w2 * paper[ii2][jj1].v + (1 - w2) * paper[ii2][jj2].v);
// }

std::vector <Point> stroke;
int stroke_count;
bool pressed;
void paint (Point p)
{
  p.x = (p.x * (PAPER_SIZE - 1));
  p.y = (p.y * (PAPER_SIZE - 1));
  for (int i = -brush_size + 1; i < brush_size; ++i)
  {
    if (p.x + i < 0 || p.x + i >= PAPER_SIZE)
      continue;

    for (int j = -brush_size + 1; j < brush_size; ++j)
    {
      if (i * i + j * j >= brush_size * brush_size)
        continue;
      if (p.y + j < 0 || p.y + j >= PAPER_SIZE)
        continue;

      stroke.push_back(Point(p.x + i, p.y + j));
      paper[(int)(p.x + i)][(int)(p.y + j)].M = 1;
    }
  }

  ++stroke_count;
}

std::vector <Point> dots;
Point dot1, dot2, dot3, new_dot;
int draw_state;
double factor;
void add_point (Point p)
{
  if (draw_state == 0)
  {
    dots.push_back(p);
    paint(p);
    draw_state = 1;
  }
  else if (draw_state == 1)
  {
    dot1 = dots.back();
    dot2 = p;

    // Linear interpolation
    for (int i = 1; i <= BEZIER; ++i)
    {
      new_dot = (i / (double) BEZIER) * dot1 + 
        ((BEZIER - i) / (double) BEZIER) * dot2;
      dots.push_back(new_dot);
      paint(new_dot);
    }

    draw_state = 2;
  }
  else
  {
    dot1 = dots[dots.size() - BEZIER];
    dot2 = dots[dots.size() - 1];
    dot3 = p;

    dot2 = dot2 - (dot3 - dot1) / 4.0;
    dot2 = 2 * dot2 - dot1;

    // Bezier curve
    for (int i = 1; i <= BEZIER; ++i)
    {
      factor = 0.5 + 0.5 * i / (double) BEZIER;
      new_dot = (1 - factor) * ((1 - factor) * dot1 + factor * dot2) +
         factor * ((1 - factor) * dot2 + factor * dot3);
      dots.push_back(new_dot);
      paint(new_dot);
    }
  }
}

void setPixel(int x, int y, GLdouble r, GLdouble g, GLdouble b) {
  glColor3f(r, g, b);
  glVertex2f(x, y);
}

void setPixel(double x, double y, GLdouble r, GLdouble g, GLdouble b)
{
  setPixel((int)(x * display_size),
    diff + (int)((1 - y) * display_size), r, g, b);
}

// Display
void display() {
  glClear(GL_COLOR_BUFFER_BIT);

  diff = window_height - display_size;

  // Fill outside the square
  glColor3f(0.5, 0.5, 0.5);
  glRecti(display_size, 0, window_width, window_height);  // right
  glRecti(0, 0, window_width, diff);                      // bottom

  // Draw dots
  glPointSize(POINT_SIZE);
  glBegin(GL_POINTS);
  for (int i = 0; i < PAPER_SIZE; ++i)
  {
    for (int j = 0; j < PAPER_SIZE; ++j)
    {
      // factor = paper[i][j].h;
      // if (paper[i][j].M == 1)
      // {
      //   factor *= min(paper[i][j].p, 1.0);
      //   setPixel((double) i / (PAPER_SIZE - 1), 
      //     (double) j / (PAPER_SIZE - 1), 
      //     factor * 0.35, factor * 0.6, factor * 0.8);
      // }
      // else
      //   setPixel((double) i / (PAPER_SIZE - 1), 
      //     (double) j / (PAPER_SIZE - 1), 
      //     factor, factor, factor);

      factor = paper[i][j].color[brush_color].d;
      setPixel((double) i / (PAPER_SIZE - 1), 
          (double) j / (PAPER_SIZE - 1), 
          factor, factor, factor);

      // if (factor != 0)
      //   cout << i << " " << j << ": " << factor << endl;

      setPixel((double) i / (PAPER_SIZE - 1), 
          (double) j / (PAPER_SIZE - 1), 
          paper[i][j].R, paper[i][j].G, paper[i][j].B);

      // if (paper[i][j].R + paper[i][j].G + paper[i][j].B > 0)
      //   cout << i << " " << j << ": " << paper[i][j].R << " " << paper[i][j].G << " " << paper[i][j].B << endl;
    }
  }

  glEnd();
  glFlush();
  glutSwapBuffers();
}

// Press/release function
void processMouse(int button, int state, int x, int y)
{
  if ((button == GLUT_LEFT_BUTTON) && (state == GLUT_DOWN))
  {
    pressed = true;
    add_point(Point((double) x / display_size, 
      (double) y / display_size));
  }
  else
  {
    pressed = false;
    dots.clear();
    draw_state = 0;
  }
}

// Drag function
void drag (int x, int y)
{
  add_point(Point((double) x / display_size, 
    (double) y / display_size));
}

void reshape(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
   // Compute aspect ratio of the new window
  if (height == 0) height = 1;                // To prevent divide by 0

  window_width = width;
  window_height = height;
  display_size = min(window_width, window_height);

  // Set the viewport to cover the new window
  glViewport(0, 0, width, height);

  // Set the aspect ratio of the clipping area to match the viewport
  glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
  glLoadIdentity();             // Reset the projection matrix
  gluOrtho2D(0, width, 0, height);
  display();
}

void Initialize() {
    glEnable(GL_MULTISAMPLE);
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glMatrixMode(GL_PROJECTION);

    window_width = glutGet(GLUT_WINDOW_WIDTH);
    window_height = glutGet(GLUT_WINDOW_HEIGHT);
    display_size = min(window_width, window_height);

    glViewport(0,0, window_width, window_height);
    glLoadIdentity();
    //Display area
    gluOrtho2D(0.0, window_width, 0.0, window_height);
}

void check()
{
  for (int i = 0; i < PAPER_SIZE; ++i)
  {
    for (int j = 0; j < PAPER_SIZE; ++j)
    {
      if (paper[i][j].u != 0 || paper[i][j].v != 0)
      {
        cout << "IN: " << paper[i][j].p << " | " << paper[i][j].u << " " << paper[i][j].v << endl;
        return;
      }
    }
  }
}

/**********************************************
 *           Watercolor algorithm
 **********************************************/

Point last_dot;
int x, y;
void ApplyStroke (double dt)
{
  if (stroke.size() == 0)
    return;

  last_dot = stroke[stroke.size() - 1];

  for (int i = 0; i < stroke.size(); ++i)
  {
    x = stroke[i].x;
    y = stroke[i].y;
    paper[x][y].p += brush_ppt * dt / stroke_count;
    paper[x][y].color[brush_color].g = brush_g;
    paper[x][y].M = 1;
  }

  stroke.clear();
  stroke_count = 0;

  if (pressed)
    paint(dots[dots.size() - 1]);
}

const double MU = 0.1;
const double KAPPA = 0.01;

double uu[PAPER_SIZE + 1][PAPER_SIZE + 1];
double vv[PAPER_SIZE + 1][PAPER_SIZE + 1];
double u1, u2, u3, u4, u5;
double v1, v2, v3, v4, v5;
double A, B;
void UpdateVelocity (double dt)
{
  // calculate gradient
  for (int i = 0; i < PAPER_SIZE; ++i)
  {
    for (int j = 0; j < PAPER_SIZE; ++j)
    {
      if (i == 0)
        vv[i][j] = (paper[i + 1][j].h - paper[i][j].h);
      else if (i == PAPER_SIZE - 1)
        vv[i][j] = (paper[i][j].h - paper[i - 1][j].h);
      else
        vv[i][j] = (paper[i + 1][j].h - paper[i - 1][j].h) / 2;

      if (j == 0)
        uu[i][j] = (paper[i][j + 1].h - paper[i][j].h);
      else if (j == PAPER_SIZE - 1)
        uu[i][j] = (paper[i][j].h - paper[i][j - 1].h);
      else
        uu[i][j] = (paper[i][j + 1].h - paper[i][j - 1].h) / 2;
    }
  }

  // apply gradient
  for (int i = 0; i < PAPER_SIZE; ++i)
  {
    for (int j = 0; j < PAPER_SIZE; ++j)
    {
      paper[i][j].u -= uu[i][j];
      paper[i][j].v -= vv[i][j];
    }
  }

  // calculate new velocity
  for (int i = 0; i < PAPER_SIZE; ++i)
  {
    for (int j = 0; j < PAPER_SIZE; ++j)
    {
      // u1 = get_u(i, j);
      u1 = paper[i][j].u;
      // u2 = get_u(i, j + 1);
      u2 = paper[i][min(j + 1, PAPER_SIZE - 1)].u;
      // u3 = get_u(i - 0.5, j + 0.5);
      u3 = (paper[max(i - 1, 0)][j].u + paper[max(i - 1, 0)][min(j + 1, PAPER_SIZE - 1)].u + paper[i][j].u + paper[i][min(j + 1, PAPER_SIZE - 1)].u) / 4;
      // v3 = get_v(i - 0.5, j + 0.5);
      v3 = u3 = (paper[max(i - 1, 0)][j].v + paper[max(i - 1, 0)][min(j + 1, PAPER_SIZE - 1)].v + paper[i][j].v + paper[i][min(j + 1, PAPER_SIZE - 1)].v) / 4;
      // u4 = get_u(i + 0.5, j + 0.5);
      u4 = (paper[i][j].u + paper[i][min(j + 1, PAPER_SIZE - 1)].u + paper[min(i + 1, PAPER_SIZE - 1)][j].u + paper[min(i + 1, PAPER_SIZE - 1)][min(j + 1, PAPER_SIZE - 1)].u) / 4;
      // v4 = get_v(i + 0.5, j + 0.5);
      v4 = (paper[i][j].v + paper[i][min(j + 1, PAPER_SIZE - 1)].v + paper[min(i + 1, PAPER_SIZE - 1)][j].v + paper[min(i + 1, PAPER_SIZE - 1)][min(j + 1, PAPER_SIZE - 1)].v) / 4;
      A = u1*u1 - u2*u2 + u3*v3 - u4*v4;

      // u1 = get_u(i, j + 1.5);
      u1 = (paper[i][min(j + 1, PAPER_SIZE - 1)].u + paper[i][min(j + 2, PAPER_SIZE - 1)].u) / 2;
      // u2 = get_u(i, j - 0.5);
      u2 = (paper[i][max(j - 1, 0)].u + paper[i][j].u) / 2;
      // u3 = get_u(i + 1, j + 0.5);
      u3 = (paper[min(i + 1, PAPER_SIZE - 1)][j].u + paper[min(i + 1, PAPER_SIZE - 1)][min(j + 1, PAPER_SIZE - 1)].u) / 2;
      // u4 = get_u(i - 1, j + 0.5);
      u4 = (paper[max(i - 1, 0)][j].u + paper[max(i - 1, 0)][min(j + 1, PAPER_SIZE - 1)].u) / 2;
      // u5 = get_u(i, j + 0.5);
      u5 = (paper[i][j].u + paper[i][min(j + 1, PAPER_SIZE - 1)].u) / 2;
      B = u1 + u2 + u3 + u4 - 4*u5;

      uu[i][j] = u5 + dt * (A - MU*B + paper[i][j].p
        - paper[i][min(j + 1, PAPER_SIZE - 1)].p - KAPPA*u5);

      // v1 = get_v(i, j);
      v1 = paper[i][j].v;
      // v2 = get_v(i + 1, j);
      v2 = paper[min(i + 1, PAPER_SIZE)][j].v;
      // u3 = get_u(i + 0.5, j - 0.5);
      u3 = (paper[i][max(j - 1, 0)].u + paper[i][j].u + paper[min(i + 1, PAPER_SIZE - 1)][max(j - 1, 0)].u + paper[min(i + 1, PAPER_SIZE - 1)][j].u) / 4;
      // v3 = get_v(i + 0.5, j - 0.5);
      v3 = (paper[i][max(j - 1, 0)].v + paper[i][j].v + paper[min(i + 1, PAPER_SIZE - 1)][max(j - 1, 0)].v + paper[min(i + 1, PAPER_SIZE - 1)][j].v) / 4;
      // u4 = get_u(i + 0.5, j + 0.5);
      u4 = (paper[i][j].u + paper[i][min(j + 1, PAPER_SIZE - 1)].u + paper[min(i + 1, PAPER_SIZE - 1)][j].u + paper[min(i + 1, PAPER_SIZE - 1)][min(j + 1, PAPER_SIZE - 1)].u) / 4;
      // v4 = get_v(i + 0.5, j + 0.5);
      v4 = (paper[i][j].v + paper[i][min(j + 1, PAPER_SIZE - 1)].v + paper[min(i + 1, PAPER_SIZE - 1)][j].v + paper[min(i + 1, PAPER_SIZE - 1)][min(j + 1, PAPER_SIZE - 1)].v) / 4;
      A = v1*v1 - v2*v2 + u3*v3 - u4*v4;

      // v1 = get_v(i + 0.5, j + 1);
      v1 = (paper[i][min(j + 1, PAPER_SIZE - 1)].v + paper[min(i + 1, PAPER_SIZE - 1)][min(j + 1, PAPER_SIZE - 1)].v) / 2;
      // v2 = get_v(i + 0.5, j - 1);
      v2 = (paper[i][max(j - 1, 0)].v + paper[min(i + 1, PAPER_SIZE - 1)][max(j - 1, 0)].v) / 2;
      // v3 = get_v(i + 1.5, j);
      v3 = (paper[min(i + 1, PAPER_SIZE - 1)][j].v + paper[min(i + 2, PAPER_SIZE - 1)][j].v) / 2;
      // v4 = get_v(i - 0.5, j);
      v4 = (paper[max(i - 1, 0)][j].v + paper[i][j].v) / 2;
      // v5 = get_v(i + 0.5, j);
      v5 = (paper[i][j].v + paper[min(i + 1, PAPER_SIZE - 1)][j].v) / 2;
      B = v1 + v2 + v3 + v4 - 4*v5;

      vv[i][j] = v5 + dt * (A - MU*B + paper[i][j].p
        - paper[min(i + 1, PAPER_SIZE - 1)][j].p - KAPPA*v5);

      cout << uu[i][j] << " " << vv[i][j] << endl;
    }
  }

  // update velocity
  for (int i = PAPER_SIZE - 1; i >= 0; --i)
  {
    for (int j = PAPER_SIZE - 1; j >= 0; --j)
    {
      if (i == PAPER_SIZE - 1)
        paper[i][j].v = vv[i][j];
      else
        paper[i][j].v = (2 * vv[i][j] - paper[i + 1][j].v);

      if (j == PAPER_SIZE - 1)
        paper[i][j].u = uu[i][j];
      else
        paper[i][j].u = (2 * uu[i][j] - paper[i][j + 1].u);
    }
  }

  // EnforceBoundaryCondition
  cout << "ENFORCE!!!!" << endl;
  for (int i = 0; i < PAPER_SIZE; ++i)
  {
    for (int j = 0; j < PAPER_SIZE; ++j)
    {
      if (paper[i][j].M == 0)
      {
        paper[i][j].u = 0;
        paper[i][j].v = 0;
      }
      else
      {
        cout << i << " " << j << ": " << paper[i][j].u << " " << paper[i][j].v << endl;
      }
      
    }
  }
}

const int N = 5; //50;
const double TAU = 0.01;
double delta, delta_max;
void RelaxDivergence ()
{
  for (int n = 0; n < N; ++n)
  {
    for (int i = 0; i < PAPER_SIZE; ++i)
    {
      for (int j = 0; j < PAPER_SIZE; ++j)
      {
        u1 = (paper[i][j].u + paper[i][min(j + 1, PAPER_SIZE - 1)].u) / 2;
        u2 = (paper[i][max(j - 1, 0)].u + paper[i][j].u) / 2;
        v3 = (paper[i][j].v + paper[min(i + 1, PAPER_SIZE - 1)][j].v) / 2;
        v4 = (paper[max(i - 1, 0)][j].v + paper[i][j].v) / 2;
        delta = KAPPA * (u1 - u2 + v3 - v4);

        uu[i][j + 1] += delta;
        uu[i][j] -= delta;
        vv[i + 1][j] += delta;
        vv[i][j] -= delta;

        delta_max = max(delta_max, abs(delta));
      }
    }

    check();

    for (int i = PAPER_SIZE - 1; i >= 0; --i)
    {
      for (int j = PAPER_SIZE - 1; j >= 0; --j)
      {
        if (i == PAPER_SIZE - 1)
          paper[i][j].v = vv[i][j];
        else
          paper[i][j].v = (2 * vv[i][j] - paper[i + 1][j].v);

        if (j == PAPER_SIZE - 1)
          paper[i][j].u = uu[i][j];
        else
          paper[i][j].u = (2 * uu[i][j] - paper[i][j + 1].u);
      }
    }

    if (delta_max <= TAU)
      break;
  }
}

const int K = 9;    // Make it odd, not even
const double ETA = 0.01;
double filter[K][K], M;
void EdgeDarkening ()
{
  for (int i = 0; i < PAPER_SIZE; ++i)
  {
    for (int j = 0; j < PAPER_SIZE; ++j)
    {
      // Apply Gaussian blur
      M = 0;
      for (int k = 0; k < K; ++k)
      {
        for (int l = 0; l < K; ++l)
        {
          if (i + k - K / 2 < 0 || i + k - K / 2 >= PAPER_SIZE || j + l - K / 2 < 0 || j + l - K / 2 >= PAPER_SIZE)
            continue;
          M += filter[k][l] * paper[i + k - K / 2][j + l - K / 2].M;
        }
      }

      // paper[i][j].p -= ETA * (1 - M) * paper[i][j].M * paper[i][j].p;
      paper[i][j].p -= ETA * (1 - M) * paper[i][j].M;
      paper[i][j].p = max(paper[i][j].p, 0);
    }
  }
}

double gg[PAPER_SIZE + 1][PAPER_SIZE + 1];
double g;
void MovingPigment (double dt)
{
  for (int i = 0; i <= PAPER_SIZE; ++i)
  {
    for (int j = 0; j <= PAPER_SIZE; ++j)
    {
      uu[i][j] = (paper[i][max(j - 1, 0)].u + paper[i][min(j, PAPER_SIZE - 1)].u) / 2;
      vv[i][j] = (paper[max(i - 1, 0)][j].v + paper[min(i, PAPER_SIZE - 1)][j].v) / 2;

      // if (uu[i][j] != 0)
      //   cout << i << " " << j << ": " << uu[i][j] << endl;
    }
  }

  for (int c = 0; c < COLOR; ++c)
  {
    for (int i = 0; i < PAPER_SIZE; ++i)
    {
      for (int j = 0; j < PAPER_SIZE; ++j)
        gg[i][j] = paper[i][j].color[c].g;
    }

    for (int i = 0; i < PAPER_SIZE; ++i)
    {
      for (int j = 0; j < PAPER_SIZE; ++j)
      {
        g = gg[i][j];
        v1 = max(0, g * uu[i][j + 1]);
        paper[i][min(j + 1, PAPER_SIZE - 1)].color[c].g += v1;
        v2 = max(0, -g * uu[i][j]);
        paper[i][max(j - 1, 0)].color[c].g += v2;
        v3 = max(0, g * vv[i + 1][j]);
        paper[min(i + 1, PAPER_SIZE - 1)][j].color[c].g += v3;
        v4 = max(0, -g * vv[i][j]);
        paper[max(i - 1, 0)][j].color[c].g += v4;
        // Suspicious typo in paper
        paper[i][j].color[c].g = paper[i][j].color[c].g - v1 - v2 - v3 - v4;
      }
    }
  }
}

double delta_up, delta_down;
void TransterPigment ()
{
  for (int c = 0; c < COLOR; ++c)
  {
    for (int i = 0; i < PAPER_SIZE; ++i)
    {
      for (int j = 0; j < PAPER_SIZE; ++j)
      {
        if (paper[i][j].M == 0)
          continue;
        
        delta_down = paper[i][j].color[c].g * (1 - paper[i][j].h * COLOR_PROPERTY[c].gamma) * COLOR_PROPERTY[c].rho;
        delta_up = paper[i][j].color[c].d * (1 + (paper[i][j].h - 1) * COLOR_PROPERTY[c].gamma) * COLOR_PROPERTY[c].rho / COLOR_PROPERTY[c].omega;

        if (paper[i][j].color[c].d + delta_down > 1)
          delta_down = max(0, 1 - paper[i][j].color[c].d);
        if (paper[i][j].color[c].g + delta_up > 1)
          delta_up = max(0, 1 - paper[i][j].color[c].g);

        // cout << i << " " << j << ": " << paper[i][j].color[c].d << " " << paper[i][j].color[c].g << endl;

        paper[i][j].color[c].d += delta_down - delta_up;
        paper[i][j].color[c].g += delta_up - delta_down;

        // cout << delta_down << " " << delta_up << endl;
        // cout << i << " " << j << ": " << paper[i][j].color[c].d << " " << paper[i][j].color[c].g << endl;
      }
    }
  }
}

const double ALPHA = 1;
const double EPSILON = 0.2;
const double DELTA = -1;
const double SIGMA = 0.7;
double ss[PAPER_SIZE + 1][PAPER_SIZE + 1];
double ds;
void SimulateCapillaryFlow ()
{
  for (int i = 0; i < PAPER_SIZE; ++i)
  {
    for (int j = 0; j < PAPER_SIZE; ++j)
    {
      if (paper[i][j].M > 0)
        paper[i][j].s += max(0, min(ALPHA, paper[i][j].c - paper[i][j].s));
      ss[i][j] = paper[i][j].s;
    }
  }

  for (int i = 0; i < PAPER_SIZE; ++i)
  {
    for (int j = 0; j < PAPER_SIZE; ++j)
    {
      if (paper[i][j].s < EPSILON)
        continue;

      if (i > 0 && paper[i][j].s > paper[i - 1][j].s && paper[i - 1][j].s > DELTA)
      {
        ds = max(0, min(ss[i][j] - ss[i - 1][j], paper[i - 1][j].c - ss[i - 1][j]) / 4);
        paper[i][j].s -= ds;
        paper[i - 1][j].s += ds;
      }

      if (i < PAPER_SIZE - 1 && paper[i][j].s > paper[i + 1][j].s && paper[i + 1][j].s > DELTA)
      {
        ds = max(0, min(ss[i][j] - ss[i + 1][j], paper[i + 1][j].c - ss[i + 1][j]) / 4);
        paper[i][j].s -= ds;
        paper[i + 1][j].s += ds;
      }

      if (j > 0 && paper[i][j].s > paper[i][j - 1].s && paper[i][j - 1].s > DELTA)
      {
        ds = max(0, min(ss[i][j] - ss[i][j - 1], paper[i][j - 1].c - ss[i][j - 1]) / 4);
        paper[i][j].s -= ds;
        paper[i][j - 1].s += ds;
      }

      if (j < PAPER_SIZE - 1 && paper[i][j].s > paper[i][j + 1].s && paper[i][j + 1].s > DELTA)
      {
        ds = max(0, min(ss[i][j] - ss[i][j + 1], paper[i][j + 1].c - ss[i][j + 1]) / 4);
        paper[i][j].s -= ds;
        paper[i][j + 1].s += ds;
      }
    }
  }

  for (int i = 0; i < PAPER_SIZE; ++i)
  {
    for (int j = 0; j < PAPER_SIZE; ++j)
    {
      if (paper[i][j].s > SIGMA)
        paper[i][j].M = 1;
    }
  }
}

// double ap_sinh (double x)
// {
//   return x + (x * x * x) / 6;
// }

// double ap_cosh (double x)
// {
//   return 1 + (x * x) / 2;
// }

#define ap_sinh(x) (x + (x * x * x) / 6)
#define ap_cosh(x) (1 + (x * x) / 2)

const double a = 1;
const double b = 1;
double cc, R, T, Rr, Rg, Rb, Tr, Tg, Tb, param_tmp, param, inv;
void CalculateColor ()
{
  for (int i = 0; i < PAPER_SIZE; ++i)
  {
    for (int j = 0; j < PAPER_SIZE; ++j)
    {
      paper[i][j].R = paper[i][j].G = paper[i][j].B = 0;
      for (int c = 0; c < COLOR; ++c)
      {
        param_tmp = b * paper[i][j].color[c].d;
        param = COLOR_PROPERTY[c].Sr * param_tmp;
        cc = a * sinh(param) + b * cosh(param);
        R = sinh(param / cc);
        inv = 1 / (1 - Rr * R);
        T = b / cc;
        Rr = Rr + Tr * Tr * R * inv;
        Tr = Tr * T * inv;

        param = COLOR_PROPERTY[c].Sg * param_tmp;
        cc = a * sinh(param) + b * cosh(param);
        R = sinh(param / cc);
        inv = 1 / (1 - Rg * R);
        T = b / cc;
        Rg = Rg + Tg * Tg * R * inv;
        Tg = Tg * T * inv;

        param = COLOR_PROPERTY[c].Sb * param_tmp;
        cc = a * sinh(param) + b * cosh(param);
        R = sinh(param / cc);
        inv = 1 / (1 - Rb * R);
        T = b / cc;
        Rb = Rb + Tb * Tb * R * inv;
        Tb = Tb * T * inv;
      }

      paper[i][j].R = Rr;
      paper[i][j].G = Rg;
      paper[i][j].B = Rb;
    }
  }
}

void preprocess ()
{
  // EdgeDarkening: Gaussian blur filter
  for (int i = 0; i < K; ++i)
  {
    filter[i][0] = filter[0][i] = 1;
    for (int j = i - 1; j >= 1; --j)
      filter[j][0] = filter[0][j] = filter[0][j - 1] + filter[0][j];
  }

  double sm = 0;
  for (int i = 0; i < K; ++i)
  {
    for (int j = 0; j < K; ++j)
    {
      filter[i][j] = filter[i][0] * filter[0][j];
      sm += filter[i][j];
    }
  }

  for (int i = 0; i < K; ++i)
  {
    for (int j = 0; j < K; ++j)
      filter[i][j] /= sm;
  }
}

double dt;
void MainLoop (int unused) {
  // Watercolor algorithm here
  dt = timer.get_gap();

  ApplyStroke(dt);
  UpdateVelocity(dt);
  RelaxDivergence();
  EdgeDarkening();
  MovingPigment(dt);
  TransterPigment();
  SimulateCapillaryFlow();
  CalculateColor();

  display();

  glutTimerFunc(1000 / 60, MainLoop, 0);
}

int main(int iArgc, char** cppArgv) {
  srand(time(NULL));
  preprocess();

  glutInit(&iArgc, cppArgv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
  glutInitWindowSize(window_width, window_height);
  glutInitWindowPosition(100, 100);
  glutCreateWindow(cppArgv[0]);

  glutReshapeFunc(reshape);
  glutMotionFunc(drag);
  glutMouseFunc(processMouse);

  MainLoop(0);

  Initialize();
  glutDisplayFunc(display);
  glutMainLoop();

  return 0;
}