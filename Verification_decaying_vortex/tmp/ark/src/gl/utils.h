#if defined(__APPLE__)
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#else
#  include <GL/gl.h>
#  include <GL/glu.h>
#endif
#include <stdio.h>
#include <stdbool.h>

void gl_write_image (FILE * fp, const GLubyte * buffer,
		     unsigned width, unsigned height, unsigned samples);
void init_gl();

void matrix_multiply (float * m, float * n);
void vector_multiply (float * v, float * m);

typedef struct {
  float m[16], p[16];
  float n[6][3];
  float d[6];
  unsigned width;
#if 0
  GList * symmetries;
  FttVector * s;
#endif
} Frustum;

void gl_get_frustum (Frustum * f);
int sphere_in_frustum (double x, double y, double z, double r, Frustum * f);
float sphere_diameter (double x, double y, double z, double r, Frustum * f);
void gl_check_error();

enum FeedbackFormat {
  FEEDBACK_GNU, FEEDBACK_OBJ, FEEDBACK_KML
};

float * gl_feedback_begin (unsigned buffersize);
bool    gl_feedback_end   (float * f,
			   FILE * fp,
			   enum FeedbackFormat format);

int polygonize (const double val[8], double isolevel, double triangles[5][3][3]);

#include "parser.h"
