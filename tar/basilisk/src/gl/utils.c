#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "utils.h"

/**

## Various helper functions

A helper function to write a PPM file from a RGB buffer. Downsampling
by averaging is performed if *samples* is larger than one. */

void gl_write_image (FILE * fp, const GLubyte * buffer,
		     unsigned width, unsigned height, unsigned samples)
{
  const GLubyte *ptr = buffer;

  if (samples < 1)
    samples = 1;
  if (samples > 4)
    samples = 4;

  width /= samples, height /= samples;
  fprintf (fp, "P6 %d %d 255\n", width, height);
  int x, y, j, k;
  for (y = height - 1; y >= 0; y--)
    for (x = 0; x < width; x++) {
      int r = 0, g = 0, b = 0;
      for (j = 0; j < samples; j++)
	for (k = 0; k < samples; k++) {
	  int i = (((y*samples + j)*width + x)*samples + k)*4;
	  r += ptr[i]; g += ptr[i+1]; b += ptr[i+2];
	}
      fputc (r/samples/samples, fp); /* write red */
      fputc (g/samples/samples, fp); /* write green */
      fputc (b/samples/samples, fp); /* write blue */
    }
}

/**
This is the basic OpenGL setup. */

void init_gl() {
  GLfloat light0_pos[4]   = { 0.0, 0.0, 50.0, 0.0 };
  GLfloat light0_color[4] = { 1., 1., 1., 1.0 }; /* white light */

  glDisable (GL_CULL_FACE);
  glEnable (GL_DEPTH_TEST);
  glEnable (GL_NORMALIZE);

  /* speedups */
  glEnable (GL_DITHER);
  glShadeModel (GL_SMOOTH);
  glHint (GL_PERSPECTIVE_CORRECTION_HINT, GL_FASTEST);
  glHint (GL_POLYGON_SMOOTH_HINT, GL_FASTEST);

  /* light */
  glLightModeli (GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  glLightfv (GL_LIGHT0, GL_POSITION, light0_pos);
  glLightfv (GL_LIGHT0, GL_DIFFUSE,  light0_color);
  glEnable (GL_LIGHT0);
  glEnable (GL_LIGHTING);

  glColorMaterial (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glEnable (GL_COLOR_MATERIAL);
}

void gl_draw_texture (GLuint id, int width, int height)
{
  glMatrixMode (GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho (0.0, width, 0.0, height, -1.0, 1.0);
  glMatrixMode (GL_MODELVIEW);
  glPushMatrix();

  glLoadIdentity();
  glDisable (GL_LIGHTING);

  glColor3f (1,1,1);
  glEnable (GL_TEXTURE_2D);
  glBindTexture (GL_TEXTURE_2D, id);

  glBegin(GL_QUADS);
  glTexCoord2f(0, 0); glVertex3f (0, 0, 0);
  glTexCoord2f(0, 1); glVertex3f (0, 100, 0);
  glTexCoord2f(1, 1); glVertex3f (100, 100, 0);
  glTexCoord2f(1, 0); glVertex3f (100, 0, 0);
  glEnd();

  glDisable(GL_TEXTURE_2D);
  glPopMatrix();

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  glMatrixMode(GL_MODELVIEW);
}

#define RC(r,c) m[(r)+(c)*4]
#define RCM(m,r,c) (m)[(r)+(c)*4]

void matrix_multiply (float * m, float * n)
{
  float o[16];
  int i;
  for (i = 0; i < 16; i++) o[i] = m[i];
  RC(0,0)=RCM(o,0,0)*RCM(n,0,0)+RCM(o,0,1)*RCM(n,1,0)+
          RCM(o,0,2)*RCM(n,2,0)+RCM(o,0,3)*RCM(n,3,0);
  RC(0,1)=RCM(o,0,0)*RCM(n,0,1)+RCM(o,0,1)*RCM(n,1,1)+
          RCM(o,0,2)*RCM(n,2,1)+RCM(o,0,3)*RCM(n,3,1);
  RC(0,2)=RCM(o,0,0)*RCM(n,0,2)+RCM(o,0,1)*RCM(n,1,2)+
          RCM(o,0,2)*RCM(n,2,2)+RCM(o,0,3)*RCM(n,3,2);
  RC(0,3)=RCM(o,0,0)*RCM(n,0,3)+RCM(o,0,1)*RCM(n,1,3)+
          RCM(o,0,2)*RCM(n,2,3)+RCM(o,0,3)*RCM(n,3,3);
  RC(1,0)=RCM(o,1,0)*RCM(n,0,0)+RCM(o,1,1)*RCM(n,1,0)+
          RCM(o,1,2)*RCM(n,2,0)+RCM(o,1,3)*RCM(n,3,0);
  RC(1,1)=RCM(o,1,0)*RCM(n,0,1)+RCM(o,1,1)*RCM(n,1,1)+
          RCM(o,1,2)*RCM(n,2,1)+RCM(o,1,3)*RCM(n,3,1);
  RC(1,2)=RCM(o,1,0)*RCM(n,0,2)+RCM(o,1,1)*RCM(n,1,2)+
          RCM(o,1,2)*RCM(n,2,2)+RCM(o,1,3)*RCM(n,3,2);
  RC(1,3)=RCM(o,1,0)*RCM(n,0,3)+RCM(o,1,1)*RCM(n,1,3)+
          RCM(o,1,2)*RCM(n,2,3)+RCM(o,1,3)*RCM(n,3,3);
  RC(2,0)=RCM(o,2,0)*RCM(n,0,0)+RCM(o,2,1)*RCM(n,1,0)+
          RCM(o,2,2)*RCM(n,2,0)+RCM(o,2,3)*RCM(n,3,0);
  RC(2,1)=RCM(o,2,0)*RCM(n,0,1)+RCM(o,2,1)*RCM(n,1,1)+
          RCM(o,2,2)*RCM(n,2,1)+RCM(o,2,3)*RCM(n,3,1);
  RC(2,2)=RCM(o,2,0)*RCM(n,0,2)+RCM(o,2,1)*RCM(n,1,2)+
          RCM(o,2,2)*RCM(n,2,2)+RCM(o,2,3)*RCM(n,3,2);
  RC(2,3)=RCM(o,2,0)*RCM(n,0,3)+RCM(o,2,1)*RCM(n,1,3)+
          RCM(o,2,2)*RCM(n,2,3)+RCM(o,2,3)*RCM(n,3,3);
  RC(3,0)=RCM(o,3,0)*RCM(n,0,0)+RCM(o,3,1)*RCM(n,1,0)+
          RCM(o,3,2)*RCM(n,2,0)+RCM(o,3,3)*RCM(n,3,0);
  RC(3,1)=RCM(o,3,0)*RCM(n,0,1)+RCM(o,3,1)*RCM(n,1,1)+
          RCM(o,3,2)*RCM(n,2,1)+RCM(o,3,3)*RCM(n,3,1);
  RC(3,2)=RCM(o,3,0)*RCM(n,0,2)+RCM(o,3,1)*RCM(n,1,2)+
          RCM(o,3,2)*RCM(n,2,2)+RCM(o,3,3)*RCM(n,3,2);
  RC(3,3)=RCM(o,3,0)*RCM(n,0,3)+RCM(o,3,1)*RCM(n,1,3)+
          RCM(o,3,2)*RCM(n,2,3)+RCM(o,3,3)*RCM(n,3,3);
}

void vector_multiply (float * v, float * m)
{
  float o[4];
  int i;
  for (i = 0; i < 4; i++) o[i] = v[i];  
  v[0]=RC(0,0)*o[0]+RC(0,1)*o[1]+RC(0,2)*o[2]+RC(0,3)*o[3];
  v[1]=RC(1,0)*o[0]+RC(1,1)*o[1]+RC(1,2)*o[2]+RC(1,3)*o[3];
  v[2]=RC(2,0)*o[0]+RC(2,1)*o[1]+RC(2,2)*o[2]+RC(2,3)*o[3];
  v[3]=RC(3,0)*o[0]+RC(3,1)*o[1]+RC(3,2)*o[2]+RC(3,3)*o[3];
}

void gl_check_error()
{
  switch (glGetError()) {
  case GL_NO_ERROR: return;
  case GL_INVALID_ENUM: fprintf (stderr, "OpenGL: invalid enum\n"); break;
  case GL_INVALID_VALUE: fprintf (stderr, "OpenGL: invalid value\n"); break;
  case GL_INVALID_OPERATION: fprintf (stderr, "OpenGL: invalid operation\n");
    break;
  case GL_INVALID_FRAMEBUFFER_OPERATION:
    fprintf (stderr, "OpenGL: invalid framebuffer operation\n"); break;
  case GL_OUT_OF_MEMORY:
    fprintf (stderr, "OpenGL: out of memory\n"); break;
  case GL_STACK_UNDERFLOW:
    fprintf (stderr, "OpenGL: stack underflow\n"); break;
  case GL_STACK_OVERFLOW:
    fprintf (stderr, "OpenGL: stack overflow\n"); break;
  }
  abort();
}

void gl_get_frustum (Frustum * f)
{
  GLint v[4];
  glGetIntegerv (GL_VIEWPORT, v);
  gl_check_error();
  f->width = v[2];
  glGetFloatv (GL_MODELVIEW_MATRIX, f->m);
  gl_check_error();
  glGetFloatv (GL_PROJECTION_MATRIX, f->p);
  gl_check_error();
  float p[16];
  int i;
  for (i = 0; i < 16; i++) p[i] = f->p[i];
  matrix_multiply (p, f->m);

  /* right */
  f->n[0][0] = p[3] - p[0];
  f->n[0][1] = p[7] - p[4];
  f->n[0][2] = p[11] - p[8];
  f->d[0]    = p[15] - p[12];
   
  /* left */
  f->n[1][0] = p[3] + p[0];
  f->n[1][1] = p[7] + p[4];
  f->n[1][2] = p[11] + p[8];
  f->d[1]    = p[15] + p[12];
  
  /* top */
  f->n[2][0] = p[3] - p[1];
  f->n[2][1] = p[7] - p[5];
  f->n[2][2] = p[11] - p[9];
  f->d[2]    = p[15] - p[13];

  /* bottom */
  f->n[3][0] = p[3] + p[1];
  f->n[3][1] = p[7] + p[5];
  f->n[3][2] = p[11] + p[9];
  f->d[3]    = p[15] + p[13];
  
  /* front */
  f->n[4][0] = p[3] + p[2];
  f->n[4][1] = p[7] + p[6];
  f->n[4][2] = p[11] + p[10];
  f->d[4]    = p[15] + p[14];
  
  /* back */
  f->n[5][0] = p[3] - p[2];
  f->n[5][1] = p[7] - p[6];
  f->n[5][2] = p[11] - p[10];
  f->d[5]    = p[15] - p[14];
  
  for (i = 0; i < 6; i++) {
    float n = sqrt(f->n[i][0]*f->n[i][0] +
		   f->n[i][1]*f->n[i][1] +
		   f->n[i][2]*f->n[i][2]);
    if (n > 0.) {
      f->n[i][0] /= n; f->n[i][1] /= n; f->n[i][2] /= n;
      f->d[i] /= n;
    }
  }
}

/*
  Returns 0 if the sphere is outside the view frustum, 1, if it is
  inside and -1 if it is partly inside. 
*/
int sphere_in_frustum (double x, double y, double z, double r, Frustum * f)
{
  int I1 = 0, i, I = 1;
  for (i = 0; i < 6; i++) {
    double d = f->n[i][0]*x + f->n[i][1]*y + f->n[i][2]*z + f->d[i];
    if (d < -r) {
      I = 0;
      break;
    }
    if (d < r)
      I = -1;
  }
  if (I == 1)
    return 1;
  if (I == -1)
    I1 = -1;
  return I1;
}

/*
  Returns the diameter (in pixels) of a sphere projected on the
  screen.
*/
float sphere_diameter (double x, double y, double z, double r, Frustum * f)
{
  float v[4];
  v[0] = x; v[1] = y; v[2] = z; v[3] = 1.;
  vector_multiply (v, f->m);
  v[0] = r;
  vector_multiply (v, f->p);
  float rp = v[3] == 0. ? 0 : v[0]*f->width/v[3];
  return rp;
}

static bool atobool (char * s)
{
  if (!strcmp (s, "true"))
    return true;
  if (!strcmp (s, "false"))
    return false;
  return atoi (s) != 0;
}
  
static bool args (Params * p, char * val)
{
  static char * name[] = { "string", "int", "unsigned",
			   "bool", "float", "double" };  
  switch (p->type) {

  case pstring:
    if (val[0] != '"') {
      fprintf (stderr, "expecting a string for '%s' got '%s'\n", p->key, val);
      return false;
    }
    if (val[strlen(val) - 1] != '"') {
      fprintf (stderr, "unterminated quoted string '%s'\n", val);
      return false;
    }
    val[strlen(val) - 1] = '\0';
    *((char **)p->val) = &val[1];
    break;

  case pint: case punsigned: case pbool: case pdouble: case pfloat:
    if (val[0] == '"') {
      fprintf (stderr, "expecting a %s for '%s' got %s\n",
	       name[p->type], p->key, val);
      return false;
    }
    if (!p->n) {
      switch (p->type) {
      case pint: *((int *)p->val) = atoi(val); break;
      case punsigned: *((unsigned *)p->val) = atoi(val); break;
      case pbool: *((bool *)p->val) = atobool(val); break;
      case pfloat: *((float *)p->val) = atof(val); break;
      case pdouble: *((double *)p->val) = atof(val); break;
      default: assert (false);
      }
    }
    else {
      if (val[0] != '{') {
	fprintf (stderr, "expecting an array for '%s' got %s\n", p->key, val);
	return false;
      }
      val++;
      int i = 0;
      char c = ',';
      while (i < p->n && c != '}') {
	char * s = strchr (val, ',');
	if (!s)
	  s = strchr (val, '}');
	if (!s) {
	  fprintf (stderr, "expecting an array for '%s' got %s\n", p->key, val);
	  return false;
	}
	c = *s;
	*s++ = '\0';
	switch (p->type) {
	case pint: ((int *)p->val)[i++] = atoi (val); break;
	case punsigned: ((unsigned *)p->val)[i++] = atoi (val); break;
	case pbool: ((bool *)p->val)[i++] = atobool (val); break;
	case pfloat: ((float *)p->val)[i++] = atof (val); break;
	case pdouble: ((double *)p->val)[i++] = atof (val); break;
	default: assert (false);
	}
	val = s;
      }
      if (c != '}') {
	fprintf (stderr, "expecting '}' for '%s' got %s\n", p->key, val);
	return false;
      }
    }
    break;

  default:
    assert (false);
  }
  return true;
}

static char * find_comma (char * s)
{
  int par = 0;
  while (*s != '\0') {
    if (*s == ',' && par == 0) {
      *s = '\0';
      return s + 1;
    }
    if (*s == '{')
      par++;
    else if (*s == '}')
      par--;
    s++;
  }
  return NULL;
}

int parse_params (Params * params)
{
  char * s;
  int i = 0, n = 0;
  Params * p = params;
  while (p->key) p++, n++;
  if (!(s = strtok (NULL, ");")) || s[0] == '\n')
    return false;
  while (s) {
    char * next = find_comma (s), * key = s;
    if ((s = strchr (key, '='))) {
      s[0] = '\0', s++;
      i = -1;
      Params * p = params;
      while (p->key && strcmp(p->key, key)) p++;
      if (!p->key) {
	fprintf (stderr, "unknown key '%s'\n", key);
	return false;
      }
      if (!args (p, s))
	return false;
    }
    else {
      if (i < 0) {
	fprintf (stderr, "anonymous value '%s' after keys\n", key);
	return false;
      }
      if (i >= n) {
	fprintf (stderr, "too many parameters: '%s' %d %d\n", key, i, n);
	return false;
      }
      if (!args (&params[i], key))
	return false;
      i++;
    }
    s = next;
  }
  return true;
}

/**
 * gl_feedback_begin:
 * @buffersize: the size of the feedback buffer.
 *
 * Returns: a newly setup openGL feedback buffer.
 */
float * gl_feedback_begin (unsigned buffersize)
{
  float * f = malloc (sizeof (float)*buffersize);
  glFeedbackBuffer (buffersize, GL_3D_COLOR, f);
  glRenderMode (GL_FEEDBACK);
  return f;
}

/**
 * gl_feedback_end:
 * @f: the #GlFeedback.
 * @sim: a #GfsSimulation.
 * @fp: a file pointer.
 * @format: the file format.
 *
 * Writes to @fp the contents of @f.
 *
 * Returns: %false if @f has overflowed, %true otherwise.
 */
bool gl_feedback_end (float * f, FILE * fp, enum FeedbackFormat format)
{
  GLint size = glRenderMode (GL_RENDER);
  if (size >= 0) {
    GLint used = size;
    GLfloat * current = f;
    GLdouble model[16], proj[16];
    GLint viewport[4];
    struct { double x, y, z; } p, lastp = { 1e30, 1e30, 1e30 };
    unsigned nps = 0, np = 0, linetoken = false;

    /* Header */
    switch (format) {
    case FEEDBACK_KML:
      fputs ("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
	     "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n"
	     "<Placemark>\n"
	     "<name>GfsView</name>\n"
	     "<MultiGeometry>\n", fp);
      break;
    default:
      break;
    }

    /* Body */
    glGetDoublev (GL_MODELVIEW_MATRIX, model);
    glGetDoublev (GL_PROJECTION_MATRIX, proj);
    glGetIntegerv (GL_VIEWPORT, viewport);
    while (used > 0)
      switch ((GLint) *current) {

      case GL_POINT_TOKEN :
	current++; used--;
	gluUnProject (current[0], current[1], current[2], model, proj, viewport,
		      &p.x, &p.y, &p.z);
	//	gfs_simulation_map_inverse (sim, &p);
	switch (format) {
	case FEEDBACK_GNU:
	  fprintf (fp, "%g %g %g\n\n", p.x, p.y, p.z); break;
	case FEEDBACK_OBJ:
	  fprintf (fp, "v %g %g %g\np -1\n", p.x, p.y, p.z); 
	  np++;
	  nps = 0;
	  break;
	case FEEDBACK_KML:
	  if (linetoken) {
	    fputs ("</coordinates>\n</LineString>\n", fp);
	    linetoken = false;
	  }
	  fprintf (fp, "<Point>\n<coordinates>%f,%f,%f</coordinates>\n"
		   "</Point>\n", p.x, p.y, p.z);
	  break;
	default:
	  assert (false);
	}
	current += 7; used -= 7;
	break;

      case GL_LINE_TOKEN :
      case GL_LINE_RESET_TOKEN :
	current++; used--;
	gluUnProject (current[0], current[1], current[2],
		      model, proj, viewport, &p.x, &p.y, &p.z);
	//	gfs_simulation_map_inverse (sim, &p);
	if (p.x != lastp.x || p.y != lastp.y || p.z != lastp.z) {
	  switch (format) {
	  case FEEDBACK_GNU:
	    fprintf (fp, "\n%g %g %g\n", p.x, p.y, p.z);
	    break;
	  case FEEDBACK_OBJ:
	    if (nps > 0) {
	      fputc ('l', fp);
	      int i;
	      for (i = nps; i <= np; i++)
		fprintf (fp, " %d", i);
	      fputc ('\n', fp);
	    }
	    nps = ++np;
	    fprintf (fp, "v %g %g %g\n", p.x, p.y, p.z);
	    break;
	  case FEEDBACK_KML:
	    if (linetoken)
	      fputs ("</coordinates>\n</LineString>\n", fp);
	    fprintf (fp, "<LineString>\n<coordinates>\n%f,%f,%f\n",
		     p.x, p.y, p.z);
	    linetoken = true;
	    break;
	  default:
	    assert (false);
	  }
	}
	current += 7; used -= 7;
	gluUnProject (current[0], current[1], current[2],
		      model, proj, viewport, &p.x, &p.y, &p.z);
	//	gfs_simulation_map_inverse (sim, &p);
	lastp = p;
	switch (format) {
	case FEEDBACK_GNU:
	  fprintf (fp, "%g %g %g\n", p.x, p.y, p.z); 
	  break;
	case FEEDBACK_OBJ:
	  fprintf (fp, "v %g %g %g\n", p.x, p.y, p.z); 
	  np++;
	  break;
	case FEEDBACK_KML:
	  fprintf (fp, "%f,%f,%f\n", p.x, p.y, p.z);
	  break;
	default:
	  assert (false);
	}
	current += 7; used -= 7;
	break;

      case GL_POLYGON_TOKEN : {
	GLint count = (GLint) current[1], vcount = 0;
	current += 2;
	used -= 2;
	if (format == FEEDBACK_KML) {
	  if (linetoken) {
	    fputs ("</coordinates>\n</LineString>\n", fp);
	    linetoken = false;
	  }
	  fputs ("<Polygon>\n<coordinates>\n", fp);
	}
	while (count > 0 && used > 0) {
	  gluUnProject (current[0], current[1], current[2],
			model, proj, viewport, 
			&p.x, &p.y, &p.z);
	  //	  gfs_simulation_map_inverse (sim, &p);
	  switch (format) {
	  case FEEDBACK_GNU:
	    fprintf (fp, "%g %g %g\n", p.x, p.y, p.z); break;
	  case FEEDBACK_OBJ:
	    fprintf (fp, "v %g %g %g\n", p.x, p.y, p.z); 
	    np++;
	    nps = 0;
	    break;
	  case FEEDBACK_KML:
	    fprintf (fp, "%f,%f,%f\n", p.x, p.y, p.z); break;
	  default:
	    assert (false);
	  }
	  current += 7; used -= 7;
	  count--;
	  vcount++;
	}
	switch (format) {
	case FEEDBACK_KML:
	  fputs ("</coordinates>\n</Polygon>\n", fp); break;
	case FEEDBACK_OBJ:
	  fprintf (fp, "f -%d", vcount--);
	  while (vcount)
	    fprintf (fp, " -%d", vcount--);
	  /* fall through */
	case FEEDBACK_GNU: fputc ('\n', fp); break;
	default:
	  assert (false);
	}
	break;
      }
      }

    /* Footer */
    switch (format) {
    case FEEDBACK_OBJ:
      if (np > nps && nps > 0) {
	fputc ('l', fp);
	int i;
	for (i = nps; i <= np; i++)
	  fprintf (fp, " %d", i);
	fputc ('\n', fp);
      }
      break;
     
    case FEEDBACK_KML:
      if (linetoken) {
	fputs ("</coordinates>\n</LineString>\n", fp);
	linetoken = false;
      }
      fputs ("</MultiGeometry>\n"
	     "</Placemark>\n"
	     "</kml>\n", fp);
      break;
    
    default:
      break;
    }
  }
  free (f);
  return (size >= 0);
}
