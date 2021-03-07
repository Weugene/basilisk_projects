typedef struct _framebuffer framebuffer;
typedef unsigned int fbdepth_t;

framebuffer * framebuffer_new (unsigned width, unsigned height);
void framebuffer_destroy (framebuffer * p);
unsigned char * framebuffer_image (framebuffer * p);
fbdepth_t * framebuffer_depth (framebuffer * p);
