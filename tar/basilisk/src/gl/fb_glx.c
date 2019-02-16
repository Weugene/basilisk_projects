#include "framebuffer.h"
#include "OffscreenContext.h"

void framebuffer_destroy (framebuffer * p)
{
  teardown_offscreen_context ((OffscreenContext *) p);
}

framebuffer * framebuffer_new (unsigned width, unsigned height)
{
  return (framebuffer *) create_offscreen_context (width, height);
}

unsigned char * framebuffer_image (framebuffer * p)
{
  return get_framebuffer_pixels ((OffscreenContext *) p);
}

fbdepth_t * framebuffer_depth (framebuffer * p)
{
  return get_framebuffer_depth ((OffscreenContext *) p);
}

