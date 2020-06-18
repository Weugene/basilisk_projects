/**
# Multiple "layers"

This adds support for multiple "layers" i.e. a constant number `nl` of
contiguous scalar fields. 

Note that codes using this header must be compiled with the
`-DLAYERS=1` compilation flag. */

int nl = 1, layer = 0;

/**
We redefine two types of block traversals. An "outer" traversal,
usable outside foreach() loops, which uses a global `layer` index... */

#undef foreach_block
@define foreach_block() for (layer = 0; layer < nl; layer++)
@define end_foreach_block() layer = 0

/**
... and an "inner" traversal, usable within foreach() loops, which
uses the `point.l` index as local layer index. */
  
@define foreach_block_inner() for (point.l = 0; point.l < nl; point.l++)
@define end_foreach_block_inner() point.l = 0

/**
The two indices are combined to access field values. In practice only
one of the indices is used. */
  
@undef _index
@def _index(a,m)
  (a.i + (layer + point.l + m < _attribute[a.i].block ? layer + point.l + m : 0))
@

@undef val
@define val(a,k,p,m) data(k,p,m)[_index(a,m)]

/**
foreach_layer() is just an alias for foreach_block(). */
  
#define foreach_layer() foreach_block()
