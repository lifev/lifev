#include "switch.hpp"

void set_switches_for_regionmesh( Switch & sw )
{
  sw.create("HAS_ALL_FACES");
  sw.create("HAS_ALL_EDGES");
  sw.create("HAS_BOUNDARY_FACES");
  sw.create("HAS_BOUNDARY_EDGES");
  sw.create("HAS_VOLUME_TO_FACES");
  sw.create("HAS_VOLUME_TO_EDGES");
  sw.create("HAS_BEEN_CHECKED");
  sw.create("FACES_HAVE_ADIACENCY");
}

