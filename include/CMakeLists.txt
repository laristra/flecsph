###############################################
# sets include file sources
# this is _mostly_ so that IDEs (code &t) can
# locate header sources, but is also the more
# 'correct' way to define header-only libraries
###############################################

# NOTE: cmake > 3.13 required for relative path handling
# TODO: I don't love how every header is placed. An alternative
#       would be to put a CMakeLists.txt + target_sources in each
#       include subdirectory, but the the call-tree of
#       add_subdirectory gets iffy.
# TODO: refactor into individual interface targets and
#       propogate to up to eliminate unnecessary include
#       directives
target_sources(flecsph::flags
    INTERFACE
        tree.h
        lattice.h
        log.h
        params.h
        space_vector.h
        diagnostic.h
        tensor.h

        tree_topology/tree_geometry.h
        tree_topology/hashtable.h
        tree_topology/tree_utils.h
        tree_topology/filling_curve.h
        tree_topology/tree_types.h
        tree_topology/tree_topology.h

        physics/integration.h
        physics/eforce.h
        physics/kernels.h
        physics/fmm.h
        physics/boundary.h
        physics/node.h
        physics/viscosity.h
        physics/body.h
        physics/wvt.h
        physics/analysis.h
        physics/default_physics.h
        physics/density_profiles.h

        physics/eos/eos.h
)
