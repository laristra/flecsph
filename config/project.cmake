cinch_minimum_required(1.0)

project(tree)

cinch_add_application_directory("app/fluid")
#cinch_add_application_directory("app/flecsi_legion")
cinch_add_application_directory("app/sodtube")
#cinch_add_application_directory("app/bns")
cinch_add_application_directory("mpisph")

#cinch_add_library_target(mpisph "mpisph")

set(CINCH_HEADER_SUFFIXES "\\.h")
